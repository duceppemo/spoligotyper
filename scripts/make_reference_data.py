#!/usr/bin/env python3
"""
Build the species and lineage reference data shipped in spoligotyper/data/, from public genomes.

markers.fasta
    100 bp chunks of the H37Rv genome (NC_000962.3), named after what they detect:
    MTBC_nn  control: found in all the MTBC genomes below and in none of the non-tuberculous mycobacteria (NTM).
             Their read depth measures how much MTBC DNA the sample contains.
    RD9_nn   region of difference 9: deleted in M. africanum and the animal-adapted lineages (incl. M. bovis).
    RD4_nn   RD4: deleted in M. bovis and BCG only.
    RD1_nn   RD1: deleted in all BCG strains.
    RD7_nn   RD7 (mce3 operon): deleted in lineage 6 (M. africanum) and the animal-adapted lineages.
    RD12_nn  RD12: deleted in M. bovis, BCG, M. caprae, M. orygis, and some M. canettii.
    The RD regions are the H37Rv segments missing from AF2122/97 (RD9, RD4, RD7, RD12) and BCG Pasteur (RD1). A chunk
    is kept only if it is found in every genome that has the region, and in none of the genomes lacking it or of the
    NTM.

lineage_snps.fasta
    For each SNP of lineage_barcode.tsv (Coll et al. 2014), the 61 bp of H37Rv centred on the SNP, with the
    reference and with the alternative allele: every 31-mer of these sequences contains the SNP. The description
    lists the NTM whose genome contains one of these 31-mers ("ntm=M. avium,M. marinum"): reads of such
    contaminants can add support to that allele.

Requirements: seal.sh (BBTools) and internet access to NCBI. Usage: python scripts/make_reference_data.py [cache_dir]
"""

import csv
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / 'spoligotyper' / 'data'
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={}'

H37RV = 'NC_000962.3'
AF2122 = 'NC_002945.4'  # M. bovis AF2122/97
BCG = 'NC_008769.1'  # M. bovis BCG Pasteur 1173P2
AFRICANUM = 'NC_015758.1'  # M. africanum GM041182, lineage 6
CANETTII = 'NC_015848.1'  # M. canettii CIPT 140010059
MTBC = [H37RV, AF2122, BCG, AFRICANUM, CANETTII]
NTM_NAMES = {'NC_010612.1': 'M. marinum',  # Strain M
             'NC_022663.1': 'M. kansasii',  # ATCC 12478
             'NC_008595.1': 'M. avium',  # 104
             'NC_010397.1': 'M. abscessus',  # ATCC 19977
             'NC_008611.1': 'M. ulcerans',  # Agy99
             'NC_008596.1': 'M. smegmatis'}  # MC2 155
NTM = list(NTM_NAMES)
# Region: (genomes lacking it, H37Rv search window, tiling step). The first genome defines the region's extent.
REGIONS = {'RD9': ([AF2122, BCG, AFRICANUM], (2326000, 2336000), 50),
           'RD4': ([AF2122, BCG], (1693000, 1712000), 100),
           'RD1': ([BCG], (4345000, 4362000), 100),
           'RD7': ([AF2122, BCG, AFRICANUM], (2206000, 2224000), 100),
           'RD12': ([AF2122, BCG, CANETTII], (3483000, 3490000), 50)}
CHUNK = 100
N_CONTROL, N_REGION = 40, 20


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))


def genome(accession, cache):
    path = cache / '{}.fasta'.format(accession)
    if not path.exists():
        print('Downloading', accession, file=sys.stderr)
        with urllib.request.urlopen(EFETCH.format(accession)) as response:
            path.write_bytes(response.read())
        time.sleep(0.5)  # NCBI allows 3 requests per second without an API key
    return ''.join(line.strip() for line in path.read_text().splitlines() if not line.startswith('>')).upper()


def contains(seq, genome_seq):
    return seq in genome_seq or reverse_complement(seq) in genome_seq


def deleted_segment(reference, other, window):
    """Longest run of the window, in 50 bp steps, whose 40-mers are missing from the other genome."""
    runs, start, previous = [], None, None
    for p in range(window[0], window[1], 50):
        if contains(reference[p:p + 40], other):
            continue
        if start is None or p - previous > 50:
            if start is not None:
                runs.append((start, previous + 40))
            start = p
        previous = p
    if start is not None:
        runs.append((start, previous + 40))
    return max(runs, key=lambda r: r[1] - r[0])


def seal_hits(chunks, genomes, cache, k=25, hdist=1):
    """{chunk name: {accession: found}} with the parameters used by spoligotyper for spacers and markers."""
    hits = {name: {} for name in chunks}
    with tempfile.TemporaryDirectory() as tmp:
        ref = Path(tmp) / 'chunks.fasta'
        ref.write_text(''.join('>{}\n{}\n'.format(name, seq) for name, seq in chunks.items()))
        for accession in genomes:
            stats = Path(tmp) / 'stats.tsv'
            subprocess.run(['seal.sh', '-Xmx1g', 'in={}'.format(cache / '{}.fasta'.format(accession)),
                            'ref={}'.format(ref), 'k={}'.format(k), 'hdist={}'.format(hdist), 'rcomp=t',
                            'maskmiddle=f', 'clearzone=999999', 'ambiguous=all', 'nzo=f', 'ow=t',
                            'stats={}'.format(stats), 'threads=4'],
                           check=True, capture_output=True)
            for line in stats.read_text().splitlines():
                if not line.startswith('#'):
                    name, count = line.split('\t')[:2]
                    hits[name.split()[0]][accession] = int(count) > 0
    return hits


def spread(names, n):
    step = max(1.0, len(names) / n)
    return [names[int(i * step)] for i in range(min(n, len(names)))]


def make_markers(cache):
    ref = genome(H37RV, cache)
    genomes = {a: genome(a, cache) for a in MTBC + NTM}
    chunks = {}
    # Control candidates: unique 100 bp segments spread over the genome (no IS6110, PE/PPE or other repeats)
    for p in range(20000, len(ref) - 20000, 14000):
        seq = ref[p:p + CHUNK]
        halves = (seq[:50], seq[50:])
        if all(ref.count(h) + ref.count(reverse_complement(h)) == 1 for h in halves):
            chunks['MTBC:{}-{}'.format(p + 1, p + CHUNK)] = seq
    for region, (lacking, window, step) in REGIONS.items():
        start, end = deleted_segment(ref, genomes[lacking[0]], window)
        print('{}: H37Rv {}-{} ({} bp)'.format(region, start + 1, end, end - start), file=sys.stderr)
        for p in range(start, end - CHUNK + 1, step):
            chunks['{}:{}-{}'.format(region, p + 1, p + CHUNK)] = ref[p:p + CHUNK]

    hits = seal_hits(chunks, MTBC + NTM, cache)
    selected = []
    for region, lacking in [('MTBC', [])] + [(r, v[0]) for r, v in REGIONS.items()]:
        names = [n for n in chunks if n.startswith(region + ':')
                 and all(hits[n][a] for a in MTBC if a not in lacking)
                 and not any(hits[n][a] for a in lacking + NTM)]
        kept, last = [], None  # No overlapping chunks
        for name in names:
            start = int(name.split(':')[1].split('-')[0])
            if last is None or start - last >= CHUNK:
                kept.append(name)
                last = start
        kept = spread(kept, N_CONTROL if region == 'MTBC' else N_REGION)
        print('{}: {} chunks pass, {} kept'.format(region, len(names), len(kept)), file=sys.stderr)
        selected += [('{}_{:02d}'.format(region, i), name.split(':')[1], chunks[name])
                     for i, name in enumerate(kept, 1)]
    with open(DATA / 'markers.fasta', 'w') as f:
        for name, coordinates, seq in selected:
            f.write('>{} H37Rv:{}\n{}\n'.format(name, coordinates, seq))


def make_lineage_snps(cache):
    ref = genome(H37RV, cache)
    with open(DATA / 'lineage_barcode.tsv') as f:
        rows = list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))
    sequences = {}
    for row in rows:
        p = int(row['position']) - 1
        window = ref[p - 30:p + 31]
        if window[30] != row['ref']:
            raise SystemExit('Reference allele mismatch for lineage {}'.format(row['lineage']))
        for allele in ('ref', 'alt'):
            sequences['{}|{}|{}'.format(row['lineage'], row['position'], allele)] = \
                window[:30] + row[allele] + window[31:]
    for accession in NTM:
        genome(accession, cache)
    hits = seal_hits(sequences, NTM, cache, k=31, hdist=0)
    with open(DATA / 'lineage_snps.fasta', 'w') as f:
        for name, seq in sequences.items():
            ntm = sorted(NTM_NAMES[a] for a in NTM if hits[name][a])
            f.write('>{}{}\n{}\n'.format(name, ' ntm={}'.format(','.join(ntm)) if ntm else '', seq))
            if ntm:
                print('{} also in {}'.format(name, ', '.join(ntm)), file=sys.stderr)
    print('{} lineage SNPs'.format(len(rows)), file=sys.stderr)


if __name__ == '__main__':
    cache_dir = Path(sys.argv[1] if len(sys.argv) > 1 else 'genomes')
    cache_dir.mkdir(parents=True, exist_ok=True)
    make_markers(cache_dir)
    make_lineage_snps(cache_dir)
