"""Synthetic genomes and reads with a known spoligotype, species and lineage."""

import gzip
import random

import pytest

from spoligotyper.lineage import LINEAGE_SNPS_FASTA
from spoligotyper.species import MARKERS_FASTA
from spoligotyper.spoligotype import SPACERS_FASTA

DR = 'GTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGAC'  # 36 bp direct repeat between spacers
SB0140 = '1101101000001110111111111111111111111100000'  # M. bovis AF2122/97
SB0120 = '1101111101111110111111111111111111111100000'  # M. bovis BCG
H37RV = '1111111111111111111001111111111100001111111'  # M. tuberculosis H37Rv, not in the M. bovis database

# Species profiles: regions of difference present, and lineage SNPs carrying the alternative allele.
# Lineages 4 and 4.9 are defined by the H37Rv (reference) allele, so other lineages carry the alternative one.
H37RV_PROFILE = (('RD9', 'RD4', 'RD1'), ())
BOVIS_PROFILE = (('RD1',), ('BOV', 'BOV_AFRI', '4', '4.9'))
BCG_PROFILE = ((), ('BOV', 'BOV_AFRI', '4', '4.9'))


def read_fasta(path):
    records, name = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                records[name] = ''
            else:
                records[name] += line.strip()
    return records


def random_seq(rng, length):
    return ''.join(rng.choice('ACGT') for _ in range(length))


def genome(binary, rng, profile=None, flank=3000):
    """
    Random sequence with a DR locus containing the spacers set to 1 in the binary pattern, and, with a profile,
    the MTBC control regions, the regions of difference present and one allele of each lineage SNP.
    """
    spacers = list(read_fasta(SPACERS_FASTA).values())
    parts = [random_seq(rng, flank), ''.join(DR + s for s, bit in zip(spacers, binary, strict=True) if bit == '1') + DR]
    if profile:
        regions, alt_lineages = profile
        for name, seq in read_fasta(MARKERS_FASTA).items():
            if name.startswith('MTBC_') or name.split('_')[0] in regions:
                parts += [random_seq(rng, 50), seq]
        for name, seq in read_fasta(LINEAGE_SNPS_FASTA).items():
            lineage, _, allele = name.split('|')
            if allele == ('alt' if lineage in alt_lineages else 'ref'):
                parts += [random_seq(rng, 50), seq]
    parts.append(random_seq(rng, flank))
    return ''.join(parts)


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def reads(seq, rng, depth, length=100):
    """Random reads over seq, on both strands."""
    for _ in range(int(depth * len(seq) / length)):
        p = rng.randint(0, len(seq) - length)
        read = seq[p:p + length]
        yield reverse_complement(read) if rng.random() < 0.5 else read


def write_reads(path, seq, rng, depth=20, length=100, mate=None, extra=None):
    """Gzipped fastq. extra: another sequence whose reads are added (mixed or contaminated samples)."""
    with gzip.open(path, 'wt') as f:
        all_reads = list(reads(seq, rng, depth, length))
        if extra:
            all_reads += list(reads(extra[0], rng, extra[1], length))
        for i, read in enumerate(all_reads):
            if mate == 2:
                read = reverse_complement(read)
            f.write('@r{}/{}\n{}\n+\n{}\n'.format(i, mate or 1, read, 'I' * length))


@pytest.fixture(scope='session')
def data(tmp_path_factory):
    rng = random.Random(1)
    folder = tmp_path_factory.mktemp('data')
    bovis = genome(SB0140, rng, BOVIS_PROFILE)
    h37rv = genome(H37RV, rng, H37RV_PROFILE)
    (folder / 'AF2122.fasta').write_text('>chr\n{}\n'.format(bovis))
    (folder / 'H37Rv.fna').write_text('>chr\n{}\n'.format(h37rv))
    (folder / 'BCG.fasta').write_text('>chr\n{}\n'.format(genome(SB0120, rng, BCG_PROFILE)))
    write_reads(folder / 'bovis_R1.fastq.gz', bovis, rng, depth=10, mate=1)
    write_reads(folder / 'bovis_R2.fastq.gz', bovis, rng, depth=10, mate=2)
    write_reads(folder / 'bovis_single.fq.gz', bovis, rng)
    write_reads(folder / 'low_coverage.fastq.gz', bovis, rng, depth=10)
    write_reads(folder / 'mixed.fastq.gz', h37rv, rng, depth=30, extra=(bovis, 15))
    (folder / 'not_mtbc.fasta').write_text('>chr\n{}\n'.format(random_seq(rng, 5000)))
    return folder
