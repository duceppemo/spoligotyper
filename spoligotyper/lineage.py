"""
Lineage from the SNP barcode of Coll et al. 2014 (Nat Commun 5:4812, https://doi.org/10.1038/ncomms5812).

lineage_snps.fasta has, for each of the 62 SNPs, the 61 bp of H37Rv centred on the SNP with each allele. Seal counts
the reads containing an exact 31-mer of each sequence, and every such 31-mer contains the SNP: the counts are the
reads supporting each allele.

A few of these sequences are conserved in non-tuberculous mycobacteria (NTM), whose reads can then add support to
one allele (see "ntm=" in the fasta descriptions). These SNPs are not used to detect mixed samples, and are ignored
when the sample looks contaminated.
"""

import csv
from dataclasses import dataclass, field

from .spoligotype import data_file

BARCODE = data_file('lineage_barcode.tsv')
LINEAGE_SNPS_FASTA = data_file('lineage_snps.fasta')
KMER_SIZE = 31
MIN_READS = 3  # Reads covering a SNP to call it, for fastq files
CALL_FRACTION = 0.8  # Fraction of the reads with the lineage allele to call the lineage
MIXED_MIN_READS = 10
MIXED_FRACTION = (0.15, 0.85)  # Lineage allele fraction suggesting a mix of lineages
CONTAMINATED_FRACTION = 0.8  # Below this fraction of MTBC reads, SNPs conserved in NTM are not used


@dataclass
class SnpCall:
    lineage: str
    position: int
    locus: str
    lineage_reads: int
    other_reads: int

    @property
    def reads(self):
        return self.lineage_reads + self.other_reads

    @property
    def fraction(self):
        return self.lineage_reads / self.reads if self.reads else 0.0


@dataclass
class LineageCall:
    lineage: str = ''  # Most specific lineage, e.g. "4.3.4.2", "BOV", or "" when none
    name: str = ''  # e.g. "Euro-American (LAM)"
    spoligotypes: str = ''  # Main spoligotype families of the lineage in Coll et al., e.g. "LAM"
    called: list = field(default_factory=list)  # All the lineages whose SNP is present
    snps: list = field(default_factory=list)  # SnpCall for every SNP covered by reads
    mixed: list = field(default_factory=list)  # SnpCall with both alleles
    conflict: bool = False  # Lineages that cannot occur together, e.g. 2 and 4


def ntm_conserved(path=LINEAGE_SNPS_FASTA):
    """Lineages whose SNP sequence (either allele) is also found in NTM genomes."""
    with open(path) as f:
        return {line[1:].split('|')[0] for line in f if line.startswith('>') and ' ntm=' in line}


def read_barcode(path=BARCODE):
    with open(path) as f:
        return list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))


def on_one_path(lineages):
    """
    True when the lineages can occur in one strain: nested human lineages (4, 4.3, 4.3.4), or one of the other
    lineages (5, 6, 7, BOV). BOV_AFRI is the clade of M. bovis and lineage 6: it only goes with BOV or 6.
    """
    lineages = set(lineages)
    human = sorted((lin for lin in lineages if lin[0] in '1234'), key=len)
    if any(not human[i + 1].startswith(human[i] + '.') for i in range(len(human) - 1)):
        return False
    groups = {'human' if lin[0] in '1234' else lin for lin in lineages}
    if 'BOV_AFRI' in groups:
        groups.discard('BOV_AFRI')
        return groups <= {'BOV'} or groups <= {'6'}
    return len(groups) <= 1


def call_lineage(counts, file_type, barcode=None, contaminated=False):
    """
    :param counts: {"<lineage>|<position>|ref" or "...|alt": reads} from Seal
    :param contaminated: the sample contains non-MTBC DNA: SNPs conserved in NTM are not used
    :return: LineageCall
    """
    barcode = barcode or read_barcode()
    conserved = ntm_conserved()
    minimum = 1 if file_type == 'fasta' else MIN_READS
    result = LineageCall()
    for row in barcode:
        key = '{}|{}|'.format(row['lineage'], row['position'])
        ref, alt = counts.get(key + 'ref', 0), counts.get(key + 'alt', 0)
        lineage_reads, other_reads = (alt, ref) if row['lineage_allele'] == 'alt' else (ref, alt)
        snp = SnpCall(row['lineage'], int(row['position']), row['locus'], lineage_reads, other_reads)
        if snp.reads == 0:
            continue
        result.snps.append(snp)
        if snp.lineage in conserved and contaminated:
            continue
        if snp.reads >= minimum and snp.fraction >= CALL_FRACTION:
            result.called.append(snp.lineage)
        if file_type == 'fastq' and snp.lineage not in conserved and snp.reads >= MIXED_MIN_READS and \
                MIXED_FRACTION[0] <= snp.fraction <= MIXED_FRACTION[1]:
            result.mixed.append(snp)
    if result.mixed:  # e.g. "mixed: 4.9 58%, BOV 31%": the lineage allele fraction of each mixed SNP
        result.lineage = 'mixed: ' + ', '.join('{} {:.0f}%'.format(s.lineage, s.fraction * 100)
                                               for s in sorted(result.mixed, key=lambda s: -s.fraction))
    if not result.called:
        return result
    result.conflict = not on_one_path(result.called)
    names = {row['lineage']: row for row in barcode}
    if 'BOV' in result.called:
        best = 'BOV'
    elif set(result.called) & {'5', '6', '7'}:
        best = sorted(set(result.called) & {'5', '6', '7'})[0]
    else:
        human = [lin for lin in result.called if lin[0] in '1234']
        best = max(human, key=lambda lin: lin.count('.')) if human else result.called[0]
    if result.mixed:
        pass
    elif result.conflict:
        result.lineage = 'mixed: ' + ', '.join(sorted(result.called))
    else:
        result.lineage = best
        result.name = names[best]['name']
        result.spoligotypes = names[best]['spoligotypes']
    return result
