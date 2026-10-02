"""
Lineages of the livestock-associated M. tuberculosis complex, after Zwyer et al. 2021 (Open Research Europe 1:100,
https://doi.org/10.12688/openreseurope.14029.2): La1 (M. bovis), La2 (M. caprae), La3 (M. orygis), and the La1
sublineages La1.1 to La1.8.

livestock_snps.fasta has, for each SNP of livestock_barcode.tsv, the 61 bp of H37Rv centred on the SNP with the
ancestral and the derived allele. Seal counts the reads containing an exact 31-mer of each sequence (in the same pass
as the SNPs of the lineage barcode of Coll et al.). As in the KvarQ test suite of the paper, a group is called when
at least 2 of its SNPs carry the derived allele.
"""

import csv
from dataclasses import dataclass, field

from .lineage import CALL_FRACTION, MIN_READS, MIXED_FRACTION, MIXED_MIN_READS, SnpCall, confirmed_mixed
from .spoligotype import data_file

BARCODE = data_file('livestock_barcode.tsv')
SNPS_FASTA = data_file('livestock_snps.fasta')
MIN_SNPS = 2  # SNPs with the derived allele to call a group (Zwyer et al. 2021)

# Group: (parent, lineage reported, former name). La1_La2 (shared by M. bovis and M. caprae) is not reported: it
# confirms La1 or La2. The groups of La1.7.X and La1.8.X are not monophyletic: they keep the names of their groups.
GROUPS = {
    'La1_La2': (None, '', ''),
    'La1': (None, 'La1', 'M. bovis'),
    'La2': (None, 'La2', 'M. caprae'),
    'La3': (None, 'La3', 'M. orygis'),
    'La1.1': ('La1', 'La1.1', 'pyrazinamide-susceptible M. bovis'),
    'La1.2': ('La1', 'La1.2', 'Eu3, unknown2'),
    'La1.2_BCG': ('La1.2', 'La1.2', 'BCG'),
    'La1.3': ('La1', 'La1.3', 'Af2'),
    'La1.4': ('La1', 'La1.4', 'unknown3'),
    'La1.5': ('La1', 'La1.5', 'unknown9'),
    'La1.6': ('La1', 'La1.6', 'Af1'),
    'La1.7': ('La1', 'La1.7', 'Eu2, unknown4, unknown5'),
    'La1.7.1': ('La1.7', 'La1.7.1', 'Eu2'),
    'La1.7.X-unk4': ('La1.7', 'La1.7.X', 'unknown4'),
    'La1.7.X-unk5': ('La1.7', 'La1.7.X', 'unknown5'),
    'La1.8': ('La1', 'La1.8', 'Eu1, unknown6, unknown7'),
    'La1.8.1': ('La1.8', 'La1.8.1', 'Eu1'),
    'La1.8.2': ('La1.8', 'La1.8.2', 'unknown7'),
    'La1.8.X-unk6': ('La1.8', 'La1.8.X', 'unknown6'),
}


@dataclass
class LivestockCall:
    lineage: str = ''  # Most specific group, e.g. "La1.8.1", or "" when none; "mixed: ..." for a mix
    name: str = ''  # Former name of the group, e.g. "Eu1"
    group: str = ''  # Group of the barcode, e.g. "La1.7.X-unk4" for lineage "La1.7.X"
    called: list = field(default_factory=list)  # Groups with at least 2 SNPs with the derived allele
    snps: list = field(default_factory=list)  # SnpCall for every SNP covered by reads (lineage = group)
    mixed: list = field(default_factory=list)  # SnpCall with both alleles
    conflict: bool = False  # Groups that cannot occur together, e.g. La1 and La3, or La1.7 and La1.8

    @property
    def main(self):
        """La1, La2 or La3, or '' when none or conflicting."""
        roots = {root(g) for g in self.called} - {''}
        return roots.pop() if len(roots) == 1 and not self.conflict else ''


def read_barcode(path=BARCODE):
    with open(path) as f:
        return list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))


def ntm_conserved_snps(path=SNPS_FASTA):
    """SNPs ("<group>|<position>") with an allele whose sequence is also found in NTM genomes."""
    with open(path) as f:
        return {line[1:].split(' ')[0].rsplit('|', 1)[0] for line in f if line.startswith('>') and ' ntm=' in line}


def ancestors(group):
    """The group and its parents, e.g. La1.8.1, La1.8, La1."""
    chain = []
    while group:
        chain.append(group)
        group = GROUPS[group][0]
    return chain


def root(group):
    return '' if group == 'La1_La2' else ancestors(group)[-1]


def on_one_path(groups):
    """True when the groups can occur in one strain: each is an ancestor of the most specific one (La1_La2 goes with
    La1 or La2)."""
    groups = set(groups)
    if 'La1_La2' in groups:
        groups.discard('La1_La2')
        if groups and {root(g) for g in groups} - {'La1', 'La2'}:
            return False
    if not groups:
        return True
    deepest = max(groups, key=lambda g: len(ancestors(g)))
    return groups <= set(ancestors(deepest))


def mixed_summary(mixed):
    """e.g. "La1 35%, La1.8 36%": the derived allele fraction of each group, over its SNPs with both alleles."""
    pooled = {}
    for s in mixed:
        pooled.setdefault(GROUPS[s.lineage][1] or s.lineage.replace('_', '/'), []).append(s)
    return ', '.join('{} {:.0f}%'.format(name, 100 * sum(s.lineage_reads for s in snps) / sum(s.reads for s in snps))
                     for name, snps in pooled.items())


def call_livestock(counts, file_type, barcode=None, contaminated=False):
    """
    :param counts: {"<group>|<position>|ancestral" or "...|derived": reads} from Seal
    :param contaminated: the sample contains non-MTBC DNA: SNPs conserved in NTM are not used
    :return: LivestockCall
    """
    barcode = barcode or read_barcode()
    conserved_snps = ntm_conserved_snps()
    minimum = 1 if file_type == 'fasta' else MIN_READS
    result = LivestockCall()
    positives = {}
    for row in barcode:
        key = '{}|{}'.format(row['group'], row['position'])
        snp = SnpCall(row['group'], int(row['position']), row['gene'], counts.get(key + '|derived', 0),
                      counts.get(key + '|ancestral', 0))
        if snp.reads == 0:
            continue
        result.snps.append(snp)
        if key in conserved_snps and contaminated:
            continue
        if snp.reads >= minimum and snp.fraction >= CALL_FRACTION:
            positives[snp.lineage] = positives.get(snp.lineage, 0) + 1
        if file_type == 'fastq' and key not in conserved_snps and snp.reads >= MIXED_MIN_READS and \
                MIXED_FRACTION[0] <= snp.fraction <= MIXED_FRACTION[1]:
            result.mixed.append(snp)
    result.called = [group for group in GROUPS if positives.get(group, 0) >= MIN_SNPS]
    result.mixed = confirmed_mixed(result.mixed, result.snps, lambda group: ancestors(group)[1:])
    if result.mixed:
        result.lineage = 'mixed: ' + mixed_summary(result.mixed)
    reported = [g for g in result.called if g != 'La1_La2']
    if not reported:
        return result
    result.conflict = not on_one_path(result.called)
    if result.mixed:
        pass
    elif result.conflict:
        result.lineage = 'mixed: ' + ', '.join(sorted({GROUPS[g][1] for g in reported}))
    else:
        result.group = max(reported, key=lambda g: len(ancestors(g)))
        _, result.lineage, result.name = GROUPS[result.group]
    return result
