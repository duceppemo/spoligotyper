"""
Sublineages of M. tuberculosis lineage 1, after Netikul et al. 2022 (Sci Rep 12:1565,
https://doi.org/10.1038/s41598-022-05524-0): L1.1 to L1.3 and their sublineages down to the fourth level (e.g.
L1.1.1.10, L1.2.2.3), from the 1,835 sublineage-specific SNPs of the paper (Supplementary Table S6). Counted in the
same Seal pass as the lineage SNPs, and called like the livestock lineages (see snp_groups).

The paper uses the revised lineage 1 nomenclature (L1.1 to L1.3), which differs from the lineage 1 sublineages of the
barcode of Coll et al. 2014: Coll's 1.2.1 is L1.2.2, and Coll's 1.2.2 is L1.3. The names are kept apart: Coll's in the
Lineage column, these in the L1Sublineage column.
"""

from . import snp_groups
from .spoligotype import data_file

BARCODE = data_file('l1_barcode.tsv')
SNPS_FASTA = data_file('l1_snps.fasta')

# Name of each sublineage: the equivalent of Coll et al. 2014, and the spoligotype families of at least 20% of the
# isolates of the sublineage in the paper (Supplementary Table S4)
NAMES = {
    'L1.1': 'Coll 1.1',
    'L1.1.1': 'Coll 1.1.1',
    'L1.1.1.1': 'typical spoligotypes: EAI4-VNM',
    'L1.1.1.2': 'typical spoligotypes: EAI5',
    'L1.1.1.3': 'typical spoligotypes: EAI5',
    'L1.1.1.4': 'typical spoligotypes: EAI5',
    'L1.1.1.5': 'typical spoligotypes: EAI5, EAI1-SOM',
    'L1.1.1.6': 'typical spoligotypes: EAI5',
    'L1.1.1.7': 'typical spoligotypes: EAI5',
    'L1.1.1.8': 'typical spoligotypes: EAI5',
    'L1.1.1.9': 'typical spoligotypes: EAI1-SOM, EAI5',
    'L1.1.1.10': 'typical spoligotypes: EAI1-SOM',
    'L1.1.1.11': 'typical spoligotypes: EAI5',
    'L1.1.2': 'Coll 1.1.2',
    'L1.1.2.1': 'typical spoligotypes: EAI5',
    'L1.1.2.2': 'typical spoligotypes: EAI3-IND, EAI5',
    'L1.1.3': 'Coll 1.1.3',
    'L1.1.3.1': 'typical spoligotypes: EAI6-BGD1',
    'L1.1.3.2': 'typical spoligotypes: EAI6-BGD1',
    'L1.1.3.3': 'typical spoligotypes: EAI6-BGD1',
    'L1.1.3.4': 'typical spoligotypes: EAI6-BGD1',
    'L1.2': '',
    'L1.2.1': 'typical spoligotypes: EAI6-BGD1',
    'L1.2.2': 'Coll 1.2.1',
    'L1.2.2.1': 'typical spoligotypes: EAI2-Manila',
    'L1.2.2.2': 'typical spoligotypes: EAI2-nonthaburi',
    'L1.2.2.3': 'typical spoligotypes: EAI2-Manila',
    'L1.2.2.4': 'typical spoligotypes: EAI2-Manila, EAI5',
    'L1.2.2.5': 'typical spoligotypes: EAI2-Manila',
    'L1.3': 'Coll 1.2.2',
    'L1.3.1': 'typical spoligotypes: EAI1-SOM',
    'L1.3.2': 'typical spoligotypes: EAI1-SOM',
}


def parent(group):
    """The closest group of the barcode containing this one, e.g. L1.1.1 for L1.1.1.10, or None."""
    parts = group.split('.')
    while len(parts) > 2:
        parts.pop()
        if '.'.join(parts) in NAMES:
            return '.'.join(parts)
    return None


# Groups have 4 to 224 SNPs: at least half of the covered SNPs of a group must carry its allele (M. canettii shares a
# few of them)
SCHEME = snp_groups.Scheme(BARCODE, SNPS_FASTA, tuple((g, (parent(g), g, name)) for g, name in NAMES.items()),
                           min_fraction=0.5)


def call_l1(counts, file_type, contaminated=False):
    """Lineage 1 sublineages: see snp_groups.Scheme.call."""
    return SCHEME.call(counts, file_type, contaminated)
