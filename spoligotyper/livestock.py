"""
Lineages of the livestock-associated M. tuberculosis complex, after Zwyer et al. 2021 (Open Research Europe 1:100,
https://doi.org/10.12688/openreseurope.14029.2): La1 (M. bovis), La2 (M. caprae), La3 (M. orygis), and the La1
sublineages La1.1 to La1.8.

livestock_snps.fasta has, for each SNP of livestock_barcode.tsv, the 61 bp of H37Rv centred on the SNP with the
ancestral and the derived allele. Seal counts the reads containing an exact 31-mer of each sequence (in the same pass
as the SNPs of the lineage barcode of Coll et al.). As in the KvarQ test suite published with the paper
(https://github.com/dbrites/LivestockAssociatedMTBC), a group is called when at least 2 of its SNPs carry the derived
allele; likewise, a group is mixed when at least 2 of its SNPs have both alleles.
"""

from . import snp_groups
from .snp_groups import GroupCall as LivestockCall  # noqa: F401 (the call of this barcode)
from .spoligotype import data_file

BARCODE = data_file('livestock_barcode.tsv')
SNPS_FASTA = data_file('livestock_snps.fasta')
MIN_SNPS = snp_groups.MIN_SNPS

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
    'La1.8': ('La1', 'La1.8', 'Eu1, unknown6, unknown7, unknown8'),
    'La1.8.1': ('La1.8', 'La1.8.1', 'Eu1'),
    'La1.8.2': ('La1.8', 'La1.8.2', 'unknown7'),
    'La1.8.X-unk6': ('La1.8', 'La1.8.X', 'unknown6'),
}

SCHEME = snp_groups.Scheme(BARCODE, SNPS_FASTA, tuple(GROUPS.items()), shared=(('La1_La2', ('La1', 'La2')),),
                           labels=(('La1_La2', 'La1/La2'), ('La1.2_BCG', 'La1.2 BCG')))
label, ancestors, root, on_one_path, mixed_summary = (SCHEME.label, SCHEME.ancestors, SCHEME.root,
                                                      SCHEME.on_one_path, SCHEME.mixed_summary)


def read_barcode(path=BARCODE):
    return snp_groups.read_barcode(str(path))


def ntm_conserved_snps(path=SNPS_FASTA):
    """SNPs ("<group>|<position>") with an allele whose sequence is also found in NTM genomes."""
    return snp_groups.ntm_conserved_snps(str(path))


def call_livestock(counts, file_type, barcode=None, contaminated=False):
    """La1 to La3 and the La1 sublineages: see snp_groups.Scheme.call."""
    return SCHEME.call(counts, file_type, contaminated, barcode)
