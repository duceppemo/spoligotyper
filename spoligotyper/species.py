"""
Species check from regions of difference (RD), and amount of MTBC DNA in the sample.

markers.fasta holds 100 bp chunks of the H37Rv genome (see scripts/make_reference_data.py):
  MTBC  control chunks, found in all MTBC genomes and in no non-tuberculous mycobacteria
  RD1   deleted in BCG and the Dassie bacillus
  RD4   deleted in M. bovis and BCG (and some M. canettii)
  RD7   deleted in lineage 6 (M. africanum) and the animal-adapted lineages
  RD9   deleted in M. africanum (lineages 5 and 6) and the animal-adapted lineages, including M. bovis
  RD12  deleted in M. bovis, BCG, M. caprae and M. orygis (and some M. canettii)

Each region's read depth is compared with the control depth: about the same when the region is present, 0 when it
is deleted, and in between for a mix of strains with and without it.
"""

import statistics
from dataclasses import dataclass, field

from .spoligotype import data_file

MARKERS_FASTA = data_file('markers.fasta')
CONTROL = 'MTBC'
REGIONS = ('RD1', 'RD4', 'RD7', 'RD9', 'RD12')
PRESENT, DELETED, PARTIAL = 'present', 'deleted', 'partial'
PRESENT_RATIO, DELETED_RATIO = 0.5, 0.1  # Region depth / control depth
MIN_CONTROL_READS = 3  # Median reads per control chunk to call the species from reads
MIN_CONTROL_FRACTION = 0.5  # Fraction of the control chunks found
CHUNK = 100


@dataclass
class SpeciesCheck:
    mtbc: bool = False  # Enough MTBC DNA to call the regions
    control_depth: float = 0.0  # Median reads (or contigs) per control chunk
    control_found: float = 0.0  # Fraction of the control chunks found
    regions: dict = field(default_factory=dict)  # {region: (state, depth ratio)}
    mtbc_fraction: float | None = None  # Estimated fraction of the reads from MTBC, reads only
    species: str = ''

    def state(self, region):
        return self.regions.get(region, ('', 0))[0]

    def summary(self):
        """e.g. "RD9 deleted, RD4 deleted, RD1 present"."""
        return ', '.join('{} {}'.format(r, self.state(r)) for r in REGIONS if r in self.regions)


def chunk_counts(counts, prefix):
    return [c for name, c in counts.items() if name.startswith(prefix + '_')]


def expected_control_reads(depth, read_length, paired):
    """
    Reads expected on a 100 bp control chunk for a pure MTBC sample: reads overlapping the chunk by 25 bp or more.
    Seal counts both reads of a pair when either contains the chunk, hence twice as many for paired-end reads.
    """
    return depth * (read_length + CHUNK - 2 * 25 + 1) / read_length * (2 if paired else 1)


def check_species(counts, file_type, depth=None, read_length=None, paired=False, lineages=(), mixed=False):
    """
    :param counts: {marker name: count} from Seal
    :param depth: estimated sequencing depth (all reads), for the MTBC fraction
    :param lineages: lineages called by the SNP barcode, to refine the species (e.g. M. africanum)
    :param mixed: the lineage SNPs indicate a mix of strains
    """
    check = SpeciesCheck()
    control = chunk_counts(counts, CONTROL)
    if not control:
        return check
    check.control_depth = statistics.median(control)
    check.control_found = sum(c > 0 for c in control) / len(control)
    if depth and read_length:
        expected = expected_control_reads(depth, read_length, paired)
        check.mtbc_fraction = min(1.0, check.control_depth / expected) if expected else None
    minimum = 1 if file_type == 'fasta' else MIN_CONTROL_READS
    check.mtbc = check.control_found >= MIN_CONTROL_FRACTION and check.control_depth >= minimum
    if not check.mtbc:
        check.species = 'MTBC not detected'
        return check

    for region in REGIONS:
        region_counts = chunk_counts(counts, region)
        if not region_counts:
            continue
        ratio = statistics.median(region_counts) / check.control_depth
        state = PRESENT if ratio >= PRESENT_RATIO else DELETED if ratio <= DELETED_RATIO else PARTIAL
        check.regions[region] = (state, ratio)
    check.species = 'MTBC, mixed sample?' if mixed else call_species(check, lineages)
    return check


def name_species(check, lineages=(), mixed=False, spacers=True):
    """Name the species once the lineage is known (check_species runs before the lineage call)."""
    if check.mtbc:
        check.species = 'MTBC, mixed sample?' if mixed else call_species(check, lineages, spacers)


# Species from the RD profile (RD1, RD4, RD7, RD9, RD12; + present, - deleted), after the classical RD PCR scheme.
# Lineage 5 (M. africanum West African 1) keeps RD7; lineage 6 (West African 2) lost it, like the animal lineages.
RD_PROFILES = {
    '+++++': 'M. tuberculosis',
    '+++-+': 'M. africanum (lineage 5)',
    '++--+': 'M. africanum (lineage 6), M. microti, M. pinnipedii or M. mungi',
    '++---': 'M. orygis or M. caprae',
    '+----': 'M. bovis',
    '-----': 'M. bovis BCG',
    '-+--+': 'Dassie bacillus',
}


def rd_profile(check):
    """e.g. "+----" for RD1 present, RD4, RD7, RD9 and RD12 deleted; "?" for a region not measured."""
    return ''.join({PRESENT: '+', DELETED: '-'}.get(check.state(r), '?') for r in REGIONS)


def call_species(check, lineages=(), spacers=True):
    """
    :param lineages: lineages called by the SNP barcode
    :param spacers: whether any of the 43 standard spacers was found (M. canettii usually has none)
    """
    main = {lineage.split('.')[0] for lineage in lineages}
    if any(check.state(r) == PARTIAL for r in REGIONS):
        return 'MTBC (mixed or unclear RD profile)'
    profile = rd_profile(check)
    rd1, rd4, rd7, rd9, rd12 = profile
    if rd7 == '+' and (rd4 == '-' or rd12 == '-'):  # RD4 or RD12 lost independently of the M. bovis lineage
        return 'M. canettii'
    species = RD_PROFILES.get(profile)
    if species is None:
        return 'MTBC (unusual RD profile: {})'.format(', '.join(
            '{}{}'.format(r, s) for r, s in zip(REGIONS, profile, strict=True)))
    if profile == '+++++':
        # Lineage 4 is the only lineage defined by the H37Rv allele, which some M. canettii strains carry: only a
        # sublineage of lineage 4, or lineages 1, 2, 3 or 7, are specific to M. tuberculosis
        specific = [lin for lin in lineages if lin[0] in '1237' or lin.startswith('4.')]
        if specific:
            return 'M. tuberculosis'
        if not spacers:
            return 'M. canettii'
        if main == {'4'}:
            return 'M. tuberculosis'
        if not main:  # Outside lineages 1 to 7 and the animal lineages
            return 'MTBC, RD9 intact, no lineage (e.g. M. canettii)'
        return 'MTBC (RD9 intact, unexpected for lineage {})'.format('/'.join(sorted(main)))
    if profile == '++--+':
        if '6' in main:
            return 'M. africanum (lineage 6)'
        if 'BOV_AFRI' in main:  # The clade of the animal lineages and lineage 6, without the lineage 6 SNP
            return 'M. microti, M. pinnipedii or M. mungi'
    return species


def consistency_warnings(check, lineages):
    """Contradictions between the RD profile and the lineage SNPs."""
    warnings = []
    main = {lineage.split('.')[0] for lineage in lineages}
    lineage_text = '/'.join(sorted(main))
    if check.species == 'M. canettii' and main:
        warnings.append('lineage SNPs ({}) found in a sample identified as M. canettii: the SNP barcode is not '
                        'designed for M. canettii'.format(lineage_text))
        return warnings
    profile = rd_profile(check)
    rd4, rd7, rd9 = check.state('RD4'), check.state('RD7'), check.state('RD9')
    expected_rd9 = DELETED if main & {'5', '6', 'BOV', 'BOV_AFRI'} else PRESENT if main else ''
    if expected_rd9 and rd9 in (PRESENT, DELETED) and rd9 != expected_rd9:
        warnings.append('RD9 is {} but the lineage SNPs indicate lineage {}'.format(rd9, lineage_text))
    if 'BOV' in main and rd4 == PRESENT and profile != '++---':  # M. caprae and M. orygis carry the BOV SNP
        warnings.append('the lineage SNPs indicate M. bovis but RD4 is present')
    if rd4 == DELETED and main - {'BOV', 'BOV_AFRI'}:
        warnings.append('RD4 is deleted (M. bovis) but the lineage SNPs indicate lineage {}'.format(lineage_text))
    expected_rd7 = PRESENT if '5' in main else DELETED if main & {'6', 'BOV', 'BOV_AFRI'} else ''
    if expected_rd7 and rd7 in (PRESENT, DELETED) and rd7 != expected_rd7:
        warnings.append('RD7 is {} but the lineage SNPs indicate lineage {}'.format(rd7, lineage_text))
    return warnings
