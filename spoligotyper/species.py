"""
Species check from regions of difference (RD), and amount of MTBC DNA in the sample.

markers.fasta holds 100 bp chunks of the H37Rv genome (see scripts/make_reference_data.py):
  MTBC  control chunks, found in all MTBC genomes and in no non-tuberculous mycobacteria
  RD9   deleted in M. africanum and the animal-adapted lineages, including M. bovis
  RD4   deleted in M. bovis and BCG only
  RD1   deleted in BCG (and, with a different deletion, in M. microti)

Each region's read depth is compared with the control depth: about the same when the region is present, 0 when it
is deleted, and in between for a mix of strains with and without it.
"""

import statistics
from dataclasses import dataclass, field

from .spoligotype import data_file

MARKERS_FASTA = data_file('markers.fasta')
CONTROL = 'MTBC'
REGIONS = ('RD9', 'RD4', 'RD1')
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


def name_species(check, lineages=(), mixed=False):
    """Name the species once the lineage is known (check_species runs before the lineage call)."""
    if check.mtbc:
        check.species = 'MTBC, mixed sample?' if mixed else call_species(check, lineages)


def call_species(check, lineages=()):
    main = {lineage.split('.')[0] for lineage in lineages}
    rd9, rd4, rd1 = (check.state(r) for r in REGIONS)
    if PARTIAL in (rd9, rd4, rd1):
        return 'MTBC (mixed or unclear RD profile)'
    if rd9 == PRESENT:
        if main & {'1', '2', '3', '4', '7'}:
            return 'M. tuberculosis'
        if not main:  # Outside lineages 1 to 7 and the animal lineages
            return 'MTBC, RD9 intact, no lineage (e.g. M. canettii)'
        return 'MTBC (RD9 intact, unexpected for lineage {})'.format('/'.join(sorted(main)))
    if rd4 == DELETED:
        return 'M. bovis BCG' if rd1 == DELETED else 'M. bovis'
    if main & {'5', '6'}:
        return 'M. africanum'
    if rd1 == DELETED:
        return 'Animal-adapted MTBC, not M. bovis (RD1 deleted, e.g. M. microti)'
    if 'BOV_AFRI' in main:  # The clade of M. bovis and lineage 6, without the SNPs of either
        return 'Animal-adapted MTBC, not M. bovis (e.g. M. caprae, M. pinnipedii)'
    return 'M. africanum or animal-adapted MTBC, not M. bovis (RD9 deleted, RD4 present)'


def consistency_warnings(check, lineages):
    """Contradictions between the RD profile and the lineage SNPs."""
    warnings = []
    main = {lineage.split('.')[0] for lineage in lineages}
    rd9, rd4 = check.state('RD9'), check.state('RD4')
    expected_rd9 = DELETED if main & {'5', '6', 'BOV', 'BOV_AFRI'} else PRESENT if main else ''
    if expected_rd9 and rd9 in (PRESENT, DELETED) and rd9 != expected_rd9:
        warnings.append('RD9 is {} but the lineage SNPs indicate lineage {}'.format(rd9, '/'.join(sorted(main))))
    if 'BOV' in main and rd4 == PRESENT:
        warnings.append('the lineage SNPs indicate M. bovis but RD4 is present')
    if rd4 == DELETED and main - {'BOV', 'BOV_AFRI'}:
        warnings.append('RD4 is deleted (M. bovis) but the lineage SNPs indicate lineage {}'.format(
            '/'.join(sorted(main))))
    return warnings
