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
from functools import cache

from .spoligotype import data_file

MARKERS_FASTA = data_file('markers.fasta')
CONTROL = 'MTBC'
REGIONS = ('RD1', 'RD4', 'RD7', 'RD9', 'RD12')
PRESENT, DELETED, PARTIAL, REDUCED = 'present', 'deleted', 'partial', 'reduced'
SEGMENT_FOUND = 0.05  # A segment is found when its depth is at least this fraction of the control depth. Deleted
# segments have no read at all, while GC-rich segments (ESX region of RD1) can drop to 10% of the control depth in
# real Illumina reads.
MISSING_TOLERANCE = 0.1  # In reads, a present region may miss this fraction of its segments (low depth), at least 1;
# a deleted region may have as many segments found (stray reads)
REDUCED_RATIO = 0.5  # Depth of a present region below this fraction of the control depth: mixed sample?
# The regions of difference on the H37Rv genome (NC_000962.3), as found by scripts/make_reference_data.py: the
# H37Rv segments missing from M. bovis AF2122/97 (RD4, RD7, RD9, RD12) and from BCG Pasteur (RD1).
REGION_EXTENTS = {'RD1': (4350251, 4359740), 'RD4': (1696001, 1708740), 'RD7': (2208001, 2220740),
                  'RD9': (2330051, 2332140), 'RD12': (3485101, 3487540)}
# M. microti lost part of RD1 with its own deletion, RD1mic (Brodin et al. 2002), which ends in Rv3876: in three
# M. microti genomes, the RD1 segments up to H37Rv 4,354,450 are missing and those from 4,354,851 are present.
RD1MIC = (4348827, 4354800)
MIN_CONTROL_READS = 3  # Median reads per control chunk to call the species from reads
MIN_CONTROL_FRACTION = 0.5  # Fraction of the control chunks found
CHUNK = 100


@dataclass
class RegionCall:
    state: str  # PRESENT, DELETED, PARTIAL (some segments missing) or REDUCED (all found, at low depth)
    found: int  # Segments found
    total: int
    ratio: float  # Median depth of the segments found, relative to the control depth (0 if none)
    missing: list = field(default_factory=list)  # H37Rv ranges of the missing segments: [(start, end)]

    @property
    def fraction(self):
        return self.found / self.total if self.total else 0.0

    @property
    def sign(self):
        """+ or - in the RD profile: a partially deleted region counts as present if most of its segments are found."""
        if self.state == PARTIAL:
            return '+' if self.fraction >= 0.5 else '-'
        return '-' if self.state == DELETED else '+'

    def describe(self):
        """e.g. "partially deleted: 11 of 20 segments found, H37Rv 4,350,651-4,354,950 missing"."""
        if self.state == PARTIAL:
            return 'partially deleted: {} of {} segments found, H37Rv {} missing'.format(
                self.found, self.total, ', '.join('{:,}-{:,}'.format(a, b) for a, b in self.missing))
        if self.state == REDUCED:
            return 'present at reduced depth ({:.2f} of the control): mixed sample?'.format(self.ratio)
        return self.state


@dataclass
class SpeciesCheck:
    mtbc: bool = False  # Enough MTBC DNA to call the regions
    control_depth: float = 0.0  # Median reads (or contigs) per control chunk
    control_found: float = 0.0  # Fraction of the control chunks found
    regions: dict = field(default_factory=dict)  # {region: RegionCall}
    mtbc_fraction: float | None = None  # Estimated fraction of the reads from MTBC, reads only
    species: str = ''

    def state(self, region):
        return self.regions[region].state if region in self.regions else ''

    def summary(self):
        """e.g. "RD9 deleted, RD4 deleted, RD1 present"."""
        return ', '.join('{} {}'.format(r, self.state(r)) for r in REGIONS if r in self.regions)


def chunk_counts(counts, prefix):
    return [c for name, c in counts.items() if name.startswith(prefix + '_')]


@cache
def segments(path=MARKERS_FASTA):
    """{marker name: (H37Rv start, end)}, from the descriptions of markers.fasta ("RD1_01 H37Rv:4350651-4350750")."""
    result = {}
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                name, location = line[1:].split()[:2]
                start, end = location.split(':')[1].split('-')
                result[name] = (int(start), int(end))
    return result


def missing_stretches(names, found):
    """
    H37Rv ranges of the missing segments: consecutive missing segments (no segment found between them) form one
    stretch, from the start of the first to the end of the last.
    """
    stretches, current = [], None
    for name in sorted(names, key=lambda n: segments()[n]):
        start, end = segments()[name]
        if name in found:
            current = None
        elif current is None:
            current = [start, end]
            stretches.append(current)
        else:
            current[1] = end
    return [tuple(s) for s in stretches]


def call_region(counts, region, control_depth, file_type='fastq'):
    """
    Presence of a region from its segments, found when their depth is at least 5% of the control depth. In an
    assembly, a missing segment is a real absence: the region is present only if all its segments are found.
    """
    names = sorted(name for name in segments() if name.startswith(region + '_'))
    if not names:
        return None
    threshold = max(1.0, SEGMENT_FOUND * control_depth)
    found = [name for name in names if counts.get(name, 0) >= threshold]
    tolerance = max(1, round(MISSING_TOLERANCE * len(names)))
    ratio = statistics.median(counts[name] for name in found) / control_depth if found and control_depth else 0.0
    missing = missing_stretches(names, set(found))
    if len(found) <= tolerance:
        state = DELETED
    elif len(found) < len(names) - (0 if file_type == 'fasta' else tolerance):
        state = PARTIAL
    else:  # Contig counts in an assembly say nothing about a mixed sample
        state = REDUCED if ratio < REDUCED_RATIO and file_type != 'fasta' else PRESENT
    return RegionCall(state, len(found), len(names), ratio, missing if state == PARTIAL else [])


def rd1mic(check):
    """
    True when RD1 lacks its segments in the RD1mic deletion of M. microti (Rv3871 to Rv3876) and has the others. As
    for the regions, one error in ten (at least one) is tolerated on each side: in real reads, a segment inside RD1mic
    can get a few stray reads, and a GC-rich segment outside it can get none.
    """
    call = check.regions.get('RD1')
    if call is None or call.state != PARTIAL:
        return False
    names = [n for n in segments() if n.startswith('RD1_')]
    inside = {n for n in names if segments()[n][0] >= RD1MIC[0] and segments()[n][1] <= RD1MIC[1]}
    missing = {n for n in names if any(start <= segments()[n][0] and segments()[n][1] <= end
                                       for start, end in call.missing)}
    outside = set(names) - inside
    return (bool(inside) and len(inside - missing) <= max(1, round(MISSING_TOLERANCE * len(inside)))
            and len(missing & outside) <= max(1, round(MISSING_TOLERANCE * len(outside))))


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
        call = call_region(counts, region, check.control_depth, file_type)
        if call is not None:
            check.regions[region] = call
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
    """
    e.g. "+----" for RD1 present, RD4, RD7, RD9 and RD12 deleted; "?" for a region not measured. A partially deleted
    region is "+" when at least half of its segments are found (e.g. RD1 of M. microti), "-" otherwise.
    """
    return ''.join(check.regions[r].sign if r in check.regions else '?' for r in REGIONS)


def call_species(check, lineages=(), spacers=True):
    """
    :param lineages: lineages called by the SNP barcode
    :param spacers: whether any of the 43 standard spacers was found (M. canettii usually has none)
    """
    main = {lineage.split('.')[0] for lineage in lineages}
    if any(check.state(r) == REDUCED for r in REGIONS):
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
        if rd1mic(check):
            return 'M. microti'
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
