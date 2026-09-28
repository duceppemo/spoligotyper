import json

import pytest

from spoligotyper import lineage, species
from spoligotyper.pipeline import Result, write_multiqc
from spoligotyper.spoligotype import closest, describe_closest, load_database

from .conftest import SB0140


def marker_counts(control=30, rd9=30, rd4=30, rd1=30, rd7=30, rd12=30):
    """Counts for every segment of markers.fasta: the same count for all the segments of a region."""
    values = {'MTBC': control, 'RD9': rd9, 'RD4': rd4, 'RD1': rd1, 'RD7': rd7, 'RD12': rd12}
    return {name: values[name.split('_')[0]] for name in species.segments()}


def profile_counts(profile):
    """Counts for an RD profile written as in the RD PCR table: RD1, RD4, RD7, RD9, RD12, "+" present, "-" deleted."""
    values = [30 if c == '+' else 0 for c in profile]
    return marker_counts(30, rd1=values[0], rd4=values[1], rd7=values[2], rd9=values[3], rd12=values[4])


@pytest.mark.parametrize('profile, lineages, spacers, expected', [
    ('+++++', ['4', '4.9'], True, 'M. tuberculosis'),
    ('+++++', ['2', '2.2', '2.2.1'], False, 'M. tuberculosis'),  # DR locus deleted: still M. tuberculosis
    ('+++++', [], True, 'MTBC, RD9 intact, no lineage (e.g. M. canettii)'),
    ('+++++', [], False, 'M. canettii'),  # No standard spacer
    ('+++++', ['4'], False, 'M. canettii'),  # Some M. canettii carry the H37Rv (lineage 4) allele
    ('+++++', ['4'], True, 'M. tuberculosis'),
    ('+++-+', ['5'], True, 'M. africanum (lineage 5)'),
    ('++--+', ['6', 'BOV_AFRI'], True, 'M. africanum (lineage 6)'),
    ('++--+', ['BOV_AFRI'], True, 'M. microti, M. pinnipedii or M. mungi'),
    ('++--+', [], True, 'M. africanum (lineage 6), M. microti, M. pinnipedii or M. mungi'),
    ('++++-', [], False, 'M. canettii'),  # RD12 deleted, RD7 present
    ('+-+++', [], False, 'M. canettii'),  # RD4 deleted, RD7 present
    ('+-+-+', [], False, 'M. canettii'),
    ('++---', ['BOV', 'BOV_AFRI'], True, 'M. orygis or M. caprae'),
    ('+----', ['BOV', 'BOV_AFRI'], True, 'M. bovis'),
    ('-----', ['BOV', 'BOV_AFRI'], True, 'M. bovis BCG'),
    ('-+--+', [], True, 'Dassie bacillus'),
    ('+-++-', [], True, 'M. canettii'),
    ('-++++', [], True, 'MTBC (unusual RD profile: RD1-, RD4+, RD7+, RD9+, RD12+)'),
])
def test_species(profile, lineages, spacers, expected):
    """The RD PCR table (RD1, RD4, RD7, RD9, RD12), refined with the lineage SNPs and the spacers."""
    check = species.check_species(profile_counts(profile), 'fastq')
    assert check.mtbc and species.rd_profile(check) == profile
    species.name_species(check, lineages, spacers=spacers)
    assert check.species == expected


def test_species_reduced_depth():
    """All the RD9 segments found, but at 30% of the control depth: a mix of strains with and without RD9."""
    check = species.check_species(marker_counts(30, rd9=9), 'fastq')
    assert check.regions['RD9'].state == species.REDUCED and check.regions['RD9'].sign == '+'
    assert check.species == 'MTBC (mixed or unclear RD profile)'


def test_short_region_tolerance():
    """RD9 has 8 segments: one missing at low depth is still present, one stray segment is still deleted."""
    rd9 = sorted(n for n in species.segments() if n.startswith('RD9_'))
    counts = marker_counts(30)
    counts[rd9[0]] = 0
    assert species.check_species(counts, 'fastq').regions['RD9'].state == species.PRESENT
    assert species.check_species(counts, 'fasta').regions['RD9'].state == species.PARTIAL
    counts = marker_counts(30, rd9=0)
    counts[rd9[0]] = 2
    assert species.check_species(counts, 'fastq').regions['RD9'].state == species.DELETED


def test_no_reduced_depth_in_assemblies():
    """Contig counts do not measure depth: a region found once with controls found 3 times is present."""
    check = species.check_species(marker_counts(3, rd9=1), 'fasta')
    assert check.regions['RD9'].state == species.PRESENT
    assert species.call_region(marker_counts(0), 'RD9', 0).ratio == 0.0


def rd1_segments():
    return sorted((n for n in species.segments() if n.startswith('RD1_')), key=lambda n: species.segments()[n])


def test_partial_deletion_and_rd1mic():
    """M. microti: the RD1 segments inside RD1mic (Rv3871 to part of Rv3876) are missing, the others present."""
    counts = marker_counts(30, rd7=0, rd9=0)
    inside = [n for n in rd1_segments() if species.segments()[n][1] <= species.RD1MIC[1]]
    counts.update({n: 0 for n in inside})
    check = species.check_species(counts, 'fastq')
    rd1 = check.regions['RD1']
    assert rd1.state == species.PARTIAL and rd1.found == 20 - len(inside) and rd1.sign == '+'
    assert rd1.missing == [(species.segments()[inside[0]][0], species.segments()[inside[-1]][1])]
    assert species.rd1mic(check) and species.rd_profile(check) == '++--+'
    species.name_species(check, ['BOV_AFRI'])
    assert check.species == 'M. microti'
    assert 'H37Rv 4,350,651-' in rd1.describe()


def test_rd1mic_in_real_reads():
    """M. microti reads (ERR027297): stray reads on one RD1mic segment, a GC-rich segment outside it with none."""
    counts = marker_counts(150, rd1=150, rd4=150, rd12=150, rd7=0, rd9=0)
    names = rd1_segments()
    inside = [n for n in names if species.segments()[n][1] <= species.RD1MIC[1]]
    counts.update({n: 0 for n in inside})
    counts[inside[4]] = 8  # Above 5% of the control depth: found
    counts[names[17]] = 0  # Outside RD1mic
    check = species.check_species(counts, 'fastq')
    assert species.rd1mic(check) and check.species == 'M. microti'
    counts[inside[5]] = 8  # Two segments of RD1mic found: not RD1mic
    assert not species.rd1mic(species.check_species(counts, 'fastq'))


def test_partial_deletion_other():
    """Three RD1 segments missing, not RD1mic (M. mungi genome): the species stays the group of three."""
    counts = marker_counts(30, rd7=0, rd9=0)
    names = rd1_segments()
    counts.update({n: 0 for n in names[3:6]})
    check = species.check_species(counts, 'fastq')
    assert check.regions['RD1'].state == species.PARTIAL and not species.rd1mic(check)
    assert check.regions['RD1'].missing == [(species.segments()[names[3]][0], species.segments()[names[5]][1])]
    species.name_species(check, ['BOV_AFRI'])
    assert check.species == 'M. microti, M. pinnipedii or M. mungi'
    # Most of RD4 missing (M. canettii ET1291): "-" in the RD profile
    counts = marker_counts(30)
    counts.update({n: 0 for n in sorted(n for n in species.segments() if n.startswith('RD4_'))[:11]})
    check = species.check_species(counts, 'fastq')
    assert check.regions['RD4'].state == species.PARTIAL and check.regions['RD4'].sign == '-'


def test_gc_rich_segments_at_low_depth():
    """Real Illumina reads (AF2122/97, ERR1744454): GC-rich RD1 segments at 10% of the control depth are found."""
    counts = marker_counts(80, rd1=80, rd4=0, rd7=0, rd9=0, rd12=0)
    counts.update({n: 8 for n in rd1_segments()[4:8]})
    check = species.check_species(counts, 'fastq')
    assert check.regions['RD1'].state == species.PRESENT
    counts.update({n: 3 for n in rd1_segments()[4:8]})  # Below 5% of the control depth: missing
    assert species.check_species(counts, 'fastq').regions['RD1'].state == species.PARTIAL


def test_assembly_needs_all_segments():
    """In an assembly, one missing segment makes the region partial; reads tolerate 10% (low depth)."""
    counts = marker_counts(1)
    counts[rd1_segments()[-1]] = 0
    assert species.check_species(counts, 'fasta').regions['RD1'].state == species.PARTIAL
    counts = marker_counts(30)
    counts[rd1_segments()[-1]] = 0
    assert species.check_species(counts, 'fastq').regions['RD1'].state == species.PRESENT


def test_region_extents():
    """Every segment lies inside its region of difference."""
    for name, (start, end) in species.segments().items():
        region = name.split('_')[0]
        if region != 'MTBC':
            assert species.REGION_EXTENTS[region][0] <= start and end <= species.REGION_EXTENTS[region][1], name


def test_species_not_mtbc():
    check = species.check_species(marker_counts(control=1), 'fastq')
    assert not check.mtbc and check.species == 'MTBC not detected' and check.regions == {}
    assert species.check_species(marker_counts(control=1), 'fasta').mtbc  # One contig is enough for an assembly
    assert species.check_species({}, 'fastq').species == ''


def test_mtbc_fraction():
    # 150 bp single-end reads at 30x: about 30 * 201 / 150 = 40 reads per 100 bp control chunk
    assert species.check_species(marker_counts(40), 'fastq', depth=30, read_length=150).mtbc_fraction == \
        pytest.approx(0.995, abs=0.01)
    assert species.check_species(marker_counts(20), 'fastq', depth=30, read_length=150).mtbc_fraction == \
        pytest.approx(0.5, abs=0.01)
    # Paired-end: both reads of a pair are counted
    assert species.check_species(marker_counts(80), 'fastq', depth=30, read_length=150, paired=True) \
        .mtbc_fraction == pytest.approx(0.995, abs=0.01)
    assert species.check_species(marker_counts(), 'fasta').mtbc_fraction is None


def test_expected_control_reads():
    assert species.expected_control_reads(30, 150, paired=False) == pytest.approx(40.2)
    assert species.expected_control_reads(30, 150, paired=True) == pytest.approx(80.4)


def test_consistency_warnings():
    check = species.check_species(profile_counts('+++++'), 'fastq', lineages=['BOV', 'BOV_AFRI'])
    warnings = species.consistency_warnings(check, ['BOV', 'BOV_AFRI'])
    assert any('RD9 is present' in w for w in warnings) and any('RD4 is present' in w for w in warnings)
    assert any('RD7 is present' in w for w in warnings)
    check = species.check_species(profile_counts('+----'), 'fastq', lineages=['4'])
    assert len(species.consistency_warnings(check, ['4'])) == 2
    assert species.consistency_warnings(check, ['BOV_AFRI']) == []  # Consistent with M. bovis
    assert species.consistency_warnings(check, ['BOV', 'BOV_AFRI']) == []
    check = species.check_species(profile_counts('++---'), 'fastq', lineages=['BOV', 'BOV_AFRI'])
    assert species.consistency_warnings(check, ['BOV', 'BOV_AFRI']) == []  # M. caprae, M. orygis carry the BOV SNP
    check = species.check_species(profile_counts('+++-+'), 'fastq', lineages=['6'])
    assert any('RD7 is present' in w for w in species.consistency_warnings(check, ['6']))
    check = species.check_species(profile_counts('+++++'), 'fastq')
    species.name_species(check, ['4'], spacers=False)
    assert any('not designed for M. canettii' in w for w in species.consistency_warnings(check, ['4']))


def snp_counts(alt_lineages=(), reads=20, mixed=None):
    """Counts as Seal reports them. mixed: {lineage: fraction of reads with the alternative allele}."""
    counts = {}
    for row in lineage.read_barcode():
        key = '{}|{}|'.format(row['lineage'], row['position'])
        fraction = (mixed or {}).get(row['lineage'], 1.0 if row['lineage'] in alt_lineages else 0.0)
        counts[key + 'alt'] = round(reads * fraction)
        counts[key + 'ref'] = reads - counts[key + 'alt']
    return counts


@pytest.mark.parametrize('alt, expected, name', [
    ((), '4.9', 'Euro-American (H37Rv-like)'),  # H37Rv: the reference alleles
    (('4', '4.9', '2', '2.2', '2.2.1'), '2.2.1', 'East-Asian'),
    (('4', '4.9', '1', '1.1', '1.1.2'), '1.1.2', 'Indo-Oceanic'),
    (('4', '4.9', 'BOV', 'BOV_AFRI'), 'BOV', 'M. bovis, M. caprae, M. orygis'),
    (('4', '4.9', '6', 'BOV_AFRI'), '6', 'West-Africa 2'),
    (('4.9', '4.3', '4.3.4', '4.3.4.2'), '4.3.4.2', 'Euro-American (LAM)'),
])
def test_lineage(alt, expected, name):
    call = lineage.call_lineage(snp_counts(alt), 'fastq')
    assert (call.lineage, call.name, call.conflict, call.mixed) == (expected, name, False, [])


def test_lineage_thresholds():
    assert lineage.call_lineage(snp_counts(reads=2), 'fastq').lineage == ''  # Too few reads
    assert lineage.call_lineage(snp_counts(reads=1), 'fasta').lineage == '4.9'  # One contig is enough
    assert lineage.call_lineage({}, 'fastq').lineage == ''


def test_lineage_mixed_and_conflict():
    call = lineage.call_lineage(snp_counts(mixed={'4': 0.3, '4.9': 0.3, 'BOV': 0.3, 'BOV_AFRI': 0.3}), 'fastq')
    assert call.lineage.startswith('mixed: ') and len(call.mixed) == 4
    call = lineage.call_lineage(snp_counts(('2', '2.2')), 'fastq')  # 2 and the H37Rv-like 4.9: impossible
    assert call.conflict and call.lineage == 'mixed: 2, 2.2, 4, 4.9'
    call = lineage.call_lineage(snp_counts(('BOV_AFRI',)), 'fastq')  # H37Rv-like with the BOV_AFRI allele
    assert call.conflict and call.lineage == 'mixed: 4, 4.9, BOV_AFRI'
    call = lineage.call_lineage(snp_counts(('4', '4.9', 'BOV_AFRI')), 'fastq')  # The BOV SNP not covered
    assert not call.conflict and call.lineage == 'BOV_AFRI'


def test_on_one_path():
    assert lineage.on_one_path(['4', '4.3', '4.3.4'])
    assert lineage.on_one_path(['BOV', 'BOV_AFRI']) and lineage.on_one_path(['6', 'BOV_AFRI'])
    assert not lineage.on_one_path(['4.1', '4.3']) and not lineage.on_one_path(['2', '4'])
    assert not lineage.on_one_path(['BOV', '6'])
    assert lineage.on_one_path(['BOV_AFRI'])
    for other in (['4', '4.9'], ['5'], ['1'], ['7'], ['2', '2.2']):
        assert not lineage.on_one_path(other + ['BOV_AFRI']), other


def test_closest():
    db = load_database()
    assert closest(SB0140, db) == []
    one_off = SB0140[:5] + '1' + SB0140[6:20] + '0' + SB0140[21:]  # Two changes from SB0140: not in the database
    matches = closest(one_off, db)
    assert matches and all(len(diff) == 1 for _, diff in matches)
    assert describe_closest([('SB0001', [7]), ('SB0002', [7, 13])]) == \
        'SB0001 (spacer 7 differs); SB0002 (spacers 7, 13 differ)'
    assert closest('0' * 20 + '1' * 23, db, max_distance=0) == []


def test_multiqc(tmp_path):
    write_multiqc([Result('S1', sb='Not in Mbovis.org', octal='000000000003771')], tmp_path / 'x_mqc.json')
    content = json.loads((tmp_path / 'x_mqc.json').read_text())
    assert content['id'] == 'spoligotyper' and content['plot_type'] == 'table'
    assert content['data'] == {'S1': {'SB': 'Not in Mbovis.org', 'SIT': '-', 'Octal': '000000000003771',
                                      'Species': '-', 'Lineage': '-', 'Status': 'ok'}}
