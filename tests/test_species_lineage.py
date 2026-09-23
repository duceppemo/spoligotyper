import json

import pytest

from spoligotyper import lineage, species
from spoligotyper.pipeline import Result, write_multiqc
from spoligotyper.spoligotype import closest, describe_closest, load_database

from .conftest import SB0140


def marker_counts(control=30, rd9=30, rd4=30, rd1=30):
    counts = {'MTBC_{:02d}'.format(i): control for i in range(1, 41)}
    for region, value in (('RD9', rd9), ('RD4', rd4), ('RD1', rd1)):
        counts.update({'{}_{:02d}'.format(region, i): value for i in range(1, 9)})
    return counts


@pytest.mark.parametrize('rd9, rd4, rd1, lineages, expected', [
    (30, 30, 30, ['4', '4.9'], 'M. tuberculosis'),
    (30, 30, 30, [], 'MTBC, RD9 intact, no lineage (e.g. M. canettii)'),
    (0, 0, 30, ['BOV', 'BOV_AFRI'], 'M. bovis'),
    (0, 0, 0, ['BOV', 'BOV_AFRI'], 'M. bovis BCG'),
    (0, 30, 30, ['6', 'BOV_AFRI'], 'M. africanum'),
    (0, 30, 30, ['BOV_AFRI'], 'Animal-adapted MTBC, not M. bovis (e.g. M. caprae, M. pinnipedii)'),
    (0, 30, 0, ['BOV_AFRI'], 'Animal-adapted MTBC, not M. bovis (RD1 deleted, e.g. M. microti)'),
    (9, 30, 30, ['4'], 'MTBC (mixed or unclear RD profile)'),
])
def test_species(rd9, rd4, rd1, lineages, expected):
    check = species.check_species(marker_counts(30, rd9, rd4, rd1), 'fastq', lineages=lineages)
    assert check.mtbc and check.species == expected


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


def test_consistency_warnings():
    check = species.check_species(marker_counts(30, 30, 30, 30), 'fastq', lineages=['BOV', 'BOV_AFRI'])
    warnings = species.consistency_warnings(check, ['BOV', 'BOV_AFRI'])
    assert any('RD9 is present' in w for w in warnings) and any('RD4 is present' in w for w in warnings)
    check = species.check_species(marker_counts(30, 0, 0, 30), 'fastq', lineages=['4'])
    assert len(species.consistency_warnings(check, ['4'])) == 2


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
    (('4', '4.9', 'BOV', 'BOV_AFRI'), 'BOV', 'M. bovis'),
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


def test_on_one_path():
    assert lineage.on_one_path(['4', '4.3', '4.3.4'])
    assert lineage.on_one_path(['BOV', 'BOV_AFRI']) and lineage.on_one_path(['6', 'BOV_AFRI'])
    assert not lineage.on_one_path(['4.1', '4.3']) and not lineage.on_one_path(['2', '4'])
    assert not lineage.on_one_path(['BOV', '6'])


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
    write_multiqc([Result('S1', spoligotype='Spoligo not found', octal='000000000003771')], tmp_path / 'x_mqc.json')
    content = json.loads((tmp_path / 'x_mqc.json').read_text())
    assert content['id'] == 'spoligotyper' and content['plot_type'] == 'table'
    assert content['data'] == {'S1': {'Spoligotype': 'Spoligo not found', 'Octal': '000000000003771',
                                      'Species': '-', 'Lineage': '-', 'Status': 'ok'}}
