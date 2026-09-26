import hashlib
import json

import pytest

from spoligotyper import sitdb
from spoligotyper.spoligotype import NOT_FOUND, closest

from .conftest import H37RV

H37RV_NO = H37RV.replace('1', 'n').replace('0', 'o')
SOURCE_CSV = '\n'.join([
    'StrainID;Spoligo Binary;Spoligo Octal;Lineage (SITVIT2);SIT;Country Distribution (SITVIT2)',
    "A1;{};'777777477760771;T-H37Rv;451;[US=1]".format(H37RV_NO),
    "A2;{};'000000000003771;Beijing;1;[CN=9]".format('o' * 34 + 'n' * 9),
    "A3;{};'777777777777771;Unknown;Orphan;[FR=1]".format('n' * 43),
    "A4;{}?;'777777777777770;Unknown;Orphan;[IT=1]".format('n' * 42),  # Unknown spacer: skipped
]) + '\n'


@pytest.fixture
def source(tmp_path, monkeypatch):
    path = tmp_path / 'Spoligo_list.csv'
    path.write_text(SOURCE_CSV)
    monkeypatch.setitem(sitdb.SOURCE, 'sha256', hashlib.sha256(SOURCE_CSV.encode()).hexdigest())
    return path


def test_convert_and_lookup(source, tmp_path):
    database = tmp_path / 'sit.tsv'
    assert sitdb.convert(source, database) == 3
    db = sitdb.load(database)
    assert db.lookup(H37RV) == ('SIT451', 'T-H37Rv')
    assert db.lookup('1' * 43) == ('Orphan', 'Unknown')
    assert db.lookup('0' * 43) == (NOT_FOUND, '')
    assert db.sits == {H37RV: 'SIT451', '0' * 34 + '1' * 9: 'SIT1'}
    assert closest(H37RV[:-1] + '0', db.sits) == [('SIT451', [43])]
    assert db.source.startswith('SIT database for spoligotyper') and len(db.sha256) == 64


def test_convert_rejects_wrong_octal(source, tmp_path):
    source.write_text(SOURCE_CSV.replace("'777777477760771", "'777777477760770"))
    with pytest.raises(sitdb.SitError, match='octal code does not match'):
        sitdb.convert(source, tmp_path / 'sit.tsv')


def test_download_falls_back_and_checks_checksum(source, tmp_path):
    wrong = tmp_path / 'wrong.csv'
    wrong.write_text('not the list\n')
    urls = ['file:///nonexistent/Spoligo_list.csv', wrong.as_uri(), source.as_uri()]
    assert sitdb.download(tmp_path / 'downloaded.csv', urls) == source.as_uri()
    with pytest.raises(sitdb.SitError, match='checksum'):
        sitdb.download(tmp_path / 'downloaded.csv', urls[:2])


def test_download_command(source, tmp_path, monkeypatch):
    monkeypatch.setattr(sitdb, 'URLS', [source.as_uri()])
    sitdb.main([])  # Default folder: $SPOLIGOTYPER_DATA (set by conftest)
    assert sitdb.default_database() is not None
    info = json.loads((sitdb.default_folder() / 'sit_database.json').read_text())
    assert info['patterns'] == 3 and info['origin'] == source.as_uri() and info['license'] == 'GPL-3.0'
    sitdb.main(['-o', str(tmp_path / 'offline'), '--source', str(source)])  # Computer without internet access
    assert sitdb.load(tmp_path / 'offline' / sitdb.DB_NAME).lookup(H37RV)[0] == 'SIT451'


def test_download_command_bad_source(source, tmp_path):
    source.write_text('tampered\n')
    with pytest.raises(SystemExit) as e:
        sitdb.main(['-o', str(tmp_path / 'db'), '--source', str(source)])
    assert e.value.code == 1


def test_no_default_database():
    assert sitdb.default_database() is None  # conftest points SPOLIGOTYPER_DATA to an empty folder


def test_load_rejects_other_files(tmp_path):
    other = tmp_path / 'other.tsv'
    other.write_text('a\tb\n')
    with pytest.raises(sitdb.SpoligoError):
        sitdb.load(other)
