import pytest

from spoligotyper.spoligotype import (
    NOT_FOUND,
    SPOLIGOTYPE_DB,
    SpoligoError,
    binary_to_hex,
    binary_to_octal,
    load_database,
    lookup,
    read_spacer_names,
    to_binary,
)

from .conftest import H37RV, SB0140


def test_spacer_names():
    names = read_spacer_names()
    assert names == ['spacer{:02d}'.format(i) for i in range(1, 44)]


def test_octal():
    assert binary_to_octal(SB0140) == '664073777777600'
    assert binary_to_octal(H37RV) == '777777477760771'
    assert binary_to_octal('0' * 43) == '0' * 15
    assert binary_to_octal('1' * 43) == '7' * 14 + '1'


def test_hex():
    assert binary_to_hex(SB0140) == '6D-03-5F-7F-FF-60'
    assert binary_to_hex(H37RV) == '7F-7F-7C-7F-F0-7F'
    assert binary_to_hex('0' * 43) == '00-00-00-00-00-00'
    assert binary_to_hex('1' * 43) == '7F-7F-7F-7F-FF-7F'
    assert binary_to_hex('0000101' + '0' * 36) == '05-00-00-00-00-00'


@pytest.mark.parametrize('binary', ['1' * 42, '1' * 44, '2' + '1' * 42, ''])
def test_invalid_binary(binary):
    with pytest.raises(SpoligoError):
        binary_to_octal(binary)
    with pytest.raises(SpoligoError):
        binary_to_hex(binary)


def test_to_binary():
    names = read_spacer_names()
    counts = {name: 5 for name in names[:3]}
    counts[names[3]] = 4
    assert to_binary(counts, names, 5) == '111' + '0' * 40
    assert to_binary(counts, names, 4) == '1111' + '0' * 39
    assert to_binary({}, names, 1) == '0' * 43


def test_bundled_database():
    """Every entry of the bundled database is valid, and its octal code matches its binary pattern."""
    db = load_database()
    assert len(db) > 1900
    assert lookup(SB0140, db) == 'SB0140'
    assert lookup(H37RV, db) == NOT_FOUND
    with open(SPOLIGOTYPE_DB) as f:
        for line in f:
            if line.strip():
                octal, name, binary = line.split()
                assert binary_to_octal(binary) == octal, name


def test_custom_database(tmp_path):
    db = tmp_path / 'db.txt'
    db.write_text('# Comment\n777777477760771\tSIT451\t{}\n\n'.format(H37RV))  # Tabs work too
    assert load_database(db) == {H37RV: 'SIT451'}


@pytest.mark.parametrize('content, message', [
    ('', 'does not contain'),
    ('664073777777600 SB0140\n', 'expected 3 columns'),
    ('777777777777771 SB0140 {}\n'.format(SB0140), 'does not match'),
    ('664073777777600 SB0140 {0}\n664073777777600 SB9999 {0}\n'.format(SB0140), 'listed as both'),
])
def test_bad_database(tmp_path, content, message):
    db = tmp_path / 'db.txt'
    db.write_text(content)
    with pytest.raises(SpoligoError, match=message):
        load_database(db)
