"""End-to-end tests on synthetic data. Skipped when Seal (BBTools) is not installed."""

import logging
import os
import subprocess
import sys

import pytest

from spoligotyper import __version__, seal
from spoligotyper.cli import main
from spoligotyper.pipeline import spoligotype

from .conftest import H37RV, SB0140

pytestmark = pytest.mark.skipif(seal.executable() is None, reason='seal.sh (BBTools) is not installed')


def run(capsys, *args):
    main([*map(str, args), '-t', '2', '--memory', '500m', '--no-pdf'])
    out = capsys.readouterr().out.splitlines()
    assert out[0].startswith('Sample\t')
    return out[1].split('\t')


def test_assembly(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'AF2122.fasta', '-o', tmp_path)
    assert row[0] == 'AF2122'
    assert row[2:10] == [SB0140, '664073777777600', '6D-03-5F-7F-FF-60', 'SB0140', 'fasta', '1', '', '1']
    assert row[10:] == ['ok', '']
    assert row[1].split(':')[:4] == ['1', '1', '0', '1']
    report = (tmp_path / 'AF2122_spoligotyping.txt').read_text().splitlines()
    assert report[1].split('\t') == row


def test_not_in_database(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'H37Rv.fna', '-o', tmp_path)
    assert row[2:6] == [H37RV, '777777477760771', '7F-7F-7C-7F-F0-7F', 'Spoligo not found']


def test_paired_end(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'bovis_R1.fastq.gz', '-r2', data / 'bovis_R2.fastq.gz', '-o', tmp_path)
    assert row[0] == 'bovis'
    assert row[2] == SB0140 and row[5] == 'SB0140'
    counts = [int(c) for c in row[1].split(':')]
    assert min(c for c, bit in zip(counts, SB0140, strict=True) if bit == '1') >= 5
    assert max(c for c, bit in zip(counts, SB0140, strict=True) if bit == '0') == 0


def test_single_end_and_sample_name(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'bovis_single.fq.gz', '-o', tmp_path, '-s', 'my_sample')
    assert row[0] == 'my_sample' and row[5] == 'SB0140'
    assert (tmp_path / 'my_sample_spoligotyping.txt').exists()


def test_min_count(data, tmp_path, capsys, caplog):
    """With too few reads, spacers are called absent and a warning points to the low counts."""
    with caplog.at_level(logging.WARNING, logger='spoligotyper'):
        row = run(capsys, '-r1', data / 'low_coverage.fastq.gz', '-o', tmp_path)
    assert row[2] != SB0140
    assert 'fewer than 5 reads' in caplog.text
    row = run(capsys, '-r1', data / 'low_coverage.fastq.gz', '-o', tmp_path, '-m', '1')
    assert row[5] == 'SB0140'


def test_not_mtbc(data, tmp_path, capsys, caplog):
    with caplog.at_level(logging.WARNING, logger='spoligotyper'):
        row = run(capsys, '-r1', data / 'not_mtbc.fasta', '-o', tmp_path)
    assert row[2] == '0' * 43
    assert 'no spacer found' in caplog.text
    assert row[5] == 'SB2277' and row[10] == 'warning'  # SB2277 is the pattern with no spacer


def test_fasta_min_count_warning(data, caplog):
    with caplog.at_level(logging.WARNING, logger='spoligotyper'):
        result = spoligotype(data / 'AF2122.fasta', min_count=5, threads=2, memory='500m')
    assert result.binary == '0' * 43
    assert '--min-count' in caplog.text
    assert result.status == 'warning' and any('--min-count' in w for w in result.warnings)


def test_seal_failure(data, tmp_path, caplog):
    with pytest.raises(SystemExit) as e:
        main(['-r1', str(data / 'AF2122.fasta'), '-o', str(tmp_path), '--memory', '1m'])
    assert e.value.code == 1
    assert 'Seal failed' in caplog.text
    assert not list(tmp_path.iterdir())


def test_user_errors(tmp_path, caplog):
    with pytest.raises(SystemExit) as e:
        main(['-r1', str(tmp_path / 'missing.fastq'), '-o', str(tmp_path / 'out')])
    assert e.value.code == 1
    assert 'Input file not found' in caplog.text
    assert not (tmp_path / 'out').exists()


@pytest.mark.parametrize('args', [['-m', '0'], ['--memory', '2'], ['-s', '../x'], ['-t', 'two']])
def test_bad_arguments(data, tmp_path, args):
    with pytest.raises(SystemExit) as e:
        main(['-r1', str(data / 'AF2122.fasta'), '-o', str(tmp_path), *args])
    assert e.value.code == 2


def test_module_entry_point():
    proc = subprocess.run([sys.executable, '-m', 'spoligotyper', '--version'], capture_output=True, text=True)
    assert proc.stdout.strip() == 'spoligotyper ' + __version__


def pdf_text(path):
    pypdf = pytest.importorskip('pypdf')
    reader = pypdf.PdfReader(path)
    return len(reader.pages), '\n'.join(page.extract_text() for page in reader.pages)


def test_pdf_single_sample(data, tmp_path, capsys):
    main(['-r1', str(data / 'bovis_R1.fastq.gz'), '-r2', str(data / 'bovis_R2.fastq.gz'), '-o', str(tmp_path),
          '-t', '2', '--memory', '500m', '--operator', 'Jane Doe'])
    capsys.readouterr()
    pdf = tmp_path / 'bovis_spoligotyping.pdf'
    assert (tmp_path / 'bovis_spoligotyping.txt').exists()
    pages, content = pdf_text(pdf)
    assert pages >= 2
    for expected in ('Spoligotyping report', 'bovis', 'SB0140', '664073777777600', 'paired-end', 'Jane Doe',
                     'Run information', 'Reads per spacer', 'BBTools', __version__, 'MD5'):
        assert expected in content, expected


def test_batch(data, tmp_path, capsys):
    """A folder of assemblies, single-end and paired-end reads, with a broken file: the others are still typed."""
    folder = tmp_path / 'input'
    (folder / 'reads').mkdir(parents=True)
    for f in data.iterdir():
        target = folder / ('reads' if '.f' in f.name and 'q' in f.suffixes[0] else '') / f.name
        os.symlink(f, target)
    (folder / 'broken.fasta').write_text('not a sequence\n')
    out = tmp_path / 'out'
    with pytest.raises(SystemExit) as e:
        main(['-i', str(folder), '-o', str(out), '-t', '2', '--memory', '500m', '--no-md5'])
    assert e.value.code == 1  # One sample failed
    rows = {line.split('\t')[0]: line.split('\t') for line in capsys.readouterr().out.splitlines()[1:]}
    assert sorted(rows) == ['AF2122', 'H37Rv', 'bovis', 'bovis_single', 'broken', 'low_coverage', 'not_mtbc']
    assert rows['bovis'][5] == rows['bovis_single'][5] == rows['AF2122'][5] == 'SB0140'
    assert rows['H37Rv'][3] == '777777477760771'
    assert rows['broken'][10] == 'failed' and 'not a fasta or fastq' in rows['broken'][11]
    assert rows['low_coverage'][10] == rows['not_mtbc'][10] == 'warning'
    tsv = (out / 'spoligotyping.tsv').read_text().splitlines()
    assert len(tsv) == 8
    pages, content = pdf_text(out / 'spoligotyping_report.pdf')
    assert all(name in content for name in rows) and 'FAILED' in content


def test_batch_arguments(data, tmp_path):
    for extra in (['-r2', 'x.fq'], ['-s', 'name']):
        with pytest.raises(SystemExit) as e:
            main(['-i', str(data), '-o', str(tmp_path), *extra])
        assert e.value.code == 2
    with pytest.raises(SystemExit) as e:
        main(['-i', str(data), '-r1', str(data / 'AF2122.fasta'), '-o', str(tmp_path)])
    assert e.value.code == 2
