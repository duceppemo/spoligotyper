"""End-to-end tests on synthetic data. Skipped when Seal (BBTools) is not installed."""

import gzip
import json
import logging
import os
import subprocess
import sys

import pytest

from spoligotyper import __version__, seal
from spoligotyper.cli import main
from spoligotyper.pipeline import spoligotype

from .conftest import H37RV, SB0120, SB0140

pytestmark = pytest.mark.skipif(seal.executable() is None, reason='seal.sh (BBTools) is not installed')


def table(output):
    """Rows of the printed table, as {column: value}."""
    lines = [line.split('\t') for line in output.splitlines()]
    assert lines[0][:2] == ['Sample', 'SpacerCount']
    return [dict(zip(lines[0], line, strict=True)) for line in lines[1:]]


def run(capsys, *args):
    main([*map(str, args), '-t', '2', '--memory', '500m', '--no-pdf'])
    return table(capsys.readouterr().out)[0]


def test_assembly(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'AF2122.fasta', '-o', tmp_path)
    assert list(row.values())[:12] == ['AF2122', row['SpacerCount'], SB0140, '664073777777600', '6D-03-5F-7F-FF-60',
                                       'SB0140', 'fasta', '1', '', '1', 'ok', '']  # The columns of version 0.3
    assert row['SpacerCount'].split(':')[:4] == ['1', '1', '0', '1']
    assert (row['Species'], row['Lineage'], row['LineageName']) == ('M. bovis', 'BOV', 'M. bovis')
    assert (row['RD9'], row['RD4'], row['RD1'], row['MTBCFraction'], row['Closest']) == \
        ('deleted', 'deleted', 'present', '', '')
    report = table((tmp_path / 'AF2122_spoligotyping.txt').read_text())[0]
    assert report == row
    data_json = json.loads((tmp_path / 'AF2122_spoligotyping.json').read_text())
    sample = data_json['samples'][0]
    assert sample['spoligotype'] == 'SB0140' and sample['species']['species'] == 'M. bovis'
    assert sample['species']['regions']['RD4']['state'] == 'deleted' and sample['lineage']['lineage'] == 'BOV'
    assert data_json['run']['software']['spoligotyper'] == __version__
    mqc = json.loads((tmp_path / 'AF2122_spoligotyping_mqc.json').read_text())
    assert mqc['data']['AF2122'] == {'Spoligotype': 'SB0140', 'Octal': '664073777777600', 'Species': 'M. bovis',
                                     'Lineage': 'BOV', 'Status': 'ok'}


def test_not_in_database(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'H37Rv.fna', '-o', tmp_path)
    assert [row[c] for c in ('Binary', 'Octal', 'Hexadecimal', 'Spoligotype')] == \
        [H37RV, '777777477760771', '7F-7F-7C-7F-F0-7F', 'Spoligo not found']
    assert (row['Species'], row['Lineage'], row['RD9']) == ('M. tuberculosis', '4.9', 'present')
    assert row['Status'] == 'ok' and row['Closest'] == ''  # Nothing within 3 spacers


def test_bcg(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'BCG.fasta', '-o', tmp_path)
    assert (row['Binary'], row['Spoligotype'], row['Species'], row['RD1']) == \
        (SB0120, 'SB0120', 'M. bovis BCG', 'deleted')


def test_no_species(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'AF2122.fasta', '-o', tmp_path, '--no-species')
    assert row['Spoligotype'] == 'SB0140' and row['Species'] == row['Lineage'] == row['RD9'] == ''


def test_mixed_sample(data, tmp_path, capsys, caplog):
    """H37Rv reads with a third of M. bovis reads."""
    with caplog.at_level(logging.WARNING, logger='spoligotyper'):
        row = run(capsys, '-r1', data / 'mixed.fastq.gz', '-o', tmp_path)
    assert row['Species'] == 'MTBC, mixed sample?' and row['Lineage'].startswith('mixed: ')
    assert 'BOV' in row['Lineage'] and '4.9' in row['Lineage']
    assert 'mixed sample?' in caplog.text and row['Status'] == 'warning'


def test_paired_end(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'bovis_R1.fastq.gz', '-r2', data / 'bovis_R2.fastq.gz', '-o', tmp_path)
    assert row['Sample'] == 'bovis'
    assert row['Binary'] == SB0140 and row['Spoligotype'] == 'SB0140'
    assert (row['Species'], row['Lineage']) == ('M. bovis', 'BOV')
    counts = [int(c) for c in row['SpacerCount'].split(':')]
    assert min(c for c, bit in zip(counts, SB0140, strict=True) if bit == '1') >= 5
    assert max(c for c, bit in zip(counts, SB0140, strict=True) if bit == '0') == 0


def test_single_end_and_sample_name(data, tmp_path, capsys):
    row = run(capsys, '-r1', data / 'bovis_single.fq.gz', '-o', tmp_path, '-s', 'my_sample')
    assert row['Sample'] == 'my_sample' and row['Spoligotype'] == 'SB0140'
    assert (tmp_path / 'my_sample_spoligotyping.txt').exists()


def test_min_count(data, tmp_path, capsys, caplog):
    """With too few reads, spacers are called absent and a warning points to the low counts."""
    with caplog.at_level(logging.WARNING, logger='spoligotyper'):
        row = run(capsys, '-r1', data / 'low_coverage.fastq.gz', '-o', tmp_path)
    assert row['Binary'] != SB0140
    assert 'fewer than 5 reads' in caplog.text
    row = run(capsys, '-r1', data / 'low_coverage.fastq.gz', '-o', tmp_path, '-m', '1')
    assert row['Spoligotype'] == 'SB0140'


def test_not_mtbc(data, tmp_path, capsys, caplog):
    with caplog.at_level(logging.WARNING, logger='spoligotyper'):
        row = run(capsys, '-r1', data / 'not_mtbc.fasta', '-o', tmp_path)
    assert row['Binary'] == '0' * 43
    assert 'no spacer found' in caplog.text
    assert row['Spoligotype'] == 'SB2277' and row['Status'] == 'warning'  # SB2277 is the pattern with no spacer
    assert row['Species'] == 'MTBC not detected' and row['Lineage'] == ''


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
                     'Run information', 'Reads per spacer', 'BBTools', __version__, 'MD5', 'Species and lineage',
                     'M. bovis', 'RD4', 'Coll F et al.', 'Lineage SNP'):
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
        main(['-i', str(folder), '-o', str(out), '-t', '4', '-j', '3', '--memory', '500m', '--no-md5'])
    assert e.value.code == 1  # One sample failed
    rows = {row['Sample']: row for row in table(capsys.readouterr().out)}
    assert list(rows) == ['AF2122', 'BCG', 'H37Rv', 'bovis', 'bovis_single', 'broken', 'low_coverage', 'mixed',
                          'not_mtbc']  # In order, although typed in parallel
    assert rows['bovis']['Spoligotype'] == rows['bovis_single']['Spoligotype'] == rows['AF2122']['Spoligotype'] \
        == 'SB0140'
    assert rows['H37Rv']['Octal'] == '777777477760771'
    assert rows['broken']['Status'] == 'failed' and 'not a fasta or fastq' in rows['broken']['Warnings']
    assert rows['low_coverage']['Status'] == rows['not_mtbc']['Status'] == 'warning'
    tsv = (out / 'spoligotyping.tsv').read_text().splitlines()
    assert len(tsv) == 10
    assert len(json.loads((out / 'spoligotyping.json').read_text())['samples']) == 9
    assert (out / 'spoligotyping_mqc.json').exists()
    pages, content = pdf_text(out / 'spoligotyping_report.pdf')
    assert all(name in content for name in rows) and 'FAILED' in content


def test_batch_arguments(data, tmp_path):
    for extra in (['-r2', 'x.fq'], ['-s', 'name']):
        with pytest.raises(SystemExit) as e:
            main(['-i', str(data), '-o', str(tmp_path), *extra])
        assert e.value.code == 2
    for args in (['-i', str(data), '-r1', str(data / 'AF2122.fasta')], ['-r1', str(data / 'AF2122.fasta'), '-j', '2']):
        with pytest.raises(SystemExit) as e:
            main([*args, '-o', str(tmp_path)])
        assert e.value.code == 2


def test_odd_file_names(data, tmp_path, capsys):
    """seal.sh splits its arguments on spaces: such paths are linked under safe names."""
    folder = tmp_path / 'my folder'
    folder.mkdir()
    odd = folder / 'my sample,1.fasta'
    odd.write_text((data / 'AF2122.fasta').read_text())
    row = run(capsys, '-r1', odd, '-o', tmp_path / 'out put')
    assert row['Sample'] == 'my sample,1' and row['Spoligotype'] == 'SB0140' and row['Species'] == 'M. bovis'
    assert (tmp_path / 'out put' / 'my sample,1_spoligotyping.txt').exists()


def test_reads_in_fasta_format(data, tmp_path, capsys):
    """Reads converted to fasta are typed as reads (minimum count 5), not as an assembly."""
    fasta = tmp_path / 'reads.fasta'
    with gzip.open(data / 'bovis_single.fq.gz', 'rt') as f, open(fasta, 'w') as out:
        lines = f.read().splitlines()
        for i in range(0, len(lines), 4):
            out.write('>{}\n{}\n'.format(lines[i][1:], lines[i + 1]))
    monkey_threshold = 'spoligotyper.pipeline.FASTA_READS_MIN_BASES'  # The synthetic reads are small
    import unittest.mock
    with unittest.mock.patch(monkey_threshold, 100_000):
        row = run(capsys, '-r1', fasta, '-o', tmp_path)
    assert row['FileType'] == 'fasta' and row['MinCount'] == '5' and row['Spoligotype'] == 'SB0140'
    assert 'typed as reads in fasta format' in row['Warnings']


def test_interleaved_fastq_is_not_paired(data, tmp_path, capsys):
    """A single fastq file is typed as single-end reads even when its reads are interleaved pairs."""
    with gzip.open(data / 'bovis_R1.fastq.gz', 'rt') as f1, gzip.open(data / 'bovis_R2.fastq.gz', 'rt') as f2:
        r1, r2 = f1.read().splitlines(), f2.read().splitlines()
    interleaved = tmp_path / 'interleaved.fastq'
    with open(interleaved, 'w') as out:
        for i in range(0, len(r1), 4):
            out.write('\n'.join(r1[i:i + 4] + r2[i:i + 4]) + '\n')
    counts = {}
    for name in ('R1', 'R2'):
        counts[name] = run(capsys, '-r1', data / 'bovis_{}.fastq.gz'.format(name), '-o', tmp_path, '-s', name)
    row = run(capsys, '-r1', interleaved, '-o', tmp_path)
    total = [int(a) + int(b) for a, b in zip(counts['R1']['SpacerCount'].split(':'),
                                               counts['R2']['SpacerCount'].split(':'), strict=True)]
    assert [int(c) for c in row['SpacerCount'].split(':')] == total


def test_md5_without_pdf(data, tmp_path, capsys):
    run(capsys, '-r1', data / 'AF2122.fasta', '-o', tmp_path)  # run() uses --no-pdf
    sample = json.loads((tmp_path / 'AF2122_spoligotyping.json').read_text())['samples'][0]
    assert len(sample['files'][0]['md5']) == 32


def test_package_path_with_space(tmp_path, data):
    """The package data (spacers, markers) is found and usable from a folder with a space."""
    import shutil

    import spoligotyper
    target = tmp_path / 'pkg dir'
    shutil.copytree(os.path.dirname(spoligotyper.__file__), target / 'spoligotyper')
    env = dict(os.environ, PYTHONPATH=str(target))
    proc = subprocess.run([sys.executable, '-m', 'spoligotyper', '-r1', str(data / 'AF2122.fasta'), '-o',
                           str(tmp_path / 'out'), '-t', '2', '--memory', '500m', '--no-pdf'],
                          capture_output=True, text=True, env=env, cwd=tmp_path)
    assert proc.returncode == 0, proc.stderr
    assert 'SB0140' in proc.stdout and str(target) in (tmp_path / 'out' / 'AF2122_spoligotyping.json').read_text()
