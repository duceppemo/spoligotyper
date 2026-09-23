import gzip
import hashlib

import pytest

from spoligotyper import seal
from spoligotyper.pipeline import InputFile, Result, check_inputs, file_md5, file_type, write_tsv
from spoligotyper.samples import sample_name
from spoligotyper.spoligotype import SpoligoError


@pytest.mark.parametrize('path, kind, expected', [
    ('S1_R1.fastq.gz', 'fastq', 'S1'),
    ('S1_S12_L001_R1_001.fastq.gz', 'fastq', 'S1_S12_L001'),
    ('/data/SRR123_1.fq', 'fastq', 'SRR123'),
    ('ERR1744454.fastq.gz', 'fastq', 'ERR1744454'),
    ('Iso_R10.fastq', 'fastq', 'Iso_R10'),
    ('Iso_R1.fasta', 'fasta', 'Iso_R1'),
    ('Iso_R10.fasta', 'fasta', 'Iso_R10'),
    ('NC_002945.4.fasta', 'fasta', 'NC_002945.4'),
    ('E.coli.v2.FNA.gz', 'fasta', 'E.coli.v2'),
    ('assembly.contigs', 'fasta', 'assembly.contigs'),
    ('_R1.fastq', 'fastq', '_R1'),
])
def test_sample_name(path, kind, expected):
    assert sample_name(path, kind) == expected


def test_file_type(tmp_path):
    fasta = tmp_path / 'a.fa'
    fasta.write_text('\n>c1\nACGT\n')
    assert file_type(fasta) == 'fasta'
    fastq = tmp_path / 'a.fq.gz'  # Detected from the content, not the extension
    with gzip.open(fastq, 'wt') as f:
        f.write('@r1\nACGT\n+\nIIII\n')
    assert file_type(fastq) == 'fastq'


@pytest.mark.parametrize('content, message', [(b'', 'empty'), (b'\n\n', 'empty'), (b'ACGT\n', 'not a fasta'),
                                              (b'\x1f\x8bnot gzip', 'Could not read'), (b'\xff\xfe>', 'Could not read')])
def test_file_type_errors(tmp_path, content, message):
    path = tmp_path / 'bad'
    path.write_bytes(content)
    with pytest.raises(SpoligoError, match=message):
        file_type(path)


def test_check_inputs(tmp_path):
    fasta = tmp_path / 'a.fasta'
    fasta.write_text('>c1\nACGT\n')
    r1, r2 = tmp_path / 'r_1.fq', tmp_path / 'r_2.fq'
    for r in (r1, r2):
        r.write_text('@r1\nACGT\n+\nIIII\n')
    assert check_inputs(fasta) == 'fasta'
    assert check_inputs(r1, r2) == 'fastq'
    with pytest.raises(SpoligoError, match='not found'):
        check_inputs(tmp_path / 'missing.fq')
    with pytest.raises(SpoligoError, match='not found'):
        check_inputs(r1, tmp_path / 'missing.fq')
    with pytest.raises(SpoligoError, match='same file'):
        check_inputs(r1, r1)
    with pytest.raises(SpoligoError, match='both R1 and R2 must be fastq'):
        check_inputs(fasta, r2)


def test_parse_stats(tmp_path):
    stats = tmp_path / 'stats.tsv'
    stats.write_text('#File\tx.fq\n#Total\t100\t15000\n#Matched\t3\t3%\n#Name\tReads\tReadsPct\n'
                     'spacer25\t2\t2%\nspacer02\t1\t1%\nspacer03\t0\t0%\n')
    parsed = seal.parse_stats(stats)
    assert parsed.counts == {'spacer25': 2, 'spacer02': 1, 'spacer03': 0}
    assert (parsed.reads, parsed.bases) == (100, 15000)
    stats.write_text('spacer25\tmany\n')
    with pytest.raises(seal.SealError, match='Unexpected line'):
        seal.parse_stats(stats)


def test_seal_command(monkeypatch):
    monkeypatch.setattr(seal, 'check_seal', lambda: 'seal.sh')
    cmd = seal.seal_command(['r1.fq', 'r2.fq'], 'spacers.fa', 'stats.tsv', 4, '2g')
    assert cmd[:4] == ['seal.sh', '-Xmx2g', 'in=r1.fq', 'in2=r2.fq']
    assert {'ref=spacers.fa', 'k=25', 'hdist=1', 'nzo=f', 'stats=stats.tsv', 'threads=4'} <= set(cmd)
    assert not any(c.startswith('in2=') for c in seal.seal_command(['a.fa'], 'spacers.fa', 'stats.tsv', 1, '1g'))


def test_seal_missing(monkeypatch):
    monkeypatch.setattr(seal, 'executable', lambda: None)
    with pytest.raises(seal.SealError, match='conda install -c bioconda bbmap'):
        seal.check_seal()


def test_write_tsv(tmp_path):
    ok = Result('S1', counts=[3, 0] + [1] * 41, binary='10' + '1' * 41, octal='x', hexadecimal='y',
                spoligotype='SB0000', file_type='fastq', min_count=1, reads=1000, bases=88_000_000)
    ok.warnings.append('a warning')
    failed = Result('S2', file_type='fasta', error='Seal failed: boom\nmore details')
    write_tsv([ok, failed], tmp_path / 'report.tsv')
    lines = [line.split('\t') for line in (tmp_path / 'report.tsv').read_text().splitlines()]
    assert lines[0] == ['Sample', 'SpacerCount', 'Binary', 'Octal', 'Hexadecimal', 'Spoligotype',
                        'FileType', 'Reads', 'Depth', 'MinCount', 'Status', 'Warnings']
    assert lines[1][:2] == ['S1', '3:0:' + ':'.join(['1'] * 41)]
    assert lines[1][6:] == ['fastq', '1000', '20', '1', 'warning', 'a warning']
    assert lines[2][0] == 'S2' and lines[2][-2:] == ['failed', 'Seal failed: boom']
    assert all(len(line) == 12 for line in lines)


def test_result_properties():
    r = Result('S', counts=[10, 0, 30] + [0] * 40, binary='101' + '0' * 40, file_type='fasta', bases=5)
    assert r.depth is None and r.median_present_count == 20 and r.status == 'ok' and not r.found
    r.file_type = 'fastq'
    assert r.depth == 5 / 4.4e6


def test_input_file(tmp_path):
    real = tmp_path / 'real.fq'
    real.write_text('@r1\nACGT\n+\nIIII\n')
    link = tmp_path / 'link.fq'
    link.symlink_to(real)
    info = InputFile.describe(link)
    assert info.path == str(link) and info.target == str(real)
    assert info.size == real.stat().st_size
    assert info.md5 == file_md5(real) == hashlib.md5(b'@r1\nACGT\n+\nIIII\n').hexdigest()
    assert InputFile.describe(real, md5=False).md5 == '' and InputFile.describe(real).target == ''


def test_seal_versions(monkeypatch):
    monkeypatch.setattr(seal, 'executable', lambda: None)
    assert seal.versions() == {'BBTools': 'unknown', 'Java': 'unknown'}
