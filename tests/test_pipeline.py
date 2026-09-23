import gzip

import pytest

from spoligotyper import seal
from spoligotyper.pipeline import Result, check_inputs, file_type, sample_name, write_report
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
    with pytest.raises(SpoligoError, match='both -r1 and -r2 must be fastq'):
        check_inputs(fasta, r2)


def test_parse_stats(tmp_path):
    stats = tmp_path / 'stats.tsv'
    stats.write_text('#File\tx.fq\n#Total\t100\n#Matched\t3\t3%\n#Name\tReads\tReadsPct\n'
                     'spacer25\t2\t2%\nspacer02\t1\t1%\nspacer03\t0\t0%\n')
    assert seal.parse_stats(stats) == {'spacer25': 2, 'spacer02': 1, 'spacer03': 0}
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


def test_write_report(tmp_path):
    result = Result('S1', [3, 0] + [1] * 41, '10' + '1' * 41, 'x', 'y', 'SB0000')
    write_report([result], tmp_path / 'report.txt')
    lines = (tmp_path / 'report.txt').read_text().splitlines()
    assert lines[0] == 'Sample\tSpacerCount\tBinary\tOctal\tHexadecimal\tSpoligotype'
    assert lines[1].split('\t')[:2] == ['S1', '3:0:' + ':'.join(['1'] * 41)]
