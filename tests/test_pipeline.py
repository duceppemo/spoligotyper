import gzip
import hashlib
from pathlib import Path

import pytest

from spoligotyper import pipeline, seal, species
from spoligotyper.pipeline import InputFile, Result, check_inputs, check_result, file_md5, file_type, write_tsv
from spoligotyper.samples import Sample, sample_name
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
                spoligotype='SB0000', file_type='fastq', data='reads', min_count=1, reads=1000, bases=88_000_000)
    ok.warnings.append('a warning')
    failed = Result('S2', file_type='fasta', error='Seal failed: boom\nmore details')
    write_tsv([ok, failed], tmp_path / 'report.tsv')
    lines = [line.split('\t') for line in (tmp_path / 'report.tsv').read_text().splitlines()]
    assert lines[0] == ['Sample', 'SpacerCount', 'Binary', 'Octal', 'Hexadecimal', 'Spoligotype',
                        'FileType', 'Reads', 'Depth', 'MinCount', 'Status', 'Warnings',
                        'Species', 'Lineage', 'LineageName', 'RD9', 'RD4', 'RD1', 'MTBCFraction', 'Closest',
                        'RD7', 'RD12', 'SIT', 'SITVIT2family', 'ClosestSIT']
    assert lines[1][:2] == ['S1', '3:0:' + ':'.join(['1'] * 41)]
    assert lines[1][6:12] == ['fastq', '1000', '20', '1', 'warning', 'a warning']
    assert lines[2][0] == 'S2' and lines[2][10:12] == ['failed', 'Seal failed: boom']
    assert all(len(line) == 25 for line in lines)


def test_result_properties():
    r = Result('S', counts=[10, 0, 30] + [0] * 40, binary='101' + '0' * 40, file_type='fasta', data='assembly',
               bases=5)
    assert r.depth is None and r.median_present_count == 20 and r.status == 'ok' and not r.found
    assert r.unit == 'contigs'
    r.data = 'reads'
    assert r.depth == r.mtbc_depth == 5 / 4.4e6 and r.unit == 'reads'


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


@pytest.mark.parametrize('lines, expected', [
    (['Input is being processed as unpaired', "Error: truncated or corrupt input for 'x.fq.gz'; data may be "
      'incomplete.', 'java.lang.Exception: ', 'Mismatch between length', '\tat jgi.Seal.main(Seal.java:74)',
      'Exception in thread "main" java.lang.RuntimeException: Seal terminated in an error state; the output may be '
      'corrupt.'], "Error: truncated or corrupt input for 'x.fq.gz'; data may be incomplete."),
    (['Input is being processed as paired', 'java.lang.AssertionError: ',
      'There appear to be different numbers of reads in the paired input files.', '\tat stream.X'],
     'java.lang.AssertionError: There appear to be different numbers of reads in the paired input files.'),
    (['Exception in thread "main" java.lang.RuntimeException: Seal terminated in an error state'],
     'Exception in thread "main" java.lang.RuntimeException: Seal terminated in an error state'),
    (['some output', 'last line'], 'last line'),
])
def test_failure_reason(lines, expected):
    assert seal.failure_reason(lines) == expected


def test_safe_paths(tmp_path):
    odd = tmp_path / 'my sample,1=x.fastq.gz'
    odd.write_text('x')
    plain = tmp_path / 'plain.fasta'
    plain.write_text('>a')
    links = tmp_path / 'links'
    links.mkdir()
    safe = seal.safe_paths([odd, plain], links, 'in')
    assert safe[1] == str(plain)
    assert safe[0] == str(links / 'in0.fastq.gz') and Path(safe[0]).resolve() == odd.resolve()
    # BBTools 40 takes any argument containing "xmx" or "xms" for a Java memory setting
    assert seal.UNSAFE.search('/tmp/spoligotyper_87mxmx1_/stats.tsv') and seal.UNSAFE.search('/d/S_XMS.fq')
    assert not seal.UNSAFE.search('/data/run1/S1_R1.fastq.gz')


@pytest.mark.parametrize('name, content, gz, expected', [
    ('dataset_1.dat', '@r\nACGT\n+\nIIII\n', True, 'in0.fastq.gz'),  # Galaxy names its files .dat
    ('dataset_2.dat', '>c\nACGT\n', False, 'in0.fasta'),
    ('reads.fastq', '@r\nACGT\n+\nIIII\n', True, 'in0.fastq.gz'),  # Gzipped without .gz
    ('reads.fq.gz', '@r\nACGT\n+\nIIII\n', True, None),  # Right extension: used as is
])
def test_extension_from_content(tmp_path, name, content, gz, expected):
    path = tmp_path / name
    if gz:
        with gzip.open(path, 'wt') as f:
            f.write(content)
    else:
        path.write_text(content)
    links = tmp_path / 'links'
    links.mkdir()
    safe = seal.safe_paths([path], links, 'in', check_extension=True)[0]
    assert safe == (str(links / expected) if expected else str(path))


def test_safe_temporary_folder(monkeypatch, tmp_path):
    names = iter(['spoligotyper_87mxmx1_', 'spoligotyper_ok'])

    def mkdtemp(prefix):
        folder = tmp_path / next(names)
        folder.mkdir()
        return str(folder)
    monkeypatch.setattr(seal.tempfile, 'mkdtemp', mkdtemp)
    assert seal.safe_temporary_folder() == str(tmp_path / 'spoligotyper_ok')
    assert not (tmp_path / 'spoligotyper_87mxmx1_').exists()


def test_single_input_not_interleaved(monkeypatch):
    monkeypatch.setattr(seal, 'check_seal', lambda: 'seal.sh')
    assert 'int=f' in seal.seal_command(['r.fq'], 'ref.fa', 'stats.tsv', 1, '1g')
    assert 'int=f' not in seal.seal_command(['r1.fq', 'r2.fq'], 'ref.fa', 'stats.tsv', 1, '1g')


def test_low_mtbc_depth_warning():
    """38x in total but 40% MTBC reads: 15x of MTBC, below the 20x threshold."""
    r = Result('S', counts=[20] * 43, binary='1' * 43, file_type='fastq', data='reads', min_count=5,
               reads=1_000_000, bases=int(38 * 4.4e6), species=species.SpeciesCheck(mtbc=True, mtbc_fraction=0.4))
    check_result(r)
    assert any('estimated MTBC depth 15.2x (38.0x in total)' in w for w in r.warnings)
    r = Result('S', counts=[20] * 43, binary='1' * 43, file_type='fastq', data='reads', min_count=5,
               reads=1_000_000, bases=int(19.6 * 4.4e6))
    check_result(r)
    assert any('estimated depth 19.6x' in w for w in r.warnings)


def test_no_sb_warning_without_spacers():
    check = species.SpeciesCheck(mtbc=True, regions={'RD9': ('present', 1.0)}, species='x')
    r = Result('S', counts=[0] * 43, binary='0' * 43, spoligotype='SB2277', file_type='fasta', data='assembly',
               min_count=1, species=check, lineage=pipeline.lineage.LineageCall())
    check_result(r)
    assert not any('SB number' in w for w in r.warnings) and any('no spacer found' in w for w in r.warnings)


def test_batch_unexpected_error(monkeypatch, tmp_path):
    """A bug in one sample must not stop the batch."""
    def fake(*files, sample=None, **kwargs):
        if sample == 'bad':
            raise KeyError('boom')
        return Result(sample, file_type='fasta', data='assembly')
    monkeypatch.setattr(pipeline, 'spoligotype', fake)
    samples = [Sample('good', 'fasta', [str(tmp_path / 'good.fasta')]),
               Sample('bad', 'fasta', [str(tmp_path / 'bad.fasta')])]
    results = pipeline.spoligotype_samples(samples, jobs=2)
    assert [r.sample for r in results] == ['good', 'bad']
    assert results[0].status == 'ok' and results[1].status == 'failed' and 'KeyError' in results[1].error


def test_pdf_order_by_spoligotype():
    from spoligotyper.pdf import by_spoligotype
    a, b = '1' * 43, '0' * 43
    results = [Result('S3', binary=b, spoligotype='Spoligo not found', octal='000'),
               Result('S2', binary=a, spoligotype='SB0001', octal='777'),
               Result('bad', error='boom'),
               Result('S1', binary=a, spoligotype='SB0001', octal='777'),
               Result('S0', binary='01' * 21 + '0', spoligotype='SB0002', octal='252')]
    order = [(group, r.sample) for group, r in by_spoligotype(results)]
    # Largest group first, then SB numbers before unnamed patterns, samples sorted within a group, failed last
    assert order == [(0, 'S1'), (0, 'S2'), (1, 'S0'), (2, 'S3'), (3, 'bad')]
