import pytest

from spoligotyper.samples import find_samples, split_extension
from spoligotyper.spoligotype import SpoligoError


def touch(folder, *names):
    for name in names:
        path = folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('x')


def test_split_extension():
    assert split_extension('S1_R1.fastq.gz') == ('S1_R1', 'fastq')
    assert split_extension('a.b.FNA') == ('a.b', 'fasta')
    assert split_extension('notes.txt') is None
    assert split_extension('.fastq') is None


def test_find_samples(tmp_path):
    touch(tmp_path, 'run1/S1_R1.fastq.gz', 'run1/S1_R2.fastq.gz',  # Paired
          'run1/S2_S7_L001_R1_001.fastq.gz', 'run1/S2_S7_L001_R2_001.fastq.gz',  # Illumina names
          'run2/SRR1_1.fq', 'run2/SRR1_2.fq',  # SRA names
          'run2/single.fastq',  # Single-end
          'run2/lonely_R2.fastq',  # R2 without R1: single-end
          'assemblies/Iso_R1.fasta', 'assemblies/Iso_R10.fna.gz',  # Fasta: no read suffix
          'notes.txt', '.hidden.fasta', '.snapshots/old.fasta')
    samples = {s.name: s for s in find_samples(tmp_path)}
    assert sorted(samples) == ['Iso_R1', 'Iso_R10', 'S1', 'S2_S7_L001', 'SRR1', 'lonely', 'single']
    assert [p.split('/')[-1] for p in samples['S1'].files] == ['S1_R1.fastq.gz', 'S1_R2.fastq.gz']
    assert samples['S1'].paired and samples['SRR1'].paired and samples['S2_S7_L001'].paired
    assert not samples['single'].paired and not samples['lonely'].paired
    assert samples['Iso_R1'].file_type == 'fasta'


def test_output_folder_excluded(tmp_path):
    touch(tmp_path, 'S1.fasta', 'results/S2.fasta')
    assert [s.name for s in find_samples(tmp_path, exclude=[tmp_path / 'results'])] == ['S1']


@pytest.mark.parametrize('names', [
    ['a/S1.fasta', 'b/S1.fasta'],  # Same name in two folders
    ['S1.fasta', 'S1_R1.fastq', 'S1_R2.fastq'],  # Fasta and fastq
    ['a/S1_R1.fastq', 'b/S1_R1.fastq'],  # Two R1
    ['S1.fastq', 'S1_R1.fastq'],  # Single-end and R1
])
def test_ambiguous_names(tmp_path, names):
    touch(tmp_path, *names)
    with pytest.raises(SpoligoError, match='same sample name'):
        find_samples(tmp_path)


def test_no_samples(tmp_path):
    touch(tmp_path, 'notes.txt')
    with pytest.raises(SpoligoError, match='No fasta or fastq'):
        find_samples(tmp_path)
    with pytest.raises(SpoligoError, match='not found'):
        find_samples(tmp_path / 'missing')


def test_symlinked_folders(tmp_path):
    touch(tmp_path, 'real/S1.fasta', 'input/S2.fasta')
    (tmp_path / 'input' / 'linked').symlink_to(tmp_path / 'real')
    (tmp_path / 'input' / 'loop').symlink_to(tmp_path / 'input')  # Must not loop forever
    assert [s.name for s in find_samples(tmp_path / 'input')] == ['S1', 'S2']
