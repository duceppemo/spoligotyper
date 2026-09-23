"""Find the samples in an input folder: fasta files, and single-end or paired-end fastq files."""

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path

from .spoligotype import SpoligoError

log = logging.getLogger(__name__)

FASTA_EXTENSIONS = ('.fasta', '.fa', '.fna', '.fas', '.fsa')
FASTQ_EXTENSIONS = ('.fastq', '.fq')
# Read suffixes of fastq file names, e.g. "S1_R1", "S1_S12_L001_R2_001", "SRR123_1". Group 1 is the mate number.
READ_SUFFIX = re.compile(r'_R?([12])(?:_001)?$')


@dataclass
class Sample:
    name: str
    file_type: str  # "fasta" or "fastq"
    files: list  # One file, or R1 and R2

    @property
    def paired(self):
        return len(self.files) == 2


def split_extension(filename):
    """(stem, "fasta" or "fastq") for a sequence file, gzipped or not, or None for any other file."""
    lower = filename.lower()
    if lower.endswith('.gz'):
        lower = lower[:-3]
    for file_type, extensions in (('fasta', FASTA_EXTENSIONS), ('fastq', FASTQ_EXTENSIONS)):
        for extension in extensions:
            if lower.endswith(extension) and len(lower) > len(extension):
                return filename[:len(lower) - len(extension)], file_type
    return None


def sample_name(path, file_type=None):
    """
    Sample name from a file name: the sequence extension is removed, and for fastq files the read suffix too.
    "S1_R1.fastq.gz" -> "S1", "Iso_R10.fasta" -> "Iso_R10", "E.coli.v2.fasta" -> "E.coli.v2".
    """
    name = Path(path).name
    parsed = split_extension(name)
    if parsed is not None:
        name, file_type = parsed[0], file_type or parsed[1]
    if file_type == 'fastq':
        name = READ_SUFFIX.sub('', name) or name
    return name


def find_samples(folder, exclude=()):
    """
    Look for fasta and fastq files in a folder and its subfolders, and group R1 and R2 fastq files per sample.

    Name collisions are errors rather than guesses: two files giving the same sample name would otherwise be
    merged or overwritten silently.

    :param exclude: folders to skip, e.g. the output folder when it is inside the input folder
    :return: list of Sample, sorted by name
    """
    if not Path(folder).is_dir():
        raise SpoligoError('Input folder not found: {}'.format(folder))
    exclude = {Path(p).resolve() for p in exclude}
    found = {}  # {sample name: [(file type, mate number or None, path)]}
    for root, dirs, filenames in os.walk(folder):
        dirs[:] = sorted(d for d in dirs if Path(root, d).resolve() not in exclude and not d.startswith('.'))
        for filename in sorted(filenames):
            parsed = split_extension(filename)
            if parsed is None or filename.startswith('.'):
                continue
            stem, file_type = parsed
            mate = None
            if file_type == 'fastq':
                match = READ_SUFFIX.search(stem)
                if match and match.start() > 0:
                    mate, stem = int(match.group(1)), stem[:match.start()]
            found.setdefault(stem, []).append((file_type, mate, os.path.join(root, filename)))

    samples, errors = [], []
    for name, files in sorted(found.items()):
        paths = ', '.join(f[2] for f in files)
        if len(files) == 1:
            file_type, mate, path = files[0]
            if mate == 2:
                log.warning('%s: R2 file without R1, typed as single-end reads: %s', name, path)
            samples.append(Sample(name, file_type, [path]))
        elif len(files) == 2 and {(f[0], f[1]) for f in files} == {('fastq', 1), ('fastq', 2)}:
            samples.append(Sample(name, 'fastq', [f[2] for f in sorted(files, key=lambda f: f[1])]))
        else:
            errors.append('"{}": {}'.format(name, paths))
    if errors:
        raise SpoligoError('Several files give the same sample name. Rename, move or remove them:\n  '
                           + '\n  '.join(errors))
    if not samples:
        raise SpoligoError('No fasta or fastq file found in {}'.format(folder))
    return samples
