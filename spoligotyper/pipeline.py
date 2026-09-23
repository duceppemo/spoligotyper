"""Spoligotype one sample: count the spacers, call them present or absent, and name the pattern."""

import gzip
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from . import seal
from .spoligotype import (
    NOT_FOUND,
    SPACERS_FASTA,
    SPOLIGOTYPE_DB,
    SpoligoError,
    binary_to_hex,
    binary_to_octal,
    load_database,
    lookup,
    read_spacer_names,
    to_binary,
)

log = logging.getLogger(__name__)

# Default minimum count for a spacer to be called present. A genome assembly contains each spacer once at most.
MIN_COUNT = {'fastq': 5, 'fasta': 1}

SEQUENCE_EXTENSIONS = ('.fastq', '.fq', '.fasta', '.fa', '.fna', '.fas', '.fsa')
# Paired-end suffixes removed from fastq file names, e.g. "S1_R1", "S1_S12_L001_R1_001", "SRR123_1"
READ_SUFFIX = re.compile(r'_R?[12](_001)?$')

REPORT_HEADER = ['Sample', 'SpacerCount', 'Binary', 'Octal', 'Hexadecimal', 'Spoligotype']


@dataclass
class Result:
    sample: str
    counts: list  # Number of reads (or contigs) matching each spacer, in spacer order
    binary: str
    octal: str
    hexadecimal: str
    spoligotype: str

    @property
    def found(self):
        return self.spoligotype != NOT_FOUND

    def row(self):
        return [self.sample, ':'.join(str(c) for c in self.counts), self.binary, self.octal, self.hexadecimal,
                self.spoligotype]


def open_text(path):
    with open(path, 'rb') as f:
        gzipped = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gzipped else open(path)


def file_type(path):
    """"fasta" or "fastq", from the first character of the (possibly gzipped) file."""
    try:
        with open_text(path) as f:
            for line in f:
                if line.strip():
                    first = line[0]
                    break
            else:
                raise SpoligoError('{} is empty'.format(path))
    except (OSError, EOFError, UnicodeDecodeError) as e:
        raise SpoligoError('Could not read {}: {}'.format(path, e)) from e
    if first == '>':
        return 'fasta'
    if first == '@':
        return 'fastq'
    raise SpoligoError('{} is not a fasta or fastq file'.format(path))


def sample_name(path, kind):
    """
    Sample name from the file name: the sequence extension is removed, and for fastq files the read suffix too.
    "S1_R1.fastq.gz" -> "S1", "Iso_R10.fasta" -> "Iso_R10", "E.coli.v2.fasta" -> "E.coli.v2".
    """
    name = Path(path).name
    if name.lower().endswith('.gz'):
        name = name[:-3]
    for extension in SEQUENCE_EXTENSIONS:
        if name.lower().endswith(extension):
            name = name[:-len(extension)]
            break
    if kind == 'fastq':
        name = READ_SUFFIX.sub('', name) or name
    return name


def check_inputs(r1, r2=None):
    """Validate the input files and return their type ("fasta" or "fastq")."""
    for path in (r1, r2):
        if path is not None and not Path(path).is_file():
            raise SpoligoError('Input file not found: {}'.format(path))
    kind = file_type(r1)
    if r2 is not None:
        if Path(r1).resolve() == Path(r2).resolve():
            raise SpoligoError('-r1 and -r2 are the same file: {}'.format(r1))
        if kind != 'fastq' or file_type(r2) != 'fastq':
            raise SpoligoError('-r2 is only for paired-end fastq files; both -r1 and -r2 must be fastq')
    return kind


def spoligotype(r1, r2=None, sample=None, min_count=None, threads=1, memory='1g', database=SPOLIGOTYPE_DB,
                spacers=SPACERS_FASTA):
    """
    Spoligotype a sample from single-end or paired-end reads, or from an assembly.

    :param min_count: minimum number of reads (or contigs) matching a spacer to call it present.
                      Default: 5 for fastq files, 1 for fasta files.
    :return: Result
    """
    kind = check_inputs(r1, r2)
    if sample is None:
        sample = sample_name(r1, kind)
    if min_count is None:
        min_count = MIN_COUNT[kind]
    elif kind == 'fasta' and min_count > 1:
        log.warning('%s is a fasta file: with --min-count %d, spacers are probably missed. '
                    'Leave --min-count unset (1 for fasta files).', r1, min_count)
    spacer_names = read_spacer_names(spacers)
    db = load_database(database)  # Before running Seal, to report a bad database right away

    log.info('Spoligotyping %s (%s, minimum count %d)', sample, kind, min_count)
    counts = seal.count_spacers([f for f in (r1, r2) if f], spacers, threads=threads, memory=memory)
    binary = to_binary(counts, spacer_names, min_count)
    result = Result(sample, [counts.get(name, 0) for name in spacer_names], binary, binary_to_octal(binary),
                    binary_to_hex(binary), lookup(binary, db))
    check_counts(result, min_count, kind)
    return result


def check_counts(result, min_count, kind):
    """Warn about results that deserve a second look."""
    if not any(result.counts):
        log.warning('No spacer found in %s. Is it a Mycobacterium tuberculosis complex sample?', result.sample)
        return
    borderline = [i + 1 for i, c in enumerate(result.counts) if 0 < c < min_count]
    if kind == 'fastq' and borderline:
        log.warning('%d spacer(s) called absent were seen in fewer than %d reads: %s. '
                    'Low coverage or contamination? See the SpacerCount column.',
                    len(borderline), min_count, ', '.join(str(i) for i in borderline))


def write_report(results, path):
    """Tab-separated report, one line per result."""
    with open(path, 'w') as f:
        f.write('\t'.join(REPORT_HEADER) + '\n')
        for result in results:
            f.write('\t'.join(result.row()) + '\n')
