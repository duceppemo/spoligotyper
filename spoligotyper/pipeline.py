"""Spoligotype samples: count the spacers, call them present or absent, and name the pattern."""

import getpass
import gzip
import hashlib
import logging
import os
import platform
import socket
import statistics
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from . import __version__, seal
from .samples import sample_name
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
GENOME_SIZE = 4.4e6  # M. tuberculosis complex, used to estimate the sequencing depth
LOW_DEPTH = 20  # Below this depth, present spacers may get fewer reads than the default minimum count

REPORT_HEADER = ['Sample', 'SpacerCount', 'Binary', 'Octal', 'Hexadecimal', 'Spoligotype',
                 'FileType', 'Reads', 'Depth', 'MinCount', 'Status', 'Warnings']


@dataclass
class InputFile:
    path: str  # Absolute path, as given (symbolic links are not resolved)
    size: int
    modified: str
    md5: str = ''
    target: str = ''  # Real file, when path is a symbolic link

    @classmethod
    def describe(cls, path, md5=True):
        stat = os.stat(path)
        absolute = os.path.abspath(path)
        real = os.path.realpath(path)
        return cls(absolute, stat.st_size, datetime.fromtimestamp(stat.st_mtime).astimezone().strftime(
            '%Y-%m-%d %H:%M:%S %Z'), file_md5(path) if md5 else '', real if real != absolute else '')


@dataclass
class Result:
    sample: str
    files: list = field(default_factory=list)  # InputFile
    file_type: str = ''
    min_count: int = 0
    counts: list = field(default_factory=list)  # Number of reads (or contigs) matching each spacer, in order
    binary: str = ''
    octal: str = ''
    hexadecimal: str = ''
    spoligotype: str = ''
    reads: int | None = None  # Reads (or contigs) in the input, from Seal
    bases: int | None = None
    warnings: list = field(default_factory=list)
    error: str = ''
    seconds: float = 0.0

    @property
    def status(self):
        if self.error:
            return 'failed'
        return 'warning' if self.warnings else 'ok'

    @property
    def found(self):
        return self.spoligotype not in ('', NOT_FOUND)

    @property
    def depth(self):
        """Estimated sequencing depth for reads, assuming all the reads are from the sample: None for assemblies."""
        if self.file_type != 'fastq' or not self.bases:
            return None
        return self.bases / GENOME_SIZE

    @property
    def median_present_count(self):
        present = [c for c, bit in zip(self.counts, self.binary, strict=False) if bit == '1']
        return statistics.median(present) if present else 0

    def warn(self, message, *args):
        message = message % args
        self.warnings.append(message)
        log.warning('%s: %s', self.sample, message)

    def row(self):
        depth = '' if self.depth is None else '{:.0f}'.format(self.depth)
        notes = self.error.splitlines()[0] if self.error else ' | '.join(self.warnings)
        return [self.sample, ':'.join(str(c) for c in self.counts), self.binary, self.octal, self.hexadecimal,
                self.spoligotype, self.file_type, '' if self.reads is None else str(self.reads), depth,
                str(self.min_count or ''), self.status, ' '.join(notes.split())]


@dataclass
class RunInfo:
    """Everything needed to trace how a report was produced."""
    command: str
    operator: str = ''
    parameters: dict = field(default_factory=dict)
    started: datetime = field(default_factory=lambda: datetime.now().astimezone())
    finished: datetime | None = None
    user: str = ''
    host: str = ''
    system: str = ''
    working_directory: str = ''
    software: dict = field(default_factory=dict)
    database: dict = field(default_factory=dict)
    spacers: dict = field(default_factory=dict)

    @classmethod
    def collect(cls, command, database=SPOLIGOTYPE_DB, spacers=SPACERS_FASTA, operator=None, parameters=None):
        try:
            user = getpass.getuser()
        except (KeyError, OSError):  # No user name, e.g. in some containers
            user = str(os.getuid()) if hasattr(os, 'getuid') else 'unknown'
        info = cls(command, operator or user, dict(parameters or {}), user=user, host=socket.gethostname(),
                   system=platform.platform(), working_directory=os.getcwd())
        info.software = {'spoligotyper': __version__,
                         'Python': '{} ({})'.format(platform.python_version(), sys.executable),
                         'Seal': seal.executable() or 'not found', **seal.versions()}
        info.database = {'path': str(Path(str(database)).resolve()), 'md5': file_md5(database),
                         'patterns': len(load_database(database))}
        info.spacers = {'path': str(Path(str(spacers)).resolve()), 'md5': file_md5(spacers),
                        'spacers': len(read_spacer_names(spacers))}
        return info


def file_md5(path):
    digest = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


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


def check_inputs(r1, r2=None):
    """Validate the input files and return their type ("fasta" or "fastq")."""
    for path in (r1, r2):
        if path is not None and not Path(path).is_file():
            raise SpoligoError('Input file not found: {}'.format(path))
    kind = file_type(r1)
    if r2 is not None:
        if Path(r1).resolve() == Path(r2).resolve():
            raise SpoligoError('R1 and R2 are the same file: {}'.format(r1))
        if kind != 'fastq' or file_type(r2) != 'fastq':
            raise SpoligoError('R2 is only for paired-end fastq files; both R1 and R2 must be fastq')
    return kind


def spoligotype(r1, r2=None, sample=None, min_count=None, threads=1, memory='1g', database=SPOLIGOTYPE_DB,
                spacers=SPACERS_FASTA, md5=False):
    """
    Spoligotype a sample from single-end or paired-end reads, or from an assembly.

    :param min_count: minimum number of reads (or contigs) matching a spacer to call it present.
                      Default: 5 for fastq files, 1 for fasta files.
    :param md5: compute the MD5 checksum of the input files, for the report
    :return: Result. Errors are raised, not stored in the result.
    """
    start = time.monotonic()
    kind = check_inputs(r1, r2)
    inputs = [f for f in (r1, r2) if f]
    result = Result(sample or sample_name(r1, kind), [InputFile.describe(f, md5) for f in inputs], kind)
    result.min_count = MIN_COUNT[kind] if min_count is None else min_count
    spacer_names = read_spacer_names(spacers)
    db = load_database(database)  # Before running Seal, to report a bad database right away

    log.info('Spoligotyping %s (%s, minimum count %d)', result.sample, kind, result.min_count)
    stats = seal.count_spacers(inputs, spacers, threads=threads, memory=memory)
    result.reads, result.bases = stats.reads, stats.bases
    result.counts = [stats.counts.get(name, 0) for name in spacer_names]
    result.binary = to_binary(stats.counts, spacer_names, result.min_count)
    result.octal, result.hexadecimal = binary_to_octal(result.binary), binary_to_hex(result.binary)
    result.spoligotype = lookup(result.binary, db)
    check_result(result)
    result.seconds = time.monotonic() - start
    log.info('%s: %s (octal %s)', result.sample, result.spoligotype, result.octal)
    return result


def check_result(result):
    """Add warnings for results that deserve a second look."""
    if result.file_type == 'fasta' and result.min_count > 1:
        result.warn('fasta file typed with minimum count %d: spacers are probably missed. '
                    'Leave --min-count unset (1 for fasta files).', result.min_count)
    if not any(result.counts):
        result.warn('no spacer found. Is it a Mycobacterium tuberculosis complex sample?')
        return
    if result.file_type != 'fastq':
        return
    if result.depth is not None and result.depth < LOW_DEPTH:
        result.warn('estimated depth %.0fx: present spacers may be missed.', result.depth)
    borderline = [i + 1 for i, c in enumerate(result.counts) if 0 < c < result.min_count]
    if borderline:
        result.warn('%d spacer(s) called absent were seen in fewer than %d reads: %s. Low depth or contamination? '
                    'See the spacer counts.', len(borderline), result.min_count,
                    ', '.join(str(i) for i in borderline))


def spoligotype_samples(samples, **kwargs):
    """
    Spoligotype several samples, one after the other. A sample that fails is reported and the others still run.

    :param samples: list of Sample
    :param kwargs: passed to spoligotype()
    :return: list of Result, in the same order
    """
    results = []
    for i, sample in enumerate(samples, 1):
        log.info('Sample %d of %d: %s', i, len(samples), sample.name)
        try:
            results.append(spoligotype(*sample.files, sample=sample.name, **kwargs))
        except (seal.SealError, SpoligoError, OSError) as e:
            log.error('%s: %s', sample.name, e)
            files = [InputFile.describe(f, md5=False) for f in sample.files if Path(f).is_file()]
            results.append(Result(sample.name, files, sample.file_type, error=str(e)))
    return results


def write_tsv(results, path):
    """Tab-separated report, one line per sample."""
    with open(path, 'w') as f:
        f.write('\t'.join(REPORT_HEADER) + '\n')
        for result in results:
            f.write('\t'.join(result.row()) + '\n')

