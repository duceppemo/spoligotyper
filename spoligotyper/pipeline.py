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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

from . import __version__, lineage, seal, species
from .lineage import LineageCall
from .samples import sample_name
from .species import SpeciesCheck
from .spoligotype import (
    NOT_FOUND,
    SPACERS_FASTA,
    SPOLIGOTYPE_DB,
    SpoligoError,
    binary_to_hex,
    binary_to_octal,
    closest,
    describe_closest,
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
LOW_MTBC_FRACTION = 0.6  # Below this estimated fraction of MTBC reads, the sample is probably contaminated
WEAK_SPACER = 0.4  # Present spacers with fewer reads than this fraction of the median: mixed sample?
WEAK_SPACER_MIN_MEDIAN = 30  # ... when the median is high enough for this not to happen by chance

REPORT_HEADER = ['Sample', 'SpacerCount', 'Binary', 'Octal', 'Hexadecimal', 'Spoligotype',
                 'FileType', 'Reads', 'Depth', 'MinCount', 'Status', 'Warnings',
                 'Species', 'Lineage', 'LineageName', 'RD9', 'RD4', 'RD1', 'MTBCFraction', 'Closest']


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
    closest: list = field(default_factory=list)  # Closest database patterns when not found: [(SB, [spacers])]
    species: SpeciesCheck | None = None
    lineage: LineageCall | None = None

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
    def paired(self):
        return len(self.files) == 2

    @property
    def read_length(self):
        return self.bases / self.reads if self.file_type == 'fastq' and self.reads and self.bases else None

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
        check, call = self.species, self.lineage
        regions = [check.state(r) if check else '' for r in species.REGIONS]
        fraction = '' if check is None or check.mtbc_fraction is None else '{:.2f}'.format(check.mtbc_fraction)
        return [self.sample, ':'.join(str(c) for c in self.counts), self.binary, self.octal, self.hexadecimal,
                self.spoligotype, self.file_type, '' if self.reads is None else str(self.reads), depth,
                str(self.min_count or ''), self.status, ' '.join(notes.split()),
                check.species if check else '', call.lineage if call else '', call.name if call else '',
                *regions, fraction, describe_closest(self.closest)]

    def to_dict(self):
        """Everything about the result, for the JSON report."""
        data = asdict(self)
        data.update(status=self.status, depth=self.depth, spacers_present=self.binary.count('1'),
                    median_present_count=self.median_present_count)
        data['closest'] = [{'spoligotype': name, 'differing_spacers': diff} for name, diff in self.closest]
        if self.species:
            data['species']['regions'] = {r: {'state': state, 'depth_ratio': round(ratio, 3)}
                                          for r, (state, ratio) in self.species.regions.items()}
        return data


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
    species_data: dict = field(default_factory=dict)

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
        info.species_data = {name: {'path': str(Path(str(path)).resolve()), 'md5': file_md5(path)}
                             for name, path in (('Species markers', species.MARKERS_FASTA),
                                                ('Lineage SNP barcode', lineage.BARCODE),
                                                ('Lineage SNP sequences', lineage.LINEAGE_SNPS_FASTA))}
        return info

    def to_dict(self):
        data = asdict(self)
        for key in ('started', 'finished'):
            data[key] = data[key].isoformat(timespec='seconds') if data[key] else None
        return data


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
                spacers=SPACERS_FASTA, md5=False, species_check=True):
    """
    Spoligotype a sample from single-end or paired-end reads, or from an assembly.

    :param min_count: minimum number of reads (or contigs) matching a spacer to call it present.
                      Default: 5 for fastq files, 1 for fasta files.
    :param md5: compute the MD5 checksum of the input files, for the report
    :param species_check: also check the species (regions of difference) and the lineage (SNP barcode)
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
    # Spacers and species markers in one pass: both are searched with 25-mers and 1 mismatch
    refs = [spacers, species.MARKERS_FASTA] if species_check else [spacers]
    stats = seal.run_seal(inputs, refs, threads=threads, memory=memory)
    result.reads, result.bases = stats.reads, stats.bases
    result.counts = [stats.counts.get(name, 0) for name in spacer_names]
    result.binary = to_binary(stats.counts, spacer_names, result.min_count)
    result.octal, result.hexadecimal = binary_to_octal(result.binary), binary_to_hex(result.binary)
    result.spoligotype = lookup(result.binary, db)
    result.closest = closest(result.binary, db) if any(result.counts) else []
    if species_check:
        result.species = species.check_species(stats.counts, kind, result.depth, result.read_length, result.paired)
        result.lineage = lineage.LineageCall()
        if result.species.mtbc:  # Some lineage SNPs are conserved in other mycobacteria: no lineage without MTBC
            # Lineage SNPs need exact matches: a second pass with 31-mers and no mismatch
            snps = seal.run_seal(inputs, lineage.LINEAGE_SNPS_FASTA, threads=threads, memory=memory,
                                 k=lineage.KMER_SIZE, hdist=0)
            fraction = result.species.mtbc_fraction
            result.lineage = lineage.call_lineage(snps.counts, kind, contaminated=fraction is not None and
                                                  fraction < lineage.CONTAMINATED_FRACTION)
            species.name_species(result.species, result.lineage.called,
                                 mixed=bool(result.lineage.mixed or result.lineage.conflict))
    check_result(result)
    result.seconds = time.monotonic() - start
    log.info('%s: %s (octal %s)%s', result.sample, result.spoligotype, result.octal,
             ', {}, lineage {}'.format(result.species.species, result.lineage.lineage or '-')
             if result.species else '')
    return result


def check_result(result):
    """Add warnings for results that deserve a second look."""
    if result.file_type == 'fasta' and result.min_count > 1:
        result.warn('fasta file typed with minimum count %d: spacers are probably missed. '
                    'Leave --min-count unset (1 for fasta files).', result.min_count)
    check_species(result)
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
    median = result.median_present_count
    weak = [i + 1 for i, (c, bit) in enumerate(zip(result.counts, result.binary, strict=True))
            if bit == '1' and c < WEAK_SPACER * median]
    if weak and median >= WEAK_SPACER_MIN_MEDIAN:
        result.warn('spacer(s) %s have less than %d%% of the median read count of the present spacers (%g): '
                    'mixed sample?', ', '.join(str(i) for i in weak), WEAK_SPACER * 100, median)


def check_species(result):
    check, call = result.species, result.lineage
    if check is None:
        return
    if not check.mtbc:
        if any(result.counts):
            result.warn('too little MTBC DNA to check the species (median %g reads on the MTBC control regions).',
                        check.control_depth)
        return
    if check.mtbc_fraction is not None and check.mtbc_fraction < LOW_MTBC_FRACTION:
        result.warn('only about %d%% of the reads appear to be from the M. tuberculosis complex: contamination?',
                    round(check.mtbc_fraction * 100))
    partial = [r for r in species.REGIONS if check.state(r) == species.PARTIAL]
    if partial:
        result.warn('%s partially deleted (depth ratio %s): mixed sample?', ', '.join(partial),
                    ', '.join('{:.2f}'.format(check.regions[r][1]) for r in partial))
    for warning in species.consistency_warnings(check, call.called):
        result.warn(warning)
    if call.conflict:
        result.warn('SNPs of several lineages (%s): mixed sample?', ', '.join(sorted(call.called)))
    if call.mixed:
        result.warn('both alleles of %d lineage SNP(s) seen (%s): mixed sample?', len(call.mixed),
                    ', '.join('{} {:.0f}%'.format(s.lineage, s.fraction * 100) for s in call.mixed))
    if result.found and check.state('RD9') == species.PRESENT:
        result.warn('%s is an SB number, but RD9 is present: SB numbers are for RD9-deleted (animal) lineages.',
                    result.spoligotype)


def spoligotype_samples(samples, jobs=1, **kwargs):
    """
    Spoligotype several samples. A sample that fails is reported and the others still run.

    :param samples: list of Sample
    :param jobs: number of samples typed at the same time. The threads are shared between them.
    :param kwargs: passed to spoligotype()
    :return: list of Result, in the same order
    """
    jobs = max(1, min(jobs, len(samples)))
    if jobs > 1:
        kwargs['threads'] = max(1, kwargs.get('threads', 1) // jobs)

    def run(numbered):
        i, sample = numbered
        log.info('Sample %d of %d: %s', i, len(samples), sample.name)
        try:
            return spoligotype(*sample.files, sample=sample.name, **kwargs)
        except (seal.SealError, SpoligoError, OSError) as e:
            log.error('%s: %s', sample.name, e)
            files = [InputFile.describe(f, md5=False) for f in sample.files if Path(f).is_file()]
            return Result(sample.name, files, sample.file_type, error=str(e))

    with ThreadPoolExecutor(max_workers=jobs) as pool:  # Threads are enough: the work is done by Seal
        return list(pool.map(run, enumerate(samples, 1)))


def write_tsv(results, path):
    """Tab-separated report, one line per sample."""
    with open(path, 'w') as f:
        f.write('\t'.join(REPORT_HEADER) + '\n')
        for result in results:
            f.write('\t'.join(result.row()) + '\n')



def write_json(results, run, path):
    """Everything in one JSON file, for pipelines."""
    import json
    data = {'spoligotyper': __version__, 'run': run.to_dict(), 'samples': [r.to_dict() for r in results]}
    with open(path, 'w') as f:
        json.dump(data, f, indent=2, default=str)
        f.write('\n')


def write_multiqc(results, path):
    """
    MultiQC custom content: a "Spoligotyping" table in MultiQC reports. JSON rather than TSV, so that octal codes
    such as 000000000003771 stay text instead of being read as numbers.
    """
    import json
    columns = ('Spoligotype', 'Octal', 'Species', 'Lineage', 'Status')
    data = {}
    for r in results:
        values = (r.spoligotype, r.octal, r.species.species if r.species else '',
                  r.lineage.lineage if r.lineage else '', r.status)
        data[r.sample] = {column: value or '-' for column, value in zip(columns, values, strict=True)}
    link = '<a href="https://github.com/duceppemo/spoligotyper">spoligotyper</a> {}'.format(__version__)
    content = {'id': 'spoligotyper', 'section_name': 'Spoligotyping',
               'description': 'In silico spoligotype, species and lineage from {}.'.format(link),
               'plot_type': 'table', 'pconfig': {'id': 'spoligotyper_table', 'namespace': 'spoligotyper'},
               'headers': {column: {'title': column} for column in columns}, 'data': data}
    with open(path, 'w') as f:
        json.dump(content, f, indent=2)
        f.write('\n')
