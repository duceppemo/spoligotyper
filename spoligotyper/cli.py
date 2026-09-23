"""Command line entry point."""

import logging
import os
import shlex
import sys
import time
from argparse import ArgumentParser, ArgumentTypeError, RawDescriptionHelpFormatter
from pathlib import Path

from . import __version__
from .pipeline import REPORT_HEADER, spoligotype, write_report
from .seal import SealError
from .spoligotype import SPOLIGOTYPE_DB, SpoligoError

log = logging.getLogger('spoligotyper')

# Expected errors are reported without a traceback
USER_ERRORS = (SealError, SpoligoError, OSError)

EPILOG = '''examples:
  spoligotyper -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/
  spoligotyper -r1 sample.fastq.gz -o results/
  spoligotyper -r1 assembly.fasta -o results/

documentation: https://github.com/duceppemo/spoligotyper/wiki'''


def setup_logging(verbose=False):
    """Configure the package logger only, so importing spoligotyper never alters the root logger."""
    for handler in list(log.handlers):
        log.removeHandler(handler)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)-7s %(message)s', datefmt='%H:%M:%S'))
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)


def positive_int(value):
    try:
        value = int(value)
    except ValueError:
        raise ArgumentTypeError('"{}" is not an integer'.format(value)) from None
    if value < 1:
        raise ArgumentTypeError('must be >= 1')
    return value


def java_memory(value):
    if not value[:-1].isdigit() or value[-1].lower() not in 'mg':
        raise ArgumentTypeError('"{}" is not a Java memory size such as 500m or 2g'.format(value))
    return value


def available_cpus():
    try:
        return len(os.sched_getaffinity(0))  # Respects cgroup/taskset limits on Linux (e.g. SLURM jobs)
    except AttributeError:
        return os.cpu_count() or 1


def build_parser():
    max_cpu = available_cpus()
    parser = ArgumentParser(prog='spoligotyper', formatter_class=RawDescriptionHelpFormatter, epilog=EPILOG,
                            description='In silico spoligotyping of Mycobacterium tuberculosis complex samples\n'
                                        'from sequencing reads (fastq) or assemblies (fasta).')
    parser.add_argument('-r1', '--r1', metavar='FILE', required=True,
                        help='Reads (single-end or R1 of paired-end) or assembly. fastq or fasta, gzipped or not.')
    parser.add_argument('-r2', '--r2', metavar='FILE',
                        help='R2 reads, for paired-end fastq files.')
    parser.add_argument('-o', '--output', metavar='FOLDER', required=True,
                        help='Folder to hold the report. Created if needed.')
    parser.add_argument('-s', '--sample', metavar='NAME',
                        help='Sample name used in the report and its file name. '
                             'Default: the input file name, without extension and read suffix (_R1, _1, ...).')
    parser.add_argument('-m', '--min-count', metavar='N', type=positive_int,
                        help='Minimum number of reads matching a spacer to call it present. '
                             'Default: 5 for fastq files, 1 for fasta files.')
    parser.add_argument('-t', '--threads', metavar='N', type=positive_int, default=max_cpu,
                        help='Number of threads. Default: all available ({}).'.format(max_cpu))
    parser.add_argument('--memory', metavar='SIZE', type=java_memory, default='1g',
                        help='Memory for Seal (Java heap size). Default: 1g.')
    parser.add_argument('--db', metavar='FILE', default=SPOLIGOTYPE_DB,
                        help='Spoligotype database: "octal SB-number binary" on each line. '
                             'Default: the Mbovis.org database included with spoligotyper.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show debug messages, including the Seal command.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.sample is not None and (not args.sample.strip() or '/' in args.sample or args.sample in ('.', '..')):
        parser.error('--sample must be a name, not a path: "{}"'.format(args.sample))
    setup_logging(args.verbose)
    log.debug('spoligotyper %s: %s', __version__, shlex.join([parser.prog, *(sys.argv[1:] if argv is None else argv)]))
    start = time.monotonic()

    try:
        result = spoligotype(args.r1, args.r2, sample=args.sample, min_count=args.min_count,
                             threads=args.threads, memory=args.memory, database=args.db)
        output = Path(args.output).expanduser()
        output.mkdir(parents=True, exist_ok=True)
        report = output / '{}_spoligotyping.txt'.format(result.sample)
        write_report([result], report)
    except USER_ERRORS as e:
        log.error('%s', e)
        log.debug('Traceback:', exc_info=True)  # Shown with --verbose
        sys.exit(1)
    except KeyboardInterrupt:
        log.error('Interrupted')
        sys.exit(130)

    # The result goes to stdout, messages to stderr: "spoligotyper ... > result.tsv" works
    print('\t'.join(REPORT_HEADER))
    print('\t'.join(result.row()))
    sys.stdout.flush()
    log.info('%s: %s (octal %s)', result.sample, result.spoligotype, result.octal)
    log.info('Report saved in %s (%.1f s)', report, time.monotonic() - start)
