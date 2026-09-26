"""Command line entry point."""

import logging
import os
import shlex
import sys
import time
from argparse import ArgumentParser, ArgumentTypeError, RawDescriptionHelpFormatter
from datetime import datetime
from pathlib import Path

from . import __version__, sitdb
from .pipeline import REPORT_HEADER, RunInfo, spoligotype, spoligotype_samples, write_json, write_multiqc, write_tsv
from .samples import find_samples
from .seal import SealError, check_seal
from .spoligotype import SPOLIGOTYPE_DB, SpoligoError

log = logging.getLogger('spoligotyper')

# Expected errors are reported without a traceback
USER_ERRORS = (SealError, SpoligoError, OSError)

EPILOG = '''examples:
  spoligotyper -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz -o results/
  spoligotyper -r1 sample.fastq.gz -o results/
  spoligotyper -r1 assembly.fasta -o results/
  spoligotyper -i folder_of_fastq_and_fasta/ -o results/

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
    inputs = parser.add_argument_group('input (one sample with -r1, or a folder of samples with -i)')
    source = inputs.add_mutually_exclusive_group(required=True)
    source.add_argument('-r1', '--r1', metavar='FILE',
                        help='Reads (single-end or R1 of paired-end) or assembly. fastq or fasta, gzipped or not.')
    source.add_argument('-i', '--input', metavar='FOLDER',
                        help='Folder of fastq and fasta files, searched recursively. Each fasta file, single-end '
                             'fastq file or pair of R1/R2 fastq files is a sample.')
    inputs.add_argument('-r2', '--r2', metavar='FILE',
                        help='R2 reads, for paired-end fastq files.')
    inputs.add_argument('-s', '--sample', metavar='NAME',
                        help='With -r1: sample name used in the reports. Default: the input file name, without '
                             'extension and read suffix (_R1, _1, ...).')

    output = parser.add_argument_group('output')
    output.add_argument('-o', '--output', metavar='FOLDER', required=True,
                        help='Folder to hold the reports. Created if needed.')
    output.add_argument('--no-pdf', action='store_true', help='Do not write the PDF report.')
    output.add_argument('--no-md5', action='store_true',
                        help='Do not compute the MD5 checksums of the input files (PDF and JSON reports).')
    output.add_argument('--operator', metavar='NAME',
                        help='Name of the person running the analysis, shown in the PDF report. Default: user name.')

    typing = parser.add_argument_group('spoligotyping')
    typing.add_argument('-m', '--min-count', metavar='N', type=positive_int,
                        help='Minimum number of reads matching a spacer to call it present. '
                             'Default: 5 for fastq files, 1 for fasta files.')
    typing.add_argument('--no-species', action='store_true',
                        help='Skip the species check (regions of difference RD1, RD4, RD7, RD9, RD12) and the lineage '
                             '(SNP barcode). Faster: one pass over the reads instead of two.')
    typing.add_argument('--sit-db', metavar='FILE',
                        help='SIT database for the SIT and SITVIT2family columns. Default: the one saved by '
                             'spoligotyper-download-sit, if any.')
    typing.add_argument('--db', metavar='FILE', default=SPOLIGOTYPE_DB,
                        help='Spoligotype database: "octal SB-number binary" on each line. '
                             'Default: the Mbovis.org database included with spoligotyper.')

    other = parser.add_argument_group('performance and other options')
    other.add_argument('-t', '--threads', metavar='N', type=positive_int, default=max_cpu,
                       help='Number of threads. Default: all available ({}).'.format(max_cpu))
    other.add_argument('-j', '--jobs', metavar='N', type=positive_int, default=1,
                       help='With -i: number of samples typed at the same time, sharing the threads. Each job uses '
                            'the --memory given to Seal. Default: 1.')
    other.add_argument('--memory', metavar='SIZE', type=java_memory, default='1g',
                       help='Memory for Seal (Java heap size), per job. Default: 1g.')
    other.add_argument('-v', '--verbose', action='store_true', help='Show debug messages, including the Seal command.')
    other.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    return parser


def check_arguments(parser, args):
    if args.input and (args.r2 or args.sample):
        parser.error('-r2 and --sample are only used with -r1')
    if args.r1 and args.jobs > 1:
        parser.error('--jobs is only used with -i')
    if args.sample is not None and (not args.sample.strip() or '/' in args.sample or args.sample in ('.', '..')):
        parser.error('--sample must be a name, not a path: "{}"'.format(args.sample))


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    check_arguments(parser, args)
    setup_logging(args.verbose)
    command = shlex.join([parser.prog, *(sys.argv[1:] if argv is None else argv)])
    log.debug('spoligotyper %s: %s', __version__, command)
    start = time.monotonic()
    output = Path(args.output).expanduser()
    try:
        check_seal()
        sit_path = args.sit_db or sitdb.default_database()
        sit_db = sitdb.load(sit_path) if sit_path else None
        if sit_db is None:
            log.info('No SIT database: run spoligotyper-download-sit once for the SIT and SITVIT2 family columns.')
        options = dict(min_count=args.min_count, threads=args.threads, memory=args.memory, database=args.db,
                       md5=not args.no_md5, species_check=not args.no_species, sit_db=sit_db)
        run = RunInfo.collect(command, database=args.db, operator=args.operator, sit_db=sit_db, parameters={
            'Input': ' '.join(os.path.abspath(f) for f in (args.input, args.r1, args.r2) if f),
            'Output folder': str(output.resolve()),
            'Minimum count': args.min_count or 'default (5 for fastq, 1 for fasta)',
            'Threads': args.threads, 'Jobs': args.jobs, 'Seal memory': args.memory,
            'Species and lineage': 'no' if args.no_species else 'yes'})
        if args.input:
            samples = find_samples(args.input, exclude=[output])
            log.info('%d sample(s) found in %s', len(samples), args.input)
            results = spoligotype_samples(samples, jobs=args.jobs, **options)
            prefix = 'spoligotyping'
            tsv, pdf = output / (prefix + '.tsv'), output / (prefix + '_report.pdf')
        else:
            results = [spoligotype(args.r1, args.r2, sample=args.sample, **options)]
            prefix = '{}_spoligotyping'.format(results[0].sample)
            tsv, pdf = output / (prefix + '.txt'), output / (prefix + '.pdf')
        run.finished = datetime.now().astimezone()

        output.mkdir(parents=True, exist_ok=True)
        write_tsv(results, tsv)
        write_json(results, run, output / (prefix + '.json'))
        write_multiqc(results, output / (prefix + '_mqc.json'))
        reports = [tsv, output / (prefix + '.json')]
        if not args.no_pdf:
            from .pdf import write_pdf  # reportlab is only imported when needed
            write_pdf(results, run, pdf)
            reports.append(pdf)
    except USER_ERRORS as e:
        log.error('%s', e)
        log.debug('Traceback:', exc_info=True)  # Shown with --verbose
        sys.exit(1)
    except KeyboardInterrupt:
        log.error('Interrupted')
        sys.exit(130)

    # The results go to stdout, messages to stderr: "spoligotyper ... > results.tsv" works
    print('\t'.join(REPORT_HEADER))
    for result in results:
        print('\t'.join(result.row()))
    sys.stdout.flush()
    for report in reports:
        log.info('Report saved in %s', report)
    log.info('Done in %.1f s', time.monotonic() - start)
    failed = [r.sample for r in results if r.error]
    if failed:
        log.error('%d of %d sample(s) failed: %s', len(failed), len(results), ', '.join(failed))
        sys.exit(1)
