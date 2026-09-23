"""Count reads (or contigs) matching each spacer with Seal, from BBTools."""

import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from functools import cache
from pathlib import Path

log = logging.getLogger(__name__)

KMER_SIZE = 25  # Spacers are 25 bp long: each spacer is a single k-mer


class SealError(Exception):
    pass


@cache
def executable():
    """
    Path of seal.sh: the one in the PATH, otherwise the one installed next to the running Python.
    The fallback makes "/path/to/env/bin/spoligotyper" work without activating the conda environment.

    :return: path to seal.sh, or None if it cannot be found
    """
    found = shutil.which('seal.sh')
    if found:
        return found
    # Do not resolve symlinks: a virtualenv's python is a symlink to the base interpreter
    candidate = Path(sys.executable).parent / 'seal.sh'
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def check_seal():
    """Make sure seal.sh is installed and return its path."""
    path = executable()
    if path is None:
        raise SealError('"seal.sh" was not found in your PATH or in "{}". Install BBTools in the environment '
                        'with "conda install -c bioconda bbmap".'.format(Path(sys.executable).parent))
    return path


def seal_command(inputs, spacers_fasta, stats_file, threads, memory):
    cmd = [check_seal(), '-Xmx{}'.format(memory), 'in={}'.format(inputs[0])]
    if len(inputs) > 1:
        cmd.append('in2={}'.format(inputs[1]))
    cmd += ['ref={}'.format(spacers_fasta),
            'k={}'.format(KMER_SIZE),
            'rcomp=t',
            'hdist=1',  # Up to 1 mismatch
            'maskmiddle=f',  # Do not treat the middle base of a k-mer as a wildcard
            'clearzone=999999',
            'ambiguous=all',  # Count a read for every spacer it matches
            'nzo=f',  # Also report spacers with no match
            'qin=33',  # Force the quality encoding: autodetection fails on some low quality nanopore reads
            'ow=t',
            'stats={}'.format(stats_file),
            'threads={}'.format(threads)]
    return cmd


@dataclass
class SealStats:
    counts: dict  # {spacer name: number of matching reads}
    reads: int | None = None  # Total number of reads (or contigs) in the input
    bases: int | None = None


def parse_stats(stats_file):
    """
    Read Seal's stats file.

    #File	sample_R1.fastq.gz
    #Total	822714	123407100
    #Matched	799	0.09712%
    #Name	Reads	ReadsPct
    spacer25	62	0.00754%
    """
    stats = SealStats({})
    with open(stats_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if not line.strip():
                continue
            try:
                if fields[0] == '#Total':
                    stats.reads = int(fields[1])
                    stats.bases = int(fields[2]) if len(fields) > 2 else None
                elif not line.startswith('#'):
                    stats.counts[fields[0]] = int(fields[1])
            except (IndexError, ValueError):
                raise SealError('Unexpected line in Seal stats file {}: {}'.format(stats_file, line.strip())) from None
    return stats


def count_spacers(inputs, spacers_fasta, threads=1, memory='1g'):
    """
    Count the reads (or contigs) that contain each spacer, allowing one mismatch.

    :param inputs: one fasta/fastq file, or two paired-end fastq files
    :param spacers_fasta: fasta file of the spacer sequences
    :param memory: Java heap size given to Seal, e.g. "1g"
    :return: SealStats
    """
    with tempfile.TemporaryDirectory(prefix='spoligotyper_') as tmp:
        stats_file = Path(tmp) / 'stats.tsv'
        cmd = seal_command(inputs, spacers_fasta, stats_file, threads, memory)
        log.debug('Running: %s', shlex.join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not stats_file.exists():
            detail = (proc.stderr or proc.stdout).strip().splitlines() or ['no output']
            reason = next((line for line in reversed(detail) if re.search(r'Exception|Error|error', line)), detail[-1])
            raise SealError('Seal failed (exit code {}): {}\nCommand: {}\n{}'.format(
                proc.returncode, reason.strip(), shlex.join(cmd), '\n'.join(detail[-10:])))
        log.debug('Seal output:\n%s', proc.stderr.strip())
        return parse_stats(stats_file)


def versions():
    """{"BBTools": version, "Java": version} of the Seal installation, "unknown" when they cannot be read."""
    found = {'BBTools': 'unknown', 'Java': 'unknown'}
    path = executable()
    if path is None:
        return found
    try:  # A small heap: without -Xmx, seal.sh reserves most of the free memory, even for --version
        out = subprocess.run([path, '-Xmx64m', '--version'], capture_output=True, text=True, timeout=60)
        match = re.search(r'(?:BBMap|BBTools) version (\S+)', out.stdout + out.stderr)
        if match:
            found['BBTools'] = match.group(1)
        java = Path(path).parent / 'java'
        java = str(java) if java.exists() else shutil.which('java')
        if java:
            out = subprocess.run([java, '-version'], capture_output=True, text=True, timeout=60)
            first = (out.stderr or out.stdout).strip().splitlines()
            if first:
                found['Java'] = first[0]
    except (OSError, subprocess.SubprocessError):
        pass
    return found
