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


def seal_command(inputs, refs, stats_file, threads, memory, k=KMER_SIZE, hdist=1):
    """
    :param refs: reference fasta file(s): a path or a list of paths
    :param k: k-mer size. Spacers are 25 bp long: each spacer is a single k-mer
    :param hdist: number of mismatches allowed in a k-mer
    """
    refs = [refs] if isinstance(refs, (str, os.PathLike)) else list(refs)
    cmd = [check_seal(), '-Xmx{}'.format(memory), 'in={}'.format(inputs[0])]
    if len(inputs) > 1:
        cmd.append('in2={}'.format(inputs[1]))
    else:
        cmd.append('int=f')  # Never guess that a single fastq file is interleaved: that would double the counts
    cmd += ['ref={}'.format(','.join(str(r) for r in refs)),
            'k={}'.format(k),
            'rcomp=t',
            'hdist={}'.format(hdist),
            'maskmiddle=f',  # Do not treat the middle base of a k-mer as a wildcard
            'clearzone=999999',
            'ambiguous=all',  # Count a read for every reference sequence it matches
            'nzo=f',  # Also report reference sequences with no match
            'qin=33',  # Force the quality encoding: autodetection fails on some low quality nanopore reads
            'ow=t',
            'stats={}'.format(stats_file),
            'threads={}'.format(threads)]
    return cmd


UNSAFE = re.compile(r'[\s,=]')  # seal.sh splits its arguments on spaces, and ref= lists on commas


def safe_paths(paths, folder, prefix):
    """
    Paths that Seal can read: those with a space, comma or "=" are linked into folder under a safe name, keeping the
    extensions that Seal uses to recognize the format (e.g. ".fastq.gz").
    """
    safe = []
    for i, path in enumerate(paths):
        path = str(path)
        if UNSAFE.search(path):
            link = Path(folder) / '{}{}_{}'.format(prefix, i, re.sub(r'[^A-Za-z0-9._-]', '_', Path(path).name))
            link.symlink_to(Path(path).resolve())
            path = str(link)
        safe.append(path)
    return safe


GENERIC_ERRORS = re.compile(r'Exception in thread "main"|terminated in an error state')


def failure_reason(lines):
    """
    The most informative line of Seal's error output: the first error that is not Seal's generic final message.
    Java exceptions often end with ":" and give their message on the next line, which is then added.
    """
    candidates = [i for i, line in enumerate(lines) if re.search(r'Exception|Error|error', line)
                  and not line.strip().startswith('at ')]
    specific = [i for i in candidates if not GENERIC_ERRORS.search(lines[i])]
    if not candidates:
        return lines[-1].strip()
    i = (specific or candidates)[0]
    reason = lines[i].strip()
    if reason.endswith(':'):
        following = [line.strip() for line in lines[i + 1:] if line.strip() and not line.strip().startswith('at ')]
        if following:
            reason += ' ' + following[0]
    return reason


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
                    stats.counts[fields[0].split()[0]] = int(fields[1])  # Name without the description
            except (IndexError, ValueError):
                raise SealError('Unexpected line in Seal stats file {}: {}'.format(stats_file, line.strip())) from None
    return stats


def run_seal(inputs, refs, threads=1, memory='1g', k=KMER_SIZE, hdist=1):
    """
    Count the reads (or contigs) that contain each reference sequence.

    :param inputs: one fasta/fastq file, or two paired-end fastq files
    :param refs: reference fasta file(s)
    :param memory: Java heap size given to Seal, e.g. "1g"
    :return: SealStats
    """
    refs = [refs] if isinstance(refs, (str, os.PathLike)) else list(refs)
    with tempfile.TemporaryDirectory(prefix='spoligotyper_') as tmp:
        if UNSAFE.search(tmp):
            raise SealError('The temporary folder "{}" contains a space or a comma, which Seal cannot handle. '
                            'Set TMPDIR to another folder.'.format(tmp))
        stats_file = Path(tmp) / 'stats.tsv'
        cmd = seal_command(safe_paths(inputs, tmp, 'in'), safe_paths(refs, tmp, 'ref'), stats_file, threads, memory,
                           k=k, hdist=hdist)
        log.debug('Running: %s', shlex.join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0 or not stats_file.exists():
            detail = (proc.stderr or proc.stdout).strip().splitlines() or ['no output']
            raise SealError('Seal failed (exit code {}): {}\nCommand: {}\n{}'.format(
                proc.returncode, failure_reason(detail), shlex.join(cmd), '\n'.join(detail[-10:])))
        log.debug('Seal output:\n%s', proc.stderr.strip())
        return parse_stats(stats_file)


def count_spacers(inputs, spacers_fasta, threads=1, memory='1g'):
    """Count the reads (or contigs) that contain each spacer, allowing one mismatch."""
    return run_seal(inputs, spacers_fasta, threads=threads, memory=memory)


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
