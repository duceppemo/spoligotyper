"""
Shared international types (SIT) and SITVIT2 spoligotype families.

The SITVIT2 database is not openly licensed, so it is not included with spoligotyper. The SpolLineages tool (Couvin
et al. 2020, GPL-3.0) publishes a list of 9,658 SITVIT2 spoligotype patterns, 3,850 of them with a SIT, all with their
SITVIT2 family. spoligotyper-download-sit downloads this list from the SpolLineages repository (or from its Zenodo
mirror), checks its SHA-256 checksum, and converts it into the SIT database used by spoligotyper.
"""

import csv
import hashlib
import json
import logging
import os
import sys
import tempfile
import urllib.request
from argparse import ArgumentParser
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .spoligotype import BINARY_PATTERN, NOT_FOUND, SpoligoError, binary_to_octal

log = logging.getLogger(__name__)

SOURCE = {
    'name': 'SpolLineages Spoligo_list.csv (SITVIT2 patterns)',
    'repository': 'https://github.com/dcouvin/SpolLineages',
    'commit': '4ef3ef464c0cf453584711e46be14f89dcdcd858',
    'sha256': '68b39726a82442391046d29e3b6d284a9aa048b1e5342206ab05944161723959',
    'license': 'GPL-3.0',
    'citation': 'Couvin D, Segretier W, Stattner E, Rastogi N. Novel methods included in SpolLineages tool for fast '
                'and precise prediction of Mycobacterium tuberculosis complex spoligotype families. Database (Oxford) '
                '2020:baaa108. https://doi.org/10.1093/database/baaa108',
}
URLS = [
    'https://raw.githubusercontent.com/dcouvin/SpolLineages/{}/Spoligo_list.csv'.format(SOURCE['commit']),
    # Zenodo mirror of the same file (same checksum), used when GitHub cannot be reached. Set once published.
]
DB_NAME = 'sit_database.tsv'
ORPHAN = 'Orphan'  # A SITVIT2 pattern without SIT (seen in one isolate only)


class SitError(Exception):
    pass


def default_folder():
    """Where spoligotyper-download-sit saves the database: $SPOLIGOTYPER_DATA, or the user's cache folder."""
    if os.environ.get('SPOLIGOTYPER_DATA'):
        return Path(os.environ['SPOLIGOTYPER_DATA'])
    return Path(os.environ.get('XDG_CACHE_HOME') or Path.home() / '.cache') / 'spoligotyper'


def default_database():
    """The SIT database to use when none is given, or None if it has not been downloaded."""
    path = default_folder() / DB_NAME
    return path if path.is_file() else None


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def download(destination, urls=None):
    """Download the source list from the first URL that works, and check its checksum."""
    errors = []
    for url in urls or URLS:
        try:
            log.info('Downloading %s', url)
            with urllib.request.urlopen(url, timeout=60) as response:
                destination.write_bytes(response.read())
        except OSError as e:
            errors.append('{}: {}'.format(url, e))
            continue
        checksum = sha256(destination)
        if checksum == SOURCE['sha256']:
            return url
        errors.append('{}: checksum {} instead of {}'.format(url, checksum, SOURCE['sha256']))
    raise SitError('Could not download the SITVIT2 pattern list:\n  ' + '\n  '.join(errors))


def convert(source, database, origin=''):
    """Convert the SpolLineages list (";"-separated, "n"/"o" binary patterns) into the spoligotyper SIT database."""
    rows, skipped = [], 0
    with open(source, newline='', encoding='utf-8-sig') as f:
        for i, row in enumerate(csv.DictReader(f, delimiter=';'), 2):
            binary = row['Spoligo Binary'].strip().replace('n', '1').replace('o', '0')
            sit = row['SIT'].strip()
            if '?' in binary and sit == ORPHAN:  # A spacer of unknown state: cannot match a pattern exactly
                skipped += 1
                continue
            if not BINARY_PATTERN.fullmatch(binary) or not (sit.isdigit() or sit == ORPHAN):
                raise SitError('{}, line {}: unexpected pattern or SIT: {}'.format(source, i, row))
            octal = binary_to_octal(binary)
            if octal != row['Spoligo Octal'].strip().lstrip("'"):
                raise SitError('{}, line {}: octal code does not match the binary pattern'.format(source, i))
            rows.append((octal, binary, sit, row['Lineage (SITVIT2)'].strip()))
    with open(database, 'w') as f:
        f.write('# SIT database for spoligotyper, converted from {}\n'.format(SOURCE['name']))
        f.write('# Source: {} (commit {}, SHA-256 {}){}\n'.format(
            SOURCE['repository'], SOURCE['commit'], SOURCE['sha256'], ', from ' + origin if origin else ''))
        f.write('# Converted by spoligotyper {} on {}. License: {}.\n'.format(
            __version__, datetime.now(timezone.utc).strftime('%Y-%m-%d'), SOURCE['license']))
        f.write('# Please cite: {}\n'.format(SOURCE['citation']))
        f.write('octal\tbinary\tsit\tfamily\n')
        for octal, binary, sit, family in rows:
            f.write('\t'.join((octal, binary, sit, family)) + '\n')
    if skipped:
        log.info('%d orphan patterns with spacers of unknown state skipped', skipped)
    return len(rows)


@dataclass
class SitDatabase:
    path: str
    patterns: dict  # {binary: (SIT or "Orphan", SITVIT2 family)}
    sha256: str = ''
    source: str = ''  # First comment line: what it was converted from

    @property
    def sits(self):
        """{binary: SIT} for the patterns with a SIT, e.g. for spoligotype.closest()."""
        return {binary: 'SIT' + sit for binary, (sit, _) in self.patterns.items() if sit != ORPHAN}

    def lookup(self, binary):
        """(SIT, family): e.g. ("SIT451", "T-H37Rv"), ("Orphan", "LAM5"), or (NOT_FOUND, "")."""
        sit, family = self.patterns.get(binary, (None, ''))
        if sit is None:
            return NOT_FOUND, ''
        return ('SIT' + sit if sit != ORPHAN else ORPHAN), family


def load(path):
    patterns, source = {}, ''
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                source = source or line[1:].strip()
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[0] == 'octal':
                continue
            if len(fields) != 4 or not BINARY_PATTERN.fullmatch(fields[1]):
                raise SpoligoError('{}: not a SIT database (run spoligotyper-download-sit)'.format(path))
            patterns[fields[1]] = (fields[2], fields[3])
    if not patterns:
        raise SpoligoError('{} does not contain any pattern'.format(path))
    return SitDatabase(str(Path(path).resolve()), patterns, sha256(path), source)


def main(argv=None):
    """spoligotyper-download-sit: download and convert the SITVIT2 pattern list of SpolLineages."""
    parser = ArgumentParser(prog='spoligotyper-download-sit',
                            description='Download the SITVIT2 spoligotype patterns published with SpolLineages '
                                        '(GPL-3.0), with their shared international type (SIT) and family, for the '
                                        'SIT and SITVIT2family columns of spoligotyper.')
    parser.add_argument('-o', '--output', metavar='FOLDER', default=str(default_folder()),
                        help='Folder for the database. Default: $SPOLIGOTYPER_DATA, or ~/.cache/spoligotyper. '
                             'spoligotyper finds it there; for another folder, use spoligotyper --sit-db.')
    parser.add_argument('--source', metavar='FILE',
                        help='Use a Spoligo_list.csv already downloaded (e.g. on a computer without internet '
                             'access); its checksum is still checked.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)
    logging.basicConfig(format='%(asctime)s %(levelname)-7s %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    folder = Path(args.output).expanduser()
    try:
        folder.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'Spoligo_list.csv'
            if args.source:
                source = Path(args.source)
                if sha256(source) != SOURCE['sha256']:
                    raise SitError('{}: checksum {} instead of {}'.format(source, sha256(source), SOURCE['sha256']))
                origin = str(source.resolve())
            else:
                origin = download(source)
            n = convert(source, folder / DB_NAME, origin)
        with open(folder / 'sit_database.json', 'w') as f:
            json.dump(dict(SOURCE, origin=origin, patterns=n, spoligotyper=__version__,
                           database_sha256=sha256(folder / DB_NAME),
                           date=datetime.now(timezone.utc).isoformat(timespec='seconds')), f, indent=2)
    except (SitError, OSError) as e:
        log.error('%s', e)
        sys.exit(1)
    log.info('%d patterns saved in %s', n, folder / DB_NAME)
    log.info('This list is licensed under %s. Please cite: %s', SOURCE['license'], SOURCE['citation'])
