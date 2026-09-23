"""Spoligotype codes: binary, octal and hexadecimal patterns, and SB numbers from the Mbovis.org database."""

import re
from importlib.resources import files

N_SPACERS = 43
HEX_BLOCKS = (7, 7, 7, 7, 8, 7)  # Spacers per hexadecimal block (6 blocks, 43 spacers)
NOT_FOUND = 'Spoligo not found'

BINARY_PATTERN = re.compile(r'[01]{%d}' % N_SPACERS)
OCTAL_PATTERN = re.compile(r'[0-7]{14}[01]')


class SpoligoError(Exception):
    pass


def data_file(name):
    """Path of a file shipped in the package's data folder."""
    return files('spoligotyper').joinpath('data', name)


SPACERS_FASTA = data_file('spoligo_spacers.fasta')
SPOLIGOTYPE_DB = data_file('spoligotype_db.txt')


def read_spacer_names(fasta=SPACERS_FASTA):
    """Names of the spacers, in the order of the fasta file (spacer01 to spacer43)."""
    with open(fasta) as f:
        names = [line[1:].split()[0] for line in f if line.startswith('>')]
    if len(names) != N_SPACERS or len(set(names)) != N_SPACERS:
        raise SpoligoError('{} must contain {} distinct spacers (found {})'.format(fasta, N_SPACERS, len(set(names))))
    return names


def check_binary(binary):
    if not BINARY_PATTERN.fullmatch(binary):
        raise SpoligoError('Invalid binary spoligotype "{}": expected {} digits, 0 or 1'.format(binary, N_SPACERS))


def to_binary(counts, spacer_names, min_count):
    """
    Binary spoligotype: "1" for each spacer seen at least min_count times, in spacer order.

    :param counts: {spacer name: number of reads (or contigs) matching it}. Missing spacers count as 0.
    """
    return ''.join('1' if counts.get(name, 0) >= min_count else '0' for name in spacer_names)


def binary_to_octal(binary):
    """
    Standard 15-digit octal code: spacers 1-42 by groups of 3, then spacer 43 alone.
    """
    check_binary(binary)
    return ''.join(str(int(binary[i:i + 3], 2)) for i in range(0, N_SPACERS, 3))


def binary_to_hex(binary):
    """Hexadecimal code: 6 blocks of 7, 7, 7, 7, 8 and 7 spacers, e.g. "6D-03-5F-7F-FF-60"."""
    check_binary(binary)
    blocks = []
    start = 0
    for size in HEX_BLOCKS:
        blocks.append('{:02X}'.format(int(binary[start:start + size], 2)))
        start += size
    return '-'.join(blocks)


def load_database(path=SPOLIGOTYPE_DB):
    """
    Read a spoligotype database: one pattern per line, "octal SB-number binary", separated by spaces or tabs.

    :return: {binary pattern: SB number}
    """
    database = {}
    with open(path) as f:
        for line_number, line in enumerate(f, 1):
            fields = line.split()
            if not fields or fields[0].startswith('#'):
                continue
            if len(fields) != 3:
                raise SpoligoError('{}, line {}: expected 3 columns (octal, SB number, binary), found {}'.format(
                    path, line_number, len(fields)))
            octal, name, binary = fields
            if not BINARY_PATTERN.fullmatch(binary) or binary_to_octal(binary) != octal:
                raise SpoligoError('{}, line {}: octal code "{}" does not match binary pattern "{}"'.format(
                    path, line_number, octal, binary))
            if binary in database and database[binary] != name:
                raise SpoligoError('{}, line {}: pattern {} is listed as both {} and {}'.format(
                    path, line_number, binary, database[binary], name))
            database[binary] = name
    if not database:
        raise SpoligoError('{} does not contain any spoligotype'.format(path))
    return database


def lookup(binary, database):
    """SB number of a binary pattern, or NOT_FOUND."""
    return database.get(binary, NOT_FOUND)
