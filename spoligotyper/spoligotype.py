"""Spoligotype codes: binary, octal and hexadecimal patterns, and SB numbers from the Mbovis.org database."""

import re
from importlib.resources import files

N_SPACERS = 43
HEX_BLOCKS = (7, 7, 7, 7, 8, 7)  # Spacers per hexadecimal block (6 blocks, 43 spacers)
# The spoligotype itself is the binary pattern and its octal and hexadecimal codes; SB numbers are names given to
# patterns by the Mbovis.org database. A pattern without SB number is only absent from that database.
NOT_FOUND = 'Not in Mbovis.org'
NOT_FOUND_CUSTOM = 'Not in database'  # With --db

BINARY_PATTERN = re.compile(r'[01]{%d}' % N_SPACERS)
OCTAL_PATTERN = re.compile(r'[0-7]{14}[01]')


class SpoligoError(Exception):
    pass


def data_file(name):
    """Path of a file shipped in the package's data folder."""
    return files('spoligotyper').joinpath('data', name)


SPACERS_FASTA = data_file('spoligo_spacers.fasta')
# Known variants of standard spacers, too different (more than 1 mismatch) to be found as the spacer itself, but
# detected by the spoligotyping membrane: e.g. spacer 3 of M. orygis, whose patterns are named with spacer 3 present
SPACER_VARIANTS = data_file('spacer_variants.fasta')
SPOLIGOTYPE_DB = data_file('spoligotype_db.txt')


def read_spacer_names(fasta=SPACERS_FASTA):
    """Names of the spacers, in the order of the fasta file (spacer01 to spacer43)."""
    with open(fasta) as f:
        names = [line[1:].split()[0] for line in f if line.startswith('>')]
    if len(names) != N_SPACERS or len(set(names)) != N_SPACERS:
        raise SpoligoError('{} must contain {} distinct spacers (found {})'.format(fasta, N_SPACERS, len(set(names))))
    return names


def read_spacer_variants(fasta=SPACER_VARIANTS):
    """
    {variant name: (spacer name, description)}, e.g. {"spacer03_v1": ("spacer03", "2 mismatches, found in M. orygis")},
    from fasta descriptions such as "spacer=spacer03 mismatches=2 found_in=M. orygis".
    """
    variants = {}
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                name, _, description = line[1:].strip().partition(' ')
                description, _, found_in = description.partition(' found_in=')  # found_in can contain spaces
                fields = dict(field.split('=', 1) for field in description.split())
                variants[name] = (fields['spacer'], '{} mismatches, found in {}'.format(fields['mismatches'], found_in))
    return variants


def add_variants(counts, variants):
    """
    Count the reads of known spacer variants as reads of their spacer: the spacer count becomes the higher of the two
    (a read between the spacer and its variant matches both). Returns {spacer: variant reads} for the spacers whose
    count comes from a variant.
    """
    used = {}
    for variant, (spacer, _) in variants.items():
        if counts.get(variant, 0) > counts.get(spacer, 0):
            counts[spacer] = counts[variant]
            used[spacer] = variant
    return used


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


def lookup(binary, database, not_found=NOT_FOUND):
    """Name (SB number) of a binary pattern in the database, or not_found."""
    return database.get(binary, not_found)


def closest(binary, database, max_distance=3, limit=3):
    """
    Closest patterns of the database, for a pattern that is not in it.

    :return: list of (SB number, [spacers that differ]), with the smallest number of differences (at most
             max_distance), up to limit patterns. Empty if the pattern is in the database or none is close enough.
    """
    if binary in database:
        return []
    best, found = max_distance + 1, []
    for pattern, name in database.items():
        differences = [i + 1 for i, (a, b) in enumerate(zip(binary, pattern, strict=True)) if a != b]
        if len(differences) < best:
            best, found = len(differences), []
        if len(differences) == best:
            found.append((name, differences))
    return sorted(found)[:limit] if best <= max_distance else []


def describe_closest(matches):
    """e.g. "SB0140 (spacer 7 differs); SB0265 (spacers 7, 13 differ)"."""
    return '; '.join('{} (spacer{} {} differ{})'.format(name, 's' if len(diff) > 1 else '', ', '.join(map(str, diff)),
                                                       '' if len(diff) > 1 else 's')
                     for name, diff in matches)
