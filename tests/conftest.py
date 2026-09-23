"""Synthetic genomes and reads with a known spoligotype."""

import gzip
import random

import pytest

from spoligotyper.spoligotype import SPACERS_FASTA

DR = 'GTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGAC'  # 36 bp direct repeat between spacers
SB0140 = '1101101000001110111111111111111111111100000'  # M. bovis AF2122/97
H37RV = '1111111111111111111001111111111100001111111'  # M. tuberculosis H37Rv, not in the M. bovis database


def spacer_sequences():
    seqs = []
    with open(SPACERS_FASTA) as f:
        for line in f:
            if not line.startswith('>'):
                seqs.append(line.strip())
    return seqs


def genome(binary, rng, flank=3000):
    """Random sequence with a DR locus containing the spacers set to 1 in the binary pattern."""
    spacers = spacer_sequences()
    locus = ''.join(DR + s for s, bit in zip(spacers, binary, strict=True) if bit == '1') + DR
    return ''.join(rng.choice('ACGT') for _ in range(flank)) + locus + ''.join(rng.choice('ACGT') for _ in range(flank))


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def write_reads(path, seq, rng, n=1500, length=100, mate=None):
    """Random reads over seq, on both strands, gzipped fastq. Mates are reverse complemented."""
    with gzip.open(path, 'wt') as f:
        for i in range(n):
            p = rng.randint(0, len(seq) - length)
            read = seq[p:p + length]
            if rng.random() < 0.5 or mate == 2:
                read = reverse_complement(read)
            f.write('@r{}/{}\n{}\n+\n{}\n'.format(i, mate or 1, read, 'I' * length))


@pytest.fixture(scope='session')
def data(tmp_path_factory):
    rng = random.Random(1)
    folder = tmp_path_factory.mktemp('data')
    bovis = genome(SB0140, rng)
    (folder / 'AF2122.fasta').write_text('>chr\n{}\n'.format(bovis))
    (folder / 'H37Rv.fna').write_text('>chr\n{}\n'.format(genome(H37RV, rng)))
    write_reads(folder / 'bovis_R1.fastq.gz', bovis, rng, mate=1)
    write_reads(folder / 'bovis_R2.fastq.gz', bovis, rng, mate=2)
    write_reads(folder / 'bovis_single.fq.gz', bovis, rng)
    write_reads(folder / 'low_coverage.fastq.gz', bovis, rng, n=600)
    (folder / 'not_mtbc.fasta').write_text('>chr\n{}\n'.format(''.join(rng.choice('ACGT') for _ in range(5000))))
    return folder
