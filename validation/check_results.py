#!/usr/bin/env python3
"""Compare the spoligotyper results with the expected values and print a Markdown report."""

import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BEIJING = '0' * 34 + '1' * 9
# Documented SITs: H37Rv SIT451, BCG SIT482, and SIT1 for the Beijing strain
DOCUMENTED_SIT = {'H37Rv': 'SIT451', 'BCG_Pasteur': 'SIT482', 'CCDC5079': 'SIT1'}

# Simulated and public reads: expected species, lineage (prefix), octal, warning and livestock lineage (* not checked)
READS = {
    'ERR1744454': ('M. bovis', 'BOV', '664073777777600', '', 'La1.8.1'),
    # H37Rv lab stock: spacer 40 lost by most cells, spacers 41-43 kept by a minority (recombination between direct
    # repeats during passage; SIT1647, T-H37Rv family). The reference genome has spacers 40-43 (SIT451)
    'SRR12006063': ('M. tuberculosis', '4.9', '777777477760731', 'mixed sample', ''),
    'ERR027297': ('M. microti', 'BOV_AFRI', '000000000000600', '', ''),  # M. microti Maus IV, 2010 GAII reads
    # M. orygis 51145: 8 reads with a 1-SNP variant of spacer 3, not in the PacBio assembly (minority population or
    # cross-contamination): spacer 3 is called present, and flagged. Octal not checked
    'SRR16643349': ('M. orygis', 'BOV', '', 'mixed sample', 'La3'),
    'ERR2383628': ('M. africanum (lineage 6)', '6', '770777777777671', '', ''),  # M. africanum RB30001
    'SRR18636082': ('M. canettii', '', '000000000000000', 'as usual for M. canettii', ''),  # M. canettii ET1291
    'SRR23035463': ('M. canettii', '', '000000000000000', 'long reads', ''),  # ET1291, nanopore
    'sim_H37Rv_30x': ('M. tuberculosis', '4.9', '777777477760771', '', ''),
    'sim_H37Rv_30x_PE': ('M. tuberculosis', '4.9', '777777477760771', '', ''),
    'sim_H37Rv_10x': ('M. tuberculosis', '4.9', '', 'depth', ''),  # Low depth: warning expected
    'sim_Beijing_30x': ('M. tuberculosis', '2.2.1', 'Beijing', '', ''),
    'sim_AF2122_97_30x': ('M. bovis', 'BOV', '664073777777600', '', 'La1.8.1'),
    'sim_mixed_H37Rv70_AF2122_30': ('MTBC, mixed sample?', 'mixed', '', 'mixed sample', '*'),
    'sim_contaminated_H37Rv15x_marinum15x': ('M. tuberculosis', '4.9', '777777477760771', 'contamination', ''),
}


def read_table(path):
    with open(path) as f:
        return {row['Sample']: row for row in csv.DictReader(f, delimiter='\t')}


def read_json(path):
    with open(path) as f:
        return {sample['sample']: sample for sample in json.load(f)['samples']}


def la_support(sample):
    """e.g. "La1 4/4, La1.8 4/4, La1.8.1 4/4": SNPs with the derived allele of each group called."""
    if not sample or not sample.get('livestock'):
        return '-'
    snps = sample['livestock']['snps']
    return ', '.join('{} {}/{}'.format(group, sum(s['lineage'] == group and s['lineage_reads'] >= 0.8 * (
        s['lineage_reads'] + s['other_reads']) for s in snps), sum(s['lineage'] == group for s in snps))
        for group in sample['livestock']['called'] if group != 'La1_La2') or '-'


def matches(value, expected):
    """Expected values ending with * are prefixes; * alone is not checked."""
    if expected.endswith('*'):
        return value.startswith(expected[:-1])
    return value == expected


def octal_ok(row, expected):
    if not expected:
        return True
    if expected == 'Beijing':
        return row['Binary'] == BEIJING
    return row['Octal'] == expected


def main():
    failures = 0
    lines = ['# Validation', '', '## Reference genomes', '',
             '| Sample | Organism | SB | SIT (family) | Octal | Species (RD1 RD4 RD7 RD9 RD12) | Lineage '
             '| Expected lineage | La lineage | Result |',
             '|---|---|---|---|---|---|---|---|---|---|']
    genomes = read_table(HERE / 'results' / 'genomes' / 'spoligotyping.tsv')
    with open(HERE / 'genomes.tsv') as f:
        expected = list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))
    for exp in expected:
        if exp['sample'] not in genomes:
            failures += 1
            lines.append('| {} | {} | not typed | | | | | | | **FAIL** |'.format(exp['sample'], exp['organism']))
            continue
        row = genomes[exp['sample']]
        # The lineage can be more specific than expected (a sublineage), but not different
        lineage_ok = matches(row['Lineage'], exp['expected_lineage']) or (
            exp['expected_lineage'] not in ('', '*') and row['Lineage'].startswith(exp['expected_lineage'] + '.'))
        sit_ok = row['SIT'] == DOCUMENTED_SIT.get(exp['sample'], row['SIT'])
        ok = (matches(row['Species'], exp['expected_species']) and lineage_ok and octal_ok(row, exp['expected_octal'])
              and sit_ok and row['LaLineage'] == exp['expected_la'])
        failures += not ok
        lines.append('| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |'.format(
            exp['sample'], exp['organism'], row['SB'],
            '{} ({})'.format(row['SIT'], row['SITVIT2family']) if row['SITVIT2family'] else row['SIT'] or '-',
            row['Octal'],
            '{} ({})'.format(row['Species'], ' '.join({'present': '+', 'deleted': '-', 'partial': 'p',
                                                        'reduced': 'r'}.get(row[r], row[r] or '?')
                                                       for r in ('RD1', 'RD4', 'RD7', 'RD9', 'RD12'))),
            row['Lineage'] or '-',
            {'': '-', '*': 'not documented'}.get(exp['expected_lineage'], exp['expected_lineage']),
            row['LaLineage'] or '-', 'OK' if ok else '**FAIL**'))
    lines += ['', '## Reads', '', '| Sample | SB | Octal | Species | Lineage | La lineage | MTBC fraction | Warnings | '
              'Result |', '|---|---|---|---|---|---|---|---|---|']
    reads = read_table(HERE / 'results' / 'reads' / 'spoligotyping.tsv')
    for sample, (species, lineage, octal, warning, la) in READS.items():
        if sample not in reads:
            failures += 1
            lines.append('| {} | not typed | | | | | | | **FAIL** |'.format(sample))
            continue
        row = reads[sample]
        ok = (row['Species'] == species and row['Lineage'].startswith(lineage) and octal_ok(row, octal)
              and (warning in row['Warnings'] if warning else row['Status'] == 'ok') and matches(row['LaLineage'], la))
        failures += not ok
        lines.append('| {} | {} | {} | {} | {} | {} | {} | {} | {} |'.format(
            sample, row['SB'], row['Octal'], row['Species'], row['Lineage'] or '-', row['LaLineage'] or '-',
            row['MTBCFraction'], row['Warnings'].replace('|', '/') or '-', 'OK' if ok else '**FAIL**'))
    # Livestock lineages: the sublineage and SB number given by Zwyer et al. 2021 for each run
    lines += ['', '## Livestock lineages', '', '| Run | Origin | Expected | La lineage | SNPs | Species | SB (paper) | '
              'SB | Result |', '|---|---|---|---|---|---|---|---|---|']
    with open(HERE / 'livestock_reads.tsv') as f:
        la_expected = list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))
    la_rows = read_table(HERE / 'results' / 'la_reads' / 'spoligotyping.tsv')
    la_json = read_json(HERE / 'results' / 'la_reads' / 'spoligotyping.json')
    for exp in la_expected:
        if exp['run'] not in la_rows:
            failures += 1
            lines.append('| {} | not typed | | | | | | | **FAIL** |'.format(exp['run']))
            continue
        row = la_rows[exp['run']]
        species = {'La2': 'M. caprae', 'La3': 'M. orygis'}.get(exp['expected_la'], 'M. bovis')
        ok = row['LaLineage'] == exp['expected_la'] and row['Species'] == species
        failures += not ok
        lines.append('| {} | {} | {} | {} | {} | {} | {} | {} | {} |'.format(
            exp['run'], '{}, {}'.format(exp['country'], exp['host']), exp['expected_la'], row['LaLineage'] or '-',
            la_support(la_json.get(exp['run'])), row['Species'], exp['expected_sb'], row['SB'],
            'OK' if ok else '**FAIL**'))
    total = len(expected) + len(READS) + len(la_expected)
    lines += ['', '**{} of {} checks passed.**'.format(total - failures, total)]
    print('\n'.join(lines))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
