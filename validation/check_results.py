#!/usr/bin/env python3
"""Compare the spoligotyper results with the expected values and print a Markdown report."""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BEIJING = '0' * 34 + '1' * 9

# Simulated and public reads: expected species, lineage (prefix), octal and warning
READS = {
    'ERR1744454': ('M. bovis', 'BOV', '664073777777600', ''),
    'sim_H37Rv_30x': ('M. tuberculosis', '4.9', '777777477760771', ''),
    'sim_H37Rv_30x_PE': ('M. tuberculosis', '4.9', '777777477760771', ''),
    'sim_H37Rv_10x': ('M. tuberculosis', '4.9', '', 'depth'),  # Low depth: warning expected
    'sim_Beijing_30x': ('M. tuberculosis', '2.2.1', 'Beijing', ''),
    'sim_AF2122_97_30x': ('M. bovis', 'BOV', '664073777777600', ''),
    'sim_mixed_H37Rv70_AF2122_30': ('MTBC, mixed sample?', 'mixed', '', 'mixed sample'),
    'sim_contaminated_H37Rv15x_marinum15x': ('M. tuberculosis', '4.9', '777777477760771', 'contamination'),
}


def read_table(path):
    with open(path) as f:
        return {row['Sample']: row for row in csv.DictReader(f, delimiter='\t')}


def octal_ok(row, expected):
    if not expected:
        return True
    if expected == 'Beijing':
        return row['Binary'] == BEIJING
    return row['Octal'] == expected


def main():
    failures = 0
    lines = ['# Validation', '', '## Reference genomes', '',
             '| Sample | Organism | Spoligotype | Octal | Species | Lineage | Expected lineage | Result |',
             '|---|---|---|---|---|---|---|---|']
    genomes = read_table(HERE / 'results' / 'genomes' / 'spoligotyping.tsv')
    with open(HERE / 'genomes.tsv') as f:
        expected = list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))
    for exp in expected:
        row = genomes[exp['sample']]
        # The lineage can be more specific than expected (a sublineage), but not different
        lineage_ok = row['Lineage'] == exp['expected_lineage'] or (
            exp['expected_lineage'] != '' and row['Lineage'].startswith(exp['expected_lineage'] + '.'))
        ok = row['Species'] == exp['expected_species'] and lineage_ok and octal_ok(row, exp['expected_octal'])
        failures += not ok
        lines.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(
            exp['sample'], exp['organism'], row['Spoligotype'], row['Octal'], row['Species'], row['Lineage'] or '-',
            exp['expected_lineage'] or '-', 'OK' if ok else '**FAIL**'))
    lines += ['', '## Reads', '', '| Sample | Spoligotype | Octal | Species | Lineage | MTBC fraction | Warnings | '
              'Result |', '|---|---|---|---|---|---|---|---|']
    reads = read_table(HERE / 'results' / 'reads' / 'spoligotyping.tsv')
    for sample, (species, lineage, octal, warning) in READS.items():
        row = reads[sample]
        ok = (row['Species'] == species and row['Lineage'].startswith(lineage) and octal_ok(row, octal)
              and (warning in row['Warnings'] if warning else row['Status'] == 'ok'))
        failures += not ok
        lines.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(
            sample, row['Spoligotype'], row['Octal'], row['Species'], row['Lineage'] or '-', row['MTBCFraction'],
            row['Warnings'].replace('|', '/') or '-', 'OK' if ok else '**FAIL**'))
    total = len(expected) + len(READS)
    lines += ['', '**{} of {} checks passed.**'.format(total - failures, total)]
    print('\n'.join(lines))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
