"""
Hierarchical SNP barcodes: lineages and sublineages ("groups") each defined by several SNPs, whose derived allele
marks the group. Used for the livestock lineages (Zwyer et al. 2021) and the lineage 1 sublineages (Netikul et al.
2022).

Each barcode has a tsv (group, position, ancestral, derived, gene) and a fasta with, for each SNP, the 61 bp of H37Rv
centred on it with each allele ("<group>|<position>|ancestral" and "...|derived"). Seal counts the reads containing
an exact 31-mer of each sequence. As in the KvarQ test suite of Zwyer et al., a group is called when at least 2 of its
SNPs carry the derived allele, and is mixed when at least 2 of its SNPs have both alleles.
"""

import csv
from collections import Counter
from dataclasses import dataclass, field
from functools import cache

from .lineage import CALL_FRACTION, MIN_READS, MIXED_FRACTION, MIXED_MIN_READS, SnpCall, confirmed_mixed

MIN_SNPS = 2  # SNPs with the derived allele to call a group, or with both alleles to call it mixed
MAX_MIXED_GROUPS = 6  # Groups listed in the text of a mix


@dataclass
class GroupCall:
    lineage: str = ''  # Most specific group reported, e.g. "La1.8.1", or "" when none; "mixed: ..." for a mix
    name: str = ''  # e.g. the former name of the group ("Eu1") or its typical spoligotype families
    group: str = ''  # Group of the barcode, e.g. "La1.7.X-unk4" for lineage "La1.7.X"
    called: list = field(default_factory=list)  # Groups with at least 2 SNPs with the derived allele
    snps: list = field(default_factory=list)  # SnpCall for every SNP covered by reads (lineage = group)
    mixed: list = field(default_factory=list)  # SnpCall with both alleles
    conflict: bool = False  # Groups that cannot occur together
    unsupported: list = field(default_factory=list)  # Parents of the group reported whose covered SNPs are ancestral
    main: str = ''  # Root group (e.g. La1, L1.2) when all the groups called share it
    mixed_within: bool = False  # A mix of subgroups of one root group (e.g. La1.7.1 and La1.8.1)


@dataclass(frozen=True)
class Scheme:
    """
    :param groups: {group: (parent, lineage reported, name)}, in the order of the reports. A group reported as ""
        (e.g. La1_La2, shared by La1 and La2) is only used to check the others.
    :param shared: {unreported group: roots it goes with}, e.g. {"La1_La2": ("La1", "La2")}
    :param labels: names of groups in mixed and conflict texts, e.g. {"La1_La2": "La1/La2"}
    :param min_fraction: fraction of the covered SNPs of a group that must carry the derived allele, besides MIN_SNPS:
        for barcodes with many SNPs per group, which a distant strain can share by chance (homoplasy)
    :param mixed_fraction: the same for a mix (SNPs with both alleles): lower, as the minority strain of a mix has
        fewer reads, some SNPs with too few of them to count; homoplasy gives one allele, not both
    :param require_parent: a group is only called when its parent group is called too
    """
    barcode: object
    fasta: object
    groups: tuple  # ((group, (parent, lineage, name)), ...): a tuple, for a hashable frozen dataclass
    shared: tuple = ()
    labels: tuple = ()
    min_fraction: float = 0.0
    mixed_fraction: float = 0.0
    require_parent: bool = False

    @property
    def group_info(self):
        return dict(self.groups)

    def label(self, group):
        return dict(self.labels).get(group, group)

    def ancestors(self, group):
        """The group and its parents, e.g. La1.8.1, La1.8, La1."""
        info, chain = self.group_info, []
        while group:
            chain.append(group)
            group = info[group][0]
        return chain

    def root(self, group):
        return '' if group in dict(self.shared) else self.ancestors(group)[-1]

    def on_one_path(self, groups):
        """True when the groups can occur in one strain: each is an ancestor of the most specific one (a shared
        group goes with its roots)."""
        groups = set(groups)
        shared = [(group, roots) for group, roots in self.shared if group in groups]
        groups -= {group for group, _ in shared}
        if any({self.root(g) for g in groups} - set(roots) for _, roots in shared):
            return False
        if not groups:
            return True
        deepest = max(groups, key=lambda g: len(self.ancestors(g)))
        return groups <= set(self.ancestors(deepest))

    def mixed_summary(self, mixed):
        """e.g. "La1 35%, La1.8 36%": the derived allele fraction of each group, over its SNPs with both alleles."""
        pooled = {}
        for s in mixed:
            pooled.setdefault(self.label(s.lineage), []).append(s)
        text = ['{} {:.0f}%'.format(name, 100 * sum(s.lineage_reads for s in snps) / sum(s.reads for s in snps))
                for name, snps in pooled.items()]
        return ', '.join(text[:MAX_MIXED_GROUPS]) + (', ...' if len(text) > MAX_MIXED_GROUPS else '')

    def call(self, counts, file_type, contaminated=False, barcode=None):
        """
        :param counts: {"<group>|<position>|ancestral" or "...|derived": reads} from Seal
        :param contaminated: the sample contains non-MTBC DNA: SNPs conserved in NTM are not used
        :param barcode: rows of the barcode (group, position, gene), by default those of the scheme's tsv
        :return: GroupCall
        """
        info = self.group_info
        conserved = ntm_conserved_snps(str(self.fasta))
        minimum = 1 if file_type == 'fasta' else MIN_READS
        result = GroupCall()
        positives, covered = {}, {}
        for row in barcode or read_barcode(str(self.barcode)):
            key = '{}|{}'.format(row['group'], row['position'])
            snp = SnpCall(row['group'], int(row['position']), row['gene'], counts.get(key + '|derived', 0),
                          counts.get(key + '|ancestral', 0))
            if snp.reads == 0:
                continue
            result.snps.append(snp)
            if key in conserved and contaminated:
                continue
            if snp.reads >= minimum:
                covered[snp.lineage] = covered.get(snp.lineage, 0) + 1
                if snp.fraction >= CALL_FRACTION:
                    positives[snp.lineage] = positives.get(snp.lineage, 0) + 1
            if file_type == 'fastq' and key not in conserved and snp.reads >= MIXED_MIN_READS and \
                    MIXED_FRACTION[0] <= snp.fraction <= MIXED_FRACTION[1]:
                result.mixed.append(snp)
        def enough(n, group, fraction):
            return n >= MIN_SNPS and n >= fraction * covered.get(group, 0)
        result.called = [group for group in info if enough(positives.get(group, 0), group, self.min_fraction)]
        if self.require_parent:
            result.called = [g for g in result.called if all(p in result.called for p in self.ancestors(g)[1:])]
        used = [s for s in result.snps if not (contaminated and '{}|{}'.format(s.lineage, s.position) in conserved)]
        mixed = confirmed_mixed(result.mixed, used, lambda group: self.ancestors(group)[1:])
        per_group = Counter(s.lineage for s in mixed)
        result.mixed = [s for s in mixed if enough(per_group[s.lineage], s.lineage, self.mixed_fraction)]
        if result.mixed:
            result.lineage = 'mixed: ' + self.mixed_summary(result.mixed)
        reported = [g for g in result.called if info[g][1]]
        result.conflict = bool(reported) and not self.on_one_path(result.called)
        roots = {self.root(g) for g in result.called} - {''}
        result.main = roots.pop() if len(roots) == 1 and not result.conflict else ''
        result.mixed_within = bool(result.mixed and result.main and all(info[s.lineage][0] for s in result.mixed))
        if not reported or result.mixed:
            return result
        if result.conflict:
            result.lineage = 'mixed: ' + ', '.join(self.label(g) for g in result.called)
        else:
            result.group = max(reported, key=lambda g: len(self.ancestors(g)))
            _, result.lineage, result.name = info[result.group]
            # A parent whose SNPs are covered but all ancestral contradicts the group reported
            result.unsupported = [g for g in self.ancestors(result.group)[1:] if positives.get(g, 0) == 0 and sum(
                s.lineage == g and s.reads >= minimum for s in used) >= MIN_SNPS]
        return result


@cache
def read_barcode(path):
    with open(path) as f:
        return tuple(csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t'))


@cache
def ntm_conserved_snps(path):
    """SNPs ("<group>|<position>") with an allele whose sequence is also found in NTM genomes."""
    with open(path) as f:
        return frozenset(line[1:].split(' ')[0].rsplit('|', 1)[0] for line in f
                         if line.startswith('>') and ' ntm=' in line)
