"""Sublineages of M. tuberculosis lineage 1 (Netikul et al. 2022)."""

from spoligotyper import l1, snp_groups
from spoligotyper.lineage import LineageCall
from spoligotyper.pipeline import Result, check_species
from spoligotyper.species import SpeciesCheck

BARCODE = snp_groups.read_barcode(str(l1.BARCODE))


def counts_for(groups, reads=20, n=None):
    """Seal counts with the derived allele of the first n SNPs of each group (all by default), ancestral otherwise."""
    counts, seen = {}, {}
    for row in BARCODE:
        key = '{}|{}|'.format(row['group'], row['position'])
        seen[row['group']] = seen.get(row['group'], 0) + 1
        derived = row['group'] in groups and (n is None or seen[row['group']] <= n)
        counts[key + 'derived'], counts[key + 'ancestral'] = (reads, 0) if derived else (0, reads)
    return counts


def test_barcode():
    groups = {row['group'] for row in BARCODE}
    assert len(BARCODE) == 1835 and groups == set(l1.NAMES) and len(groups) == 32
    assert min(sum(r['group'] == g for r in BARCODE) for g in groups) == 4
    assert sum(r['barcode'] == 'yes' for r in BARCODE) == 125
    assert l1.parent('L1.1.1.10') == 'L1.1.1' and l1.parent('L1.3') is None and l1.parent('L1.2.2.5') == 'L1.2.2'


def test_call():
    call = l1.call_l1(counts_for({'L1.3', 'L1.3.2'}), 'fastq')
    assert (call.lineage, call.name, call.main) == ('L1.3.2', 'typical spoligotypes: EAI1-SOM', 'L1.3')
    call = l1.call_l1(counts_for({'L1.1', 'L1.1.2', 'L1.1.2.2'}, reads=1), 'fasta')
    assert call.lineage == 'L1.1.2.2' and call.unsupported == []
    assert l1.call_l1(counts_for(set()), 'fastq').lineage == ''


def test_half_of_the_snps():
    """L1.1.3.2 has 224 SNPs: a strain sharing 2 of them (M. canettii) is not called; a real one carries most."""
    parents = {'L1.1', 'L1.1.3'}
    assert l1.call_l1(counts_for(parents | {'L1.1.3.2'}, n=2), 'fastq').called == []
    counts = counts_for(parents)
    counts.update({k: v for k, v in counts_for({'L1.1.3.2'}, n=150).items() if k.startswith('L1.1.3.2|')})
    assert l1.call_l1(counts, 'fastq').lineage == 'L1.1.3.2'


def test_parent_required():
    """L1.2.1 has 4 SNPs: 2 stray derived alleles would be enough, but L1.2 must be called too."""
    assert l1.call_l1(counts_for({'L1.2.1'}, n=2), 'fastq').called == []
    assert l1.call_l1(counts_for({'L1.2', 'L1.2.1'}), 'fastq').lineage == 'L1.2.1'


def test_mix_of_sublineages():
    """L1.1.1.2 + L1.3.1 reads: a few SNPs with both alleles in each group (minority strain at low depth)."""
    counts = counts_for({'L1.3', 'L1.3.1'})
    for group, n in (('L1.1', 8), ('L1.1.1', 8), ('L1.1.1.2', 8)):
        rows = [r for r in BARCODE if r['group'] == group][:n]
        for row in rows:
            key = '{}|{}|'.format(row['group'], row['position'])
            counts[key + 'derived'], counts[key + 'ancestral'] = 4, 16
    call = l1.call_l1(counts, 'fastq')
    assert call.mixed and call.lineage.startswith('mixed: L1.1 20%')


def test_coll_lineage_conflict_warning():
    """L1 sublineage SNPs in a sample that the barcode of Coll et al. calls lineage 4."""
    r = Result('S', counts=[20] * 43, binary='1' * 43, file_type='fasta', data='assembly', min_count=1,
               species=SpeciesCheck(mtbc=True, species='M. tuberculosis'), lineage=LineageCall(called=['4', '4.9']),
               l1=l1.call_l1(counts_for({'L1.3', 'L1.3.2'}, reads=1), 'fasta'))
    check_species(r)
    assert any('lineage 1 sublineage SNPs (L1.3.2) but the lineage SNPs indicate lineage 4' in w for w in r.warnings)
    r.warnings, r.lineage = [], LineageCall()  # No lineage SNP at all
    check_species(r)
    assert any('L1.3.2) but no lineage 1 SNP of the lineage barcode' in w for w in r.warnings)
