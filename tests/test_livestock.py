"""Lineages of the livestock-associated MTBC (Zwyer et al. 2021)."""

from spoligotyper import livestock, species

from .test_species_lineage import profile_counts

BARCODE = livestock.read_barcode()


def counts_for(groups, reads=20, other=0):
    """Seal counts with the derived allele of every SNP of the groups, and the ancestral allele of the others."""
    counts = {}
    for row in BARCODE:
        key = '{}|{}|'.format(row['group'], row['position'])
        derived = row['group'] in groups
        counts[key + 'derived'] = reads if derived else other
        counts[key + 'ancestral'] = other if derived else reads
    return counts


def test_barcode():
    groups = {row['group'] for row in BARCODE}
    assert groups == set(livestock.GROUPS)
    assert len(BARCODE) == 88 and all(sum(r['group'] == g for r in BARCODE) >= 4 for g in groups)
    assert all(r['ancestral'] != r['derived'] for r in BARCODE)


def test_most_specific_group():
    call = livestock.call_livestock(counts_for({'La1_La2', 'La1', 'La1.8', 'La1.8.1'}), 'fastq')
    assert (call.lineage, call.name, call.main, call.conflict) == ('La1.8.1', 'Eu1', 'La1', False)
    call = livestock.call_livestock(counts_for({'La1_La2', 'La1', 'La1.2', 'La1.2_BCG'}, reads=1), 'fasta')
    assert (call.lineage, call.name) == ('La1.2', 'BCG')
    call = livestock.call_livestock(counts_for({'La1', 'La1.7', 'La1.7.X-unk4'}), 'fastq')
    assert (call.lineage, call.name, call.group) == ('La1.7.X', 'unknown4', 'La1.7.X-unk4')
    call = livestock.call_livestock(counts_for({'La3'}), 'fastq')
    assert (call.lineage, call.name, call.main) == ('La3', 'M. orygis', 'La3')


def test_needs_two_snps():
    counts = counts_for({'La1_La2', 'La2'})
    la2 = [r for r in BARCODE if r['group'] == 'La2']
    for row in la2[1:]:  # Only one La2 SNP with the derived allele
        key = '{}|{}|'.format(row['group'], row['position'])
        counts[key + 'derived'], counts[key + 'ancestral'] = 0, 20
    call = livestock.call_livestock(counts, 'fastq')
    assert call.lineage == '' and call.called == ['La1_La2'] and call.main == ''


def test_no_call_for_other_lineages():
    call = livestock.call_livestock(counts_for(set()), 'fastq')
    assert (call.lineage, call.called, call.mixed) == ('', [], [])
    assert livestock.call_livestock({}, 'fastq').snps == []


def test_conflict_and_mixed():
    call = livestock.call_livestock(counts_for({'La1', 'La1.7', 'La1.7.1', 'La1.8', 'La1.8.1'}), 'fastq')
    assert call.conflict and call.lineage == 'mixed: La1, La1.7, La1.7.1, La1.8, La1.8.1' and call.main == ''
    call = livestock.call_livestock(counts_for({'La2', 'La3'}, reads=1), 'fasta')
    assert call.conflict and call.main == ''
    call = livestock.call_livestock(counts_for({'La1', 'La1.8'}, reads=12, other=8), 'fastq')  # 60% derived
    assert call.mixed and call.lineage.startswith('mixed: La1')


def test_contaminated_ignores_ntm_conserved_snps():
    conserved = livestock.ntm_conserved_snps()
    assert 'La1_La2|195566' in conserved and len(conserved) == 8
    counts = counts_for({'La1', 'La1.8'}, reads=12, other=8)
    call = livestock.call_livestock(counts, 'fastq', contaminated=True)
    assert not any('{}|{}'.format(s.lineage, s.position) in conserved for s in call.mixed)


def test_species_from_livestock_lineage():
    """La2 and La3 tell M. caprae from M. orygis, which share the RD profile."""
    for la, expected in (('La2', 'M. caprae'), ('La3', 'M. orygis'), ('', 'M. orygis or M. caprae')):
        check = species.check_species(profile_counts('++---'), 'fastq')
        species.name_species(check, ['BOV', 'BOV_AFRI'], livestock=la)
        assert check.species == expected
        assert species.consistency_warnings(check, ['BOV', 'BOV_AFRI'], la) == []
    check = species.check_species(profile_counts('+----'), 'fastq')
    species.name_species(check, ['BOV', 'BOV_AFRI'], livestock='La1')
    assert check.species == 'M. bovis' and species.consistency_warnings(check, ['BOV'], 'La1') == []
    assert species.consistency_warnings(check, ['BOV'], 'La3') == [
        'RD4 is deleted but the livestock lineage SNPs indicate La3 (M. orygis)']
    check = species.check_species(profile_counts('+++++'), 'fastq')
    assert species.consistency_warnings(check, [], 'La1') == [
        'RD9 is present but the livestock lineage SNPs indicate La1 (M. bovis)',
        'RD4 is present but the livestock lineage SNPs indicate La1 (M. bovis)']


def set_snps(counts, group, derived, ancestral, n=None):
    """Set the counts of the first n SNPs of a group (all by default)."""
    rows = [r for r in BARCODE if r['group'] == group][:n]
    for row in rows:
        key = '{}|{}|'.format(row['group'], row['position'])
        counts[key + 'derived'], counts[key + 'ancestral'] = derived, ancestral
    return counts


def test_conflict_text_names_the_groups():
    call = livestock.call_livestock(counts_for({'La1_La2', 'La3'}, reads=1), 'fasta')
    assert call.conflict and call.lineage == 'mixed: La1/La2, La3'
    call = livestock.call_livestock(counts_for({'La1', 'La1.7', 'La1.7.X-unk4', 'La1.7.X-unk5'}, reads=1), 'fasta')
    assert call.lineage == 'mixed: La1, La1.7, La1.7.X-unk4, La1.7.X-unk5'


def test_one_odd_site_is_not_a_mix():
    counts = set_snps(counts_for({'La1_La2', 'La1', 'La1.8', 'La1.8.1'}), 'La1.3', 6, 14, n=1)
    call = livestock.call_livestock(counts, 'fastq')
    assert call.mixed == [] and call.lineage == 'La1.8.1'
    call = livestock.call_livestock(set_snps(counts, 'La1.3', 6, 14, n=2), 'fastq')
    assert [s.lineage for s in call.mixed] == ['La1.3', 'La1.3']


def test_mix_of_sublineages_keeps_the_species():
    """La1.7.1 and La1.8.1 at 50% each: M. bovis (La1 in all the reads), mixed sample."""
    counts = counts_for({'La1_La2', 'La1'})
    for group in ('La1.7', 'La1.7.1', 'La1.8', 'La1.8.1'):
        set_snps(counts, group, 10, 10)
    call = livestock.call_livestock(counts, 'fastq')
    assert call.mixed_within and call.main == 'La1' and call.lineage.startswith('mixed: La1.7 50%')
    check = species.check_species(profile_counts('+----'), 'fastq')
    species.name_species(check, ['BOV', 'BOV_AFRI'], livestock=call.main, mixed_within=call.mixed_within)
    assert check.species == 'M. bovis, mixed sample?'


def test_parent_without_its_snps():
    counts = counts_for({'La1.8', 'La1.8.1'})  # La1 SNPs covered, all ancestral
    call = livestock.call_livestock(counts, 'fastq')
    assert call.lineage == 'La1.8.1' and call.unsupported == ['La1']
    assert livestock.call_livestock(counts_for({'La1', 'La1.8', 'La1.8.1'}), 'fastq').unsupported == []
