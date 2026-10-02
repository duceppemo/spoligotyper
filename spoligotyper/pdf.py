"""PDF report: results, the evidence behind each call, and everything needed to trace how they were produced."""

from datetime import datetime
from xml.sax.saxutils import escape

from reportlab.graphics.shapes import Drawing, Rect, String
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import (
    Image,
    KeepInFrame,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from . import DOI, __version__, l1, lineage, livestock, sitdb, snp_groups, species
from .seal import KMER_SIZE
from .spoligotype import N_SPACERS, data_file, describe_closest

# Colours of the logo
NAVY = colors.HexColor('#1C2541')
BLUE = colors.HexColor('#2B5FA8')
PALE_BLUE = colors.HexColor('#9DB3D4')
MAGENTA = colors.HexColor('#C8215F')
GREY = colors.HexColor('#56627A')
LIGHT = colors.HexColor('#EEF2F8')
GROUP_SHADE = colors.HexColor('#F0F3F8')  # Every other group of samples with the same spoligotype
AMBER = colors.HexColor('#F6D8A8')
STATUS_COLOURS = {'ok': colors.HexColor('#2E7D32'), 'warning': colors.HexColor('#B26A00'),
                  'failed': colors.HexColor('#C62828')}

PAGE_WIDTH, PAGE_HEIGHT = letter
MARGIN = 0.6 * inch
WIDTH = PAGE_WIDTH - 2 * MARGIN
LOGO_ASPECT = 500 / 1560

_styles = getSampleStyleSheet()
BODY = ParagraphStyle('body', parent=_styles['BodyText'], fontSize=8.5, leading=11, textColor=NAVY)
SMALL = ParagraphStyle('small', parent=BODY, fontSize=7.5, leading=9.5)
WRAP = ParagraphStyle('wrap', parent=SMALL, wordWrap='CJK')  # Breaks anywhere: long paths and checksums
MONO = ParagraphStyle('mono', parent=SMALL, fontName='Courier', wordWrap='CJK')
H1 = ParagraphStyle('h1', parent=_styles['Heading1'], fontSize=18, leading=22, textColor=NAVY, spaceAfter=2)
H2 = ParagraphStyle('h2', parent=_styles['Heading2'], fontSize=12.5, leading=16, textColor=NAVY, spaceBefore=10,
                    spaceAfter=5)
H3 = ParagraphStyle('h3', parent=_styles['Heading3'], fontSize=9.5, leading=12, textColor=BLUE, spaceBefore=6,
                    spaceAfter=3)
SUBTITLE = ParagraphStyle('subtitle', parent=BODY, textColor=GREY, fontSize=9)

GRID = [('GRID', (0, 0), (-1, -1), 0.4, PALE_BLUE),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('TOPPADDING', (0, 0), (-1, -1), 1.5),  # Compact rows: each sample fits on one page
        ('BOTTOMPADDING', (0, 0), (-1, -1), 1.5),
        ('LEFTPADDING', (0, 0), (-1, -1), 4),
        ('RIGHTPADDING', (0, 0), (-1, -1), 4)]


def text(value, style=BODY):
    return Paragraph(escape(str(value)).replace('\n', '<br/>'), style)


def plural(unit, count):
    """"reads" -> "read" when count is 1."""
    return unit[:-1] if count == 1 else unit


def timestamp(moment):
    return moment.strftime('%Y-%m-%d %H:%M:%S %Z').strip()


def human_size(size):
    for unit in ('B', 'KB', 'MB', 'GB'):
        if size < 1024 or unit == 'GB':
            return '{:.0f} {}'.format(size, unit) if unit == 'B' else '{:.1f} {}'.format(size, unit)
        size /= 1024


def _hex(colour):
    return '#' + colour.hexval()[2:]


def status_label(status):
    return Paragraph('<font color="{}"><b>{}</b></font>'.format(_hex(STATUS_COLOURS[status]), status.upper()), SMALL)


def pattern(binary, square=3.6, numbers=False):
    """The spoligotype as 43 squares: filled when the spacer is present."""
    gap = square * 0.25
    height = square + (square + 4 if numbers else 0)
    drawing = Drawing(N_SPACERS * (square + gap), height)
    for i, bit in enumerate(binary or '0' * N_SPACERS):
        x = i * (square + gap)
        drawing.add(Rect(x, height - square, square, square, strokeWidth=0.4, strokeColor=NAVY,
                         fillColor=NAVY if bit == '1' else colors.white))
        if numbers and (i == 0 or (i + 1) % 5 == 0 or i == N_SPACERS - 1):
            drawing.add(String(x + square / 2, 0, str(i + 1), fontName='Helvetica', fontSize=6, fillColor=GREY,
                               textAnchor='middle'))
    return drawing


def key_values(rows, key_width=1.45 * inch):
    table = Table([[text(k, SMALL), v if not isinstance(v, str) else text(v, SMALL)] for k, v in rows],
                  colWidths=[key_width, WIDTH - key_width])
    table.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (0, -1), LIGHT)]))
    return table


class NumberedCanvas(canvas.Canvas):
    """Canvas that knows the total number of pages, for "Page x of y" footers."""

    footer = ''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._pages = []

    def showPage(self):
        self._pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        total = len(self._pages)
        for page in self._pages:
            self.__dict__.update(page)
            self.setFont('Helvetica', 7)
            self.setFillColor(GREY)
            self.setStrokeColor(PALE_BLUE)
            self.line(MARGIN, 0.5 * inch, PAGE_WIDTH - MARGIN, 0.5 * inch)
            self.drawString(MARGIN, 0.36 * inch, self.footer)
            self.drawRightString(PAGE_WIDTH - MARGIN, 0.36 * inch, 'Page {} of {}'.format(self._pageNumber, total))
            super().showPage()
        super().save()


def by_spoligotype(results):
    """
    Samples grouped by spoligotype (identical 43-spacer pattern): the largest groups first, then SB numbers before
    patterns not in the database, and by sample name within a group. Failed samples come last.

    :return: [(group number, result)]
    """
    groups = {}
    for r in results:
        groups.setdefault('' if r.error else r.binary, []).append(r)

    def order(members):
        first = members[0]
        return (bool(first.error), -len(members), not first.found, first.sb if first.found else first.octal)

    return [(i, r) for i, members in enumerate(sorted(groups.values(), key=order))
            for r in sorted(members, key=lambda r: r.sample)]


def summary_section(results, run):
    counts = {s: sum(r.status == s for r in results) for s in ('ok', 'warning', 'failed')}
    story = [
        Image(str(data_file('logo.png')), width=3.2 * inch, height=3.2 * inch * LOGO_ASPECT, hAlign='LEFT'),
        Spacer(1, 6),
        Paragraph('Spoligotyping report', H1),
        text('{} · {} sample{} · {} ok, {} with warnings, {} failed · operator: {}'.format(
            run.started.strftime('%Y-%m-%d %H:%M %Z').strip(), len(results), 's' * (len(results) != 1), counts['ok'],
            counts['warning'], counts['failed'], run.operator), SUBTITLE),
        Paragraph('Summary', H2),
    ]
    default_db = run.database.get('default', True)
    names_label = 'SB / SIT' if default_db else 'Name / SIT'
    header = [text(h, SMALL) for h in ('Sample', 'Spoligotype (octal)', names_label, 'Species', 'Lineage',
                                       'Pattern (spacers 1 to 43)', 'Status')]
    rows, shading = [header], []
    for row, (group, r) in enumerate(by_spoligotype(results), 1):
        if group % 2:
            shading.append(('BACKGROUND', (0, row), (-1, row), GROUP_SHADE))
        names = '\n'.join(x for x in (r.sb if r.found else '', r.sit if r.sit.startswith('SIT') else '')
                           if x)
        rows.append([text(r.sample, WRAP), text(r.octal or '-', MONO), text(names or '-', SMALL),
                     text(r.species.species if r.species else '-', SMALL),
                     text((r.lineage.lineage if r.lineage else '') or '-', WRAP),
                     pattern(r.binary, square=2.4) if not r.error else text('-', SMALL), status_label(r.status)])
    widths = [1.0 * inch, 1.05 * inch, 0.8 * inch, 1.15 * inch, 0.7 * inch, 1.95 * inch, WIDTH - 6.65 * inch]
    table = Table(rows, colWidths=widths, repeatRows=1)
    table.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (-1, 0), LIGHT)] + shading))
    story += [table, text('Samples with the same spoligotype are grouped; the sample pages follow the same order. '
                          '{}: names of the pattern in the {} and SITVIT2 databases, "-" if the pattern has none. '
                          'See "Definitions and methods".'.format(names_label, 'Mbovis.org' if default_db else
                                                                  'spoligotype (--db)'), SMALL)]

    notes = [(r.sample, r.error.splitlines()[0] if r.error else w) for _, r in by_spoligotype(results)
             for w in ([r.error] if r.error else r.warnings)]
    if notes:
        story.append(Paragraph('Warnings and errors', H2))
        table = Table([[text(s, WRAP), text(n, SMALL)] for s, n in notes], colWidths=[1.55 * inch, WIDTH - 1.55 * inch])
        table.setStyle(TableStyle(GRID))
        story.append(table)

    review_block = [Paragraph('Review', H2),
                    text('Spoligotypes are called from whole genome sequencing data. The evidence for each call (reads '
                         'per spacer) is given in the sample sections, and the software, parameters and input files in '
                         'the run information section.', SMALL),
                    Spacer(1, 8)]
    review = Table([[text('Reviewed by', SMALL), '', text('Date', SMALL), '', text('Signature', SMALL), '']],
                   colWidths=[0.9 * inch, 1.9 * inch, 0.5 * inch, 1.2 * inch, 0.8 * inch, WIDTH - 5.3 * inch],
                   rowHeights=[0.35 * inch])
    review.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (0, 0), LIGHT), ('BACKGROUND', (2, 0), (2, 0), LIGHT),
                                       ('BACKGROUND', (4, 0), (4, 0), LIGHT)]))
    story.append(KeepTogether(review_block + [review]))  # Never the heading on one page and the box on the next
    return story


def spacer_table(result):
    """Spacer numbers and read counts, 22 then 21 spacers per line. Present spacers are shaded."""
    rows, style = [], list(GRID)
    for block, (start, end) in enumerate(((0, 22), (22, N_SPACERS))):
        numbers = [text('Spacer', SMALL)] + [text(i + 1, SMALL) for i in range(start, end)]
        counts = [text('Count', SMALL)] + [text(result.counts[i], SMALL) for i in range(start, end)]
        rows += [numbers + [''] * (23 - len(numbers)), counts + [''] * (23 - len(counts))]
        row = 2 * block + 1
        style.append(('BACKGROUND', (0, row - 1), (-1, row - 1), LIGHT))
        for i in range(start, end):
            column = i - start + 1
            if result.binary[i] == '1':
                style.append(('BACKGROUND', (column, row), (column, row), PALE_BLUE))
            elif result.counts[i] > 0:
                style.append(('BACKGROUND', (column, row), (column, row), AMBER))
    style += [('ALIGN', (1, 0), (-1, -1), 'CENTER'), ('LEFTPADDING', (0, 0), (-1, -1), 1.5),
              ('RIGHTPADDING', (0, 0), (-1, -1), 1.5), ('LINEBELOW', (0, 1), (-1, 1), 1.2, BLUE)]
    first = 0.5 * inch
    table = Table(rows, colWidths=[first] + [(WIDTH - first) / 22] * 22)
    table.setStyle(TableStyle(style))
    return table


def sample_section(result, number=1, total=1, db_label='SB number (Mbovis.org)'):
    story = [text('Sample {} of {}'.format(number, total), SUBTITLE),
             Paragraph('{} &nbsp; <font size="9" color="{}">{}</font>'.format(
                 escape(result.sample), _hex(STATUS_COLOURS[result.status]), result.status.upper()), H2)]
    files = [text('{}{}\n{} · modified {}{}'.format(f.path, '\n→ {}'.format(f.target) if f.target else '',
                                                    human_size(f.size), f.modified,
                                                    ' · MD5 {}'.format(f.md5) if f.md5 else ''), WRAP)
             for f in result.files]
    file_label = 'Input files' if len(result.files) > 1 else 'Input file'
    if result.error:
        story.append(key_values([(file_label, files or '-'), ('Error', text(result.error, MONO))]))
        return story

    kind = ('reads ({}, {})'.format(result.file_type, 'paired-end' if result.paired else 'single-end')
            if result.is_reads else 'assembly ({})'.format(result.file_type))
    unit = result.unit
    rows = [('Spoligotype (octal)', Paragraph('<b><font face="Courier">{}</font></b>'.format(result.octal), BODY)),
            ('Hexadecimal', text(result.hexadecimal, MONO)),
            ('Binary', text(result.binary, MONO)),
            ('Pattern', pattern(result.binary, square=7, numbers=True)),
            (db_label, result.sb if result.found else
             '{}: the pattern has no name in this database'.format(result.sb)),
            ('Data', kind),
            (file_label, files)]
    if result.reads is not None:
        rows.append(('Input size', '{:,} {}{}'.format(result.reads, plural(unit, result.reads),
                                                      ', {:,} bases'.format(result.bases) if result.bases else '')))
    if result.depth is not None:
        rows.append(('Estimated depth', '{:.0f}x (all bases / 4.4 Mb genome)'.format(result.depth)))
    rows += [('Minimum count', '{} {} per spacer to call it present'.format(result.min_count,
                                                                             plural(unit, result.min_count))),
             ('Present spacers', '{} of {}, median count {:g}'.format(result.binary.count('1'), N_SPACERS,
                                                                      result.median_present_count)),
             ('Run time', '{:.1f} s'.format(result.seconds))]
    position = 5  # After the database name
    if result.closest:
        rows.insert(position, ('Closest names', describe_closest(result.closest)))
        position += 1
    if result.sit:
        sb_note = (', and SITVIT2 lacks many patterns of the animal-adapted lineages: use the SB number'
                   if result.has_sb else '')
        sit = {sitdb.ORPHAN: 'Orphan: a SITVIT2 pattern without SIT',
               sitdb.NOT_FOUND: 'Not in SITVIT2 list: the pattern is not among the SITVIT2 patterns of the list'
               }.get(result.sit, result.sit) + (sb_note if result.sit in (sitdb.ORPHAN, sitdb.NOT_FOUND) else '')
        family = ' · SITVIT2 family {}'.format(result.sit_family) if result.sit_family else ''
        closest_sit = ' · closest: {}'.format(describe_closest(result.closest_sit)) if result.closest_sit else ''
        rows.insert(position, ('SIT (SITVIT2)', sit + family + closest_sit))
    story += [key_values(rows)]
    if result.species:
        story.append(KeepTogether(species_block(result)))
    variant_notes = [text('Spacer {}: counted from its known variant {} ({}): {} {} (standard sequence: {}).'.format(
        int(spacer[-2:]), variant, description, reads,
        plural(unit, reads), spacer_reads), SMALL)
        for spacer, (variant, reads, spacer_reads, description) in result.spacer_variants.items()]
    story += [KeepTogether([Paragraph('{} per spacer'.format(unit.capitalize()), H3), spacer_table(result),
                            text('Blue: present (count ≥ {}). Orange: called absent but seen in some {}.'.format(
                                result.min_count, unit), SMALL)] + variant_notes)]
    if result.warnings:
        story += [Paragraph('Warnings', H3)] + [text('• ' + w, SMALL) for w in result.warnings]
    return [KeepTogether(story[:4])] + story[4:]


DELETED_IN = {'RD1': 'BCG, Dassie bacillus', 'RD4': 'M. bovis, BCG (some M. canettii)',
              'RD7': 'M. africanum lineage 6, animal lineages',
              'RD9': 'M. africanum (lineages 5, 6), animal lineages',
              'RD12': 'M. bovis, BCG, M. caprae, M. orygis (some M. canettii)'}


def group_support(la, scheme):
    """
    e.g. " · SNPs with the derived allele: La1 4 of 4, La1.8 4 of 4", for the groups with derived alleles. Nothing for a
    mixed sample, whose lineage text gives the fraction of reads with each group's derived allele.
    """
    if la.mixed:
        return ''
    groups = {}
    for s in la.snps:
        groups.setdefault(s.lineage, []).append(s)
    support = ['{} {} of {}'.format(scheme.label(g), sum(s.fraction >= lineage.CALL_FRACTION for s in groups[g]),
                                    len(groups[g]))
               for g in scheme.group_info if any(s.fraction >= 0.1 for s in groups.get(g, []))]
    return ' · SNPs with the derived allele: ' + ', '.join(support) if support else ''


def species_block(result):
    check, call = result.species, result.lineage
    rows = [('Species', Paragraph('<b>{}</b>'.format(escape(check.species)), BODY))]
    if call.lineage:
        rows.append(('Lineage', '{}{}{}'.format(call.lineage, ' · {}'.format(call.name) if call.name else '',
                                              ' · typical spoligotypes: {}'.format(call.spoligotypes.replace(';', ', '))
                                              if call.spoligotypes else '')))
    elif check.mtbc:
        rows.append(('Lineage', 'no lineage SNP found (lineages 1 to 7 and animal lineages are not detected)'))
    la = result.livestock
    if la and la.lineage:
        rows.append(('Livestock lineage', '{}{} (Zwyer et al. 2021){}'.format(
            la.lineage, ' · {}'.format(la.name) if la.name else '', group_support(la, livestock.SCHEME))))
    sub = result.l1
    if sub and sub.lineage:
        rows.append(('L1 sublineage', '{}{} (Netikul et al. 2022){}'.format(
            sub.lineage, ' · {}'.format(sub.name) if sub.name else '', group_support(sub, l1.SCHEME))))
    unit = result.unit
    fraction = check.mtbc_fraction
    rows.append(('MTBC DNA', 'median {:g} {} per control region, {:.0f}% of the control regions found{}'.format(
        check.control_depth, plural(unit, check.control_depth), check.control_found * 100,
        '' if fraction is None else '; about {:.0f}% of the reads are MTBC'.format(fraction * 100))))
    story = [Paragraph('Species and lineage', H3), key_values(rows)]
    if check.regions:
        table = [[text(h, SMALL) for h in ('Region', 'H37Rv region', 'Usually deleted in', 'Segments found',
                                          'Depth vs control', 'Result')]]
        for region, region_call in check.regions.items():
            start, end = species.REGION_EXTENTS[region]
            # Short: the coordinates of a partial deletion are in the warnings, except for RD1mic (no warning)
            state = {species.PARTIAL: 'partially deleted', species.REDUCED: 'present at reduced depth: mixed sample?'
                     }.get(region_call.state, region_call.state)
            if region == 'RD1' and check.species == 'M. microti':
                state = region_call.describe() + ' (RD1mic of M. microti)'
            result_text = '{} {}'.format(region_call.sign.replace('-', '\u2212'), state)
            table.append([text(region, SMALL), text('{:,}-{:,}'.format(start, end), SMALL),
                          text(DELETED_IN[region], SMALL), text('{} of {}'.format(region_call.found, region_call.total),
                                                               SMALL),
                          text('{:.2f}'.format(region_call.ratio), SMALL), text(result_text, WRAP)])
        # One line per region (except the RD1mic of M. microti)
        t = Table(table, colWidths=[0.45 * inch, 1.2 * inch, 2.75 * inch, 0.65 * inch, 0.6 * inch, WIDTH - 5.65 * inch])
        t.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (-1, 0), LIGHT)]))
        story += [Spacer(1, 4), t,
                  text('+ / \u2212: DNA of the region present in / absent from the sample, from 100 bp segments inside '
                       'the region (not from amplicon sizes). See "Definitions and methods".', SMALL)]
    informative = [s for s in call.snps if s.fraction >= 0.1]  # Not the odd read with a sequencing error
    if informative or call.mixed:
        snps = sorted({s.position: s for s in informative + call.mixed}.values(), key=lambda s: s.position)
        table = [[text(h, SMALL) for h in ('Lineage SNP', 'Position (H37Rv)', 'Gene', 'Reads with lineage allele',
                                          'Other reads')]]
        for s in snps:
            table.append([text(s.lineage, SMALL), text('{:,}'.format(s.position), SMALL), text(s.locus, SMALL),
                          text(s.lineage_reads, SMALL), text(s.other_reads, SMALL)])
        t = Table(table, colWidths=[1.2 * inch, 1.3 * inch, 1.2 * inch, 1.8 * inch, WIDTH - 5.5 * inch], repeatRows=1)
        t.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (-1, 0), LIGHT)]))
        story += [Spacer(1, 4), t]
    return story


def definitions_section(run):
    """What the reported codes, names and regions of difference mean, and how they were determined."""
    sit = run.sit_database
    markers = species.segments()
    rd_rows = [[text(h, SMALL) for h in ('Region', 'H37Rv region (NC_000962.3)', 'Segments', 'Usually deleted in')]]
    for region in species.REGIONS:
        start, end = species.REGION_EXTENTS[region]
        n = sum(name.startswith(region + '_') for name in markers)
        rd_rows.append([text(region, SMALL), text('{:,}-{:,} ({:,} bp)'.format(start, end, end - start + 1), SMALL),
                        text(n, SMALL), text(DELETED_IN[region], SMALL)])
    rd_table = Table(rd_rows, colWidths=[0.6 * inch, 2.1 * inch, 0.7 * inch, WIDTH - 3.4 * inch])
    rd_table.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (-1, 0), LIGHT)]))
    if run.database.get('default', True):
        sb = ('The name of a pattern in the Mbovis.org database of spoligotypes of the RD9-deleted lineages '
              '(M. bovis and other animal-adapted lineages). "Not in Mbovis.org" only means that this database has no '
              'name for the pattern, as for all M. tuberculosis patterns: it does not mean that the spoligotype is '
              'invalid or new.')
    else:
        sb = ('The name of the pattern in the spoligotype database given with --db. "Not in database" only means '
              'that this database has no name for the pattern.')
    if sit:
        sit_text = ('Shared international type of the SITVIT2 database (Institut Pasteur de Guadeloupe), from the '
                    'SITVIT2 patterns published with SpolLineages ({:,} patterns, {:,} SITs, 2022 list). "Orphan": a '
                    'SITVIT2 pattern without SIT. "Not in SITVIT2 list": the pattern is not in this list, which does '
                    'not include the SITs created since 2022.'.format(sit['patterns'], sit['sits']))
        if sit.get('sb_patterns'):
            sit_text += (' SITVIT2 lacks many patterns of the animal-adapted lineages: only {:,} of the {:,} '
                         'Mbovis.org patterns have a SIT in this list. For these lineages, the SB number is the '
                         'reference name, and no closest SIT is given for a pattern with an SB number.'
                         ).format(sit['sb_with_sit'], sit['sb_patterns'])
    else:
        sit_text = 'Not reported: the SIT database was not installed for this run (spoligotyper-download-sit).'
    definitions = [
        ('Spoligotype',
         'The presence (1, filled square) or absence (0, empty square) of the 43 spacers of the direct repeat (DR) '
         'locus, written as a 43-digit binary pattern, a 15-digit octal code (Dale et al. 2001) or a hexadecimal '
         'code. These codes are universal: they identify the pattern itself, and are the ones to use to exchange '
         'spoligotypes.'),
        ('SB number', sb),
        ('SIT', sit_text),
        ('Spoligotype family',
         'SITVIT2 family of the pattern (e.g. Beijing, LAM3, BOV_1), and the spoligotype families typical of the '
         'lineage (Coll et al. 2014). Families are labels of patterns, not phylogenetic lineages.'),
        ('Lineage',
         'From the 62 SNPs of the barcode of Coll et al. (2014): the reads carrying each allele are counted with '
         'exact 31-mers; a lineage is called when at least 80% of the reads (at least 3, or 1 contig) carry its '
         'allele.'),
        ('Livestock lineage',
         'Lineage of the livestock-associated MTBC after Zwyer et al. (2021): La1 (M. bovis), La2 (M. caprae), La3 '
         '(M. orygis), and the La1 sublineages La1.1 to La1.8, with their former names (e.g. La1.8.1: clonal complex '
         'Eu1; BCG belongs to La1.2). From {} of the 89 marker SNPs of the extended data of the paper (4 or 5 '
         'per group), counted as above; a group is called when at least 2 of its SNPs carry the derived allele, as in '
         'the KvarQ test suite of the paper. La2 and La3 tell M. caprae from M. orygis, which have the same RD '
         'profile.'.format(len(livestock.read_barcode()))),
        ('L1 sublineage',
         'Sublineage of lineage 1 (East-African-Indian) after Netikul et al. (2022), in the revised nomenclature of '
         'lineage 1 (L1.1 to L1.3, down to e.g. L1.1.1.10), from the {:,} sublineage-specific SNPs of the paper (4 to '
         '224 per sublineage), counted as above; a sublineage is called when at least 2 of its SNPs, and at least '
         'half of those covered, carry the derived allele, and its parent group is called. The names of Coll et al. '
         'correspond to these groups: 1.2.1 is L1.2 (L1.2.1 and L1.2.2), 1.2.2 is L1.3, and 1.1.2 is L1.1.2.2 '
         '(L1.1.2.1 strains are 1.1).'.format(len(snp_groups.read_barcode(str(l1.BARCODE))))),
    ]
    rd_method = (
        'Regions of difference (RD) are determined in silico, from the sequencing data, without PCR. Each region is '
        'the part of the H37Rv genome missing from M. bovis AF2122/97 (RD4, RD7, RD9, RD12) or from BCG Pasteur '
        '(RD1), and is represented by 100 bp segments inside it (table below), chosen to be found in all the MTBC '
        'genomes that have the region and in none of those lacking it or of six non-tuberculous mycobacteria. A '
        'segment is found when Seal (25-mers, 1 mismatch) finds it in the reads at a depth of at least 5% of the '
        'depth of MTBC-specific control regions, or in the assembly. A region is present when all its segments are '
        'found (in reads, up to 10% of them, at least 1, may be missing: low depth); deleted when at most 10% of its '
        'segments (at least 1: stray reads) are found; partially deleted in between, with the H37Rv coordinates of '
        'the missing segments; present at reduced depth when, in reads, its segments are found at less than 50% of '
        'the control depth (mixed sample?).')
    rd_pcr = (
        'The + and \u2212 of the RD profile therefore mean that the DNA of the region is present or absent, like a '
        'PCR with primers inside the region (amplification = present). spoligotyper does not measure amplicon sizes '
        'or deletion junctions: assays that distinguish RDs by the size of an amplicon spanning the region may '
        'report partial or strain-specific deletions differently. A partially deleted region counts as + when at '
        'least half of its segments are found: for example, M. microti lacks the part of RD1 inside its own RD1mic '
        'deletion (reported as partially deleted) and is RD1 + in the RD profile, as in the classical RD PCR scheme. '
        'When the profile then matches no species, but would with a partially deleted region present, that region '
        'counts as present for the species (a deletion of one strain), as the warning says. '
        'The species is read from the RD profile, refined with the lineage SNPs and the spacers.')
    spacers = (
        'Spacers are counted with Seal (BBTools): each of the {n} spacers is searched as a single {k}-mer on both '
        'strands, allowing 1 mismatch (k={k}, hdist=1, rcomp=t, maskmiddle=f, ambiguous=all). A spacer is present when '
        'it is found in at least the minimum count of reads (contigs for assemblies). Known variants of standard '
        'spacers, too different to be found as the spacer (more than 1 mismatch) but named as present in the '
        'spoligotype databases, are searched the same way and counted for their spacer (e.g. spacer 3 of M. orygis, '
        '2 mismatches): the sample page says so. The fraction of MTBC reads is '
        'the depth of the control regions divided by the depth expected from the number of bases.'
    ).format(n=N_SPACERS, k=KMER_SIZE)
    return [
        PageBreak(),
        Paragraph('Definitions and methods', H2),
        key_values([(term, text(definition, SMALL)) for term, definition in definitions]),
        Paragraph('Regions of difference', H3),
        text(rd_method, SMALL), Spacer(1, 4), rd_table, Spacer(1, 4), text(rd_pcr, SMALL),
        Paragraph('Spacers and codes', H3),
        text(spacers, SMALL),
    ]


def run_section(run):
    params = run.parameters
    duration = (run.finished - run.started).total_seconds() if run.finished else 0
    return [
        PageBreak(),
        Paragraph('Run information', H2),
        key_values([('Operator', run.operator), ('User', run.user), ('Computer', run.host),
                    ('Operating system', run.system), ('Started', timestamp(run.started)),
                    ('Finished', timestamp(run.finished) if run.finished else '-'),
                    ('Duration', '{:.1f} s'.format(duration)), ('Working directory', text(run.working_directory, WRAP)),
                    ('Command', text(run.command, MONO))]),
        Paragraph('Parameters', H3),
        key_values([(k, str(v)) for k, v in params.items()]),
        Paragraph('Software', H3),
        key_values([(k, text(v, WRAP)) for k, v in run.software.items()]),
        Paragraph('Reference data', H3),
        key_values([('Spoligotype database', text('{path}\n{patterns:,} patterns · MD5 {md5}'.format(**run.database),
                                                  WRAP)),
                    ('Spacer sequences', text('{path}\n{spacers} spacers · MD5 {md5}'.format(**run.spacers), WRAP))]
                   + ([('Spacer variants', text('{path}\n{variants} known variant(s) · MD5 {md5}'.format(
                       **run.spacers['variants']), WRAP))] if run.spacers.get('variants') else [])
                   + [(name, text('{path}\nMD5 {md5}'.format(**info), WRAP))
                      for name, info in run.species_data.items()]
                   + [('SIT database', text('{path}\n{patterns:,} patterns, {sits:,} SITs · SHA-256 {sha256}\n'
                                            '{source}'.format(**run.sit_database), WRAP)
                       if run.sit_database else 'not installed (spoligotyper-download-sit)')]),
        Paragraph('References', H3),
        text('Kamerbeek J et al. Simultaneous detection and strain differentiation of Mycobacterium tuberculosis for '
             'diagnosis and epidemiology. J Clin Microbiol 35:907-914 (1997). '
             'doi:10.1128/jcm.35.4.907-914.1997', SMALL),
        text('Dale JW et al. Spacer oligonucleotide typing of bacteria of the Mycobacterium tuberculosis complex: '
             'recommendations for standardised nomenclature. Int J Tuberc Lung Dis 5:216-219 (2001).', SMALL),
        text('Smith NH, Upton P. Naming spoligotype patterns for the RD9-deleted lineage of the Mycobacterium '
             'tuberculosis complex; www.Mbovis.org. Infect Genet Evol 12:873-876 (2012). '
             'doi:10.1016/j.meegid.2011.08.002', SMALL),
        text('Brosch R et al. A new evolutionary scenario for the Mycobacterium tuberculosis complex. Proc Natl Acad '
             'Sci USA 99:3684-3689 (2002). doi:10.1073/pnas.052548299', SMALL),
        text('Couvin D, Segretier W, Stattner E, Rastogi N. Novel methods included in SpolLineages tool for fast and '
             'precise prediction of Mycobacterium tuberculosis complex spoligotype families. Database (Oxford) '
             '2020:baaa108. doi:10.1093/database/baaa108 (SIT database, from SITVIT2)', SMALL),
        text('Coll F et al. A robust SNP barcode for typing Mycobacterium tuberculosis complex strains. Nat Commun '
             '5:4812 (2014). doi:10.1038/ncomms5812', SMALL),
        text('Netikul T et al. Whole-genome single nucleotide variant phylogenetic analysis of Mycobacterium '
             'tuberculosis Lineage 1 in endemic regions of Asia and Africa. Sci Rep 12:1565 (2022). '
             'doi:10.1038/s41598-022-05524-0', SMALL),
        text('Zwyer M et al. A new nomenclature for the livestock-associated Mycobacterium tuberculosis complex based '
             'on phylogenomics. Open Res Eur 1:100 (2021). doi:10.12688/openreseurope.14029.2', SMALL),
        text('Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/', SMALL),
        text('Duceppe M-O. spoligotyper {}: in silico spoligotyping of Mycobacterium tuberculosis complex genomes. '
             'Zenodo. https://doi.org/{}'.format(__version__, DOI), SMALL),
    ]


def unwrap(flowables):
    """The flowables with their KeepTogether blocks unwrapped, for a KeepInFrame (which keeps them together)."""
    flat = []
    for f in flowables:
        flat += unwrap(f._content) if isinstance(f, KeepTogether) else [f]
    return flat


def write_pdf(results, run, path):
    """Write the PDF report of a run."""
    producer = 'spoligotyper {}'.format(__version__)
    doc = SimpleDocTemplate(str(path), pagesize=letter, leftMargin=MARGIN, rightMargin=MARGIN, topMargin=MARGIN,
                            bottomMargin=0.75 * inch, title='Spoligotyping report', author=run.operator,
                            subject=producer, creator=producer)
    story = summary_section(results, run)
    ordered = by_spoligotype(results)
    for i, (_, result) in enumerate(ordered, 1):
        story.append(PageBreak())  # One sample per page: its tables are never split
        section = sample_section(result, i, len(ordered), 'SB number (Mbovis.org)' if run.database.get('default', True)
                                 else 'Name (spoligotype database)')
        # A section too long for one page (long paths, many warnings) is scaled down slightly to fit it
        story.append(KeepInFrame(doc.width, doc.height - 12, unwrap(section), mode='shrink'))
    story += definitions_section(run)
    story += run_section(run)

    class Canvas(NumberedCanvas):
        footer = 'spoligotyper {} · Spoligotyping report · generated {} by {} on {}'.format(
            __version__, (run.finished or datetime.now().astimezone()).strftime('%Y-%m-%d %H:%M %Z').strip(), run.user,
            run.host)

    doc.build(story, canvasmaker=Canvas)

