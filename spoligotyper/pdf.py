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
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from . import DOI, __version__
from .seal import KMER_SIZE
from .spoligotype import N_SPACERS, NOT_FOUND, data_file, describe_closest

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
        ('TOPPADDING', (0, 0), (-1, -1), 2.5),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 2.5),
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
        return (bool(first.error), -len(members), not first.found, first.spoligotype if first.found else first.octal)

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
    header = [text(h, SMALL) for h in ('Sample', 'Spoligotype', 'Octal', 'Species', 'Lineage',
                                       'Pattern (spacers 1 to 43)', 'Status')]
    rows, shading = [header], []
    for row, (group, r) in enumerate(by_spoligotype(results), 1):
        if group % 2:
            shading.append(('BACKGROUND', (0, row), (-1, row), GROUP_SHADE))
        spoligotype = '\n'.join(x for x in (r.spoligotype, r.sit if r.sit.startswith('SIT') else '') if x)
        rows.append([text(r.sample, WRAP), text(spoligotype or '-', SMALL), text(r.octal or '-', MONO),
                     text(r.species.species if r.species else '-', SMALL),
                     text((r.lineage.lineage if r.lineage else '') or '-', WRAP),
                     pattern(r.binary, square=2.4) if not r.error else text('-', SMALL), status_label(r.status)])
    widths = [1.0 * inch, 0.8 * inch, 1.05 * inch, 1.15 * inch, 0.7 * inch, 1.95 * inch, WIDTH - 6.65 * inch]
    table = Table(rows, colWidths=widths, repeatRows=1)
    table.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (-1, 0), LIGHT)] + shading))
    story += [table, text('Samples with the same spoligotype are grouped; the sample pages follow the same order.',
                          SMALL)]

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


def sample_section(result, number=1, total=1):
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
    rows = [('Spoligotype', Paragraph('<b>{}</b>'.format(escape(result.spoligotype)), BODY)),
            ('Octal', text(result.octal, MONO)),
            ('Hexadecimal', text(result.hexadecimal, MONO)),
            ('Binary', text(result.binary, MONO)),
            ('Pattern', pattern(result.binary, square=7, numbers=True)),
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
    if result.closest:
        rows.insert(1, ('Closest patterns', describe_closest(result.closest)))
    if result.sit:
        sit = {'Orphan': 'no SIT (SITVIT2 pattern seen once: orphan)',
               NOT_FOUND: 'not in the SITVIT2 list'}.get(result.sit, result.sit)
        family = ' · SITVIT2 family {}'.format(result.sit_family) if result.sit_family else ''
        closest_sit = ' · closest: {}'.format(describe_closest(result.closest_sit)) if result.closest_sit else ''
        rows.insert(1 + bool(result.closest), ('SIT', sit + family + closest_sit))
    story += [key_values(rows)]
    if result.species:
        story.append(KeepTogether(species_block(result)))
    story += [KeepTogether([Paragraph('{} per spacer'.format(unit.capitalize()), H3), spacer_table(result),
                            text('Blue: present (count ≥ {}). Orange: called absent but seen in some {}.'.format(
                                result.min_count, unit), SMALL)])]
    if result.warnings:
        story += [Paragraph('Warnings', H3)] + [text('• ' + w, SMALL) for w in result.warnings]
    return [KeepTogether(story[:4])] + story[4:]


def species_block(result):
    check, call = result.species, result.lineage
    rows = [('Species', Paragraph('<b>{}</b>'.format(escape(check.species)), BODY))]
    if call.lineage:
        rows.append(('Lineage', '{}{}{}'.format(call.lineage, ' · {}'.format(call.name) if call.name else '',
                                              ' · typical spoligotypes: {}'.format(call.spoligotypes.replace(';', ', '))
                                              if call.spoligotypes else '')))
    elif check.mtbc:
        rows.append(('Lineage', 'no lineage SNP found (lineages 1 to 7 and animal lineages are not detected)'))
    unit = result.unit
    fraction = check.mtbc_fraction
    rows.append(('MTBC DNA', 'median {:g} {} per control region, {:.0f}% of the control regions found{}'.format(
        check.control_depth, plural(unit, check.control_depth), check.control_found * 100,
        '' if fraction is None else '; about {:.0f}% of the reads are MTBC'.format(fraction * 100))))
    story = [Paragraph('Species and lineage', H3), key_values(rows)]
    if check.regions:
        table = [[text(h, SMALL) for h in ('Region', 'Deleted in', 'Depth relative to MTBC control', 'Call')]]
        deleted_in = {'RD1': 'BCG, Dassie bacillus', 'RD4': 'M. bovis, BCG (some M. canettii)',
                      'RD7': 'M. africanum lineage 6, animal lineages',
                      'RD9': 'M. africanum (lineages 5, 6), animal lineages',
                      'RD12': 'M. bovis, BCG, M. caprae, M. orygis (some M. canettii)'}
        for region, (state, ratio) in check.regions.items():
            table.append([text(region, SMALL), text(deleted_in[region], SMALL), text('{:.2f}'.format(ratio), SMALL),
                          text(state, SMALL)])
        t = Table(table, colWidths=[0.7 * inch, 3.0 * inch, 2.0 * inch, WIDTH - 5.7 * inch])
        t.setStyle(TableStyle(GRID + [('BACKGROUND', (0, 0), (-1, 0), LIGHT)]))
        story += [Spacer(1, 4), t]
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


def run_section(run):
    params = run.parameters
    method = ('Spacers were counted with Seal (BBTools): each of the {n} spacers is searched as a single {k}-mer '
              'on both strands, allowing 1 mismatch (k={k}, hdist=1, rcomp=t, maskmiddle=f, ambiguous=all). A '
              'spacer is present when it is found in at least the minimum count of reads (contigs for assemblies). '
              'The binary pattern is converted to the octal code (Dale et al. 2001) and the hexadecimal code, and '
              'looked up in the spoligotype database for its SB number. '
              'Species: the read depth of regions of difference RD1, RD4, RD7, RD9 and RD12 (100 bp segments, same '
              'Seal parameters) is compared with the depth of MTBC-specific control regions; a region is deleted when '
              'its relative depth is at most 0.1, present when it is at least 0.5, and the species is read from the '
              'RD profile as in the RD PCR scheme, refined with the lineage SNPs and the spacers. The fraction of '
              'MTBC reads is the control depth divided by the depth expected from the number of bases. Lineage: '
              'reads carrying each '
              'allele of the 62 SNPs of the Coll et al. (2014) barcode are counted with exact 31-mers (k=31, '
              'hdist=0); a lineage is called when at least 80% of the reads (and at least 3, or 1 contig) carry its '
              'allele.').format(n=N_SPACERS, k=KMER_SIZE)
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
                   + [(name, text('{path}\nMD5 {md5}'.format(**info), WRAP))
                      for name, info in run.species_data.items()]
                   + [('SIT database', text('{path}\n{patterns:,} patterns, {sits:,} SITs · SHA-256 {sha256}\n'
                                            '{source}'.format(**run.sit_database), WRAP)
                       if run.sit_database else 'not installed (spoligotyper-download-sit)')]),
        Paragraph('Method', H3),
        text(method, SMALL),
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
        text('Bushnell B. BBTools. https://sourceforge.net/projects/bbmap/', SMALL),
        text('Duceppe M-O. spoligotyper {}: in silico spoligotyping of Mycobacterium tuberculosis complex genomes. '
             'Zenodo. https://doi.org/{}'.format(__version__, DOI), SMALL),
    ]


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
        story += sample_section(result, i, len(ordered))
    story += run_section(run)

    class Canvas(NumberedCanvas):
        footer = 'spoligotyper {} · Spoligotyping report · generated {} by {} on {}'.format(
            __version__, (run.finished or datetime.now().astimezone()).strftime('%Y-%m-%d %H:%M %Z').strip(), run.user,
            run.host)

    doc.build(story, canvasmaker=Canvas)

