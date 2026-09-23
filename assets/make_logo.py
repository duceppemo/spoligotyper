#!/usr/bin/env python3
"""
Build the logo files from logo_source.svg, whose text needs fonts that are rarely installed.

The text is converted to outlines so the logo looks the same everywhere (GitHub, PDF viewers):
  logo.svg       text as paths, for light backgrounds
  logo_dark.svg  same, with light text colours, for dark backgrounds
  ../spoligotyper/data/logo.png  bitmap used in the PDF report

Requirements (development only): pip install fonttools uharfbuzz cairosvg
Fonts: Sora (SIL Open Font License, https://github.com/google/fonts/tree/main/ofl/sora) for the name,
DejaVu Sans for the tagline. Usage: python make_logo.py Sora[wght].ttf
"""

import re
import sys
from pathlib import Path

import cairosvg
import uharfbuzz as hb
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

HERE = Path(__file__).resolve().parent
DEJAVU = Path('/usr/share/fonts/truetype/dejavu')
NAME_WEIGHT = 600  # Sora SemiBold
DARK_COLOURS = {'#1C2541': '#E6EDF3', '#56627A': '#9DA7B3'}


class Font:
    def __init__(self, path, weight=None):
        self.tt = TTFont(path)
        blob = hb.Blob.from_file_path(str(path))
        self.hb = hb.Font(hb.Face(blob))
        self.variations = {'wght': weight} if weight else {}
        if self.variations:
            self.hb.set_variations(self.variations)
        self.upem = self.tt['head'].unitsPerEm
        self.glyph_set = self.tt.getGlyphSet(location=self.variations or None)
        self.names = self.tt.getGlyphOrder()

    def outline(self, text, x, y, size, letter_spacing=0.0):
        """SVG path data of text drawn with its baseline starting at (x, y). Returns (path data, end x)."""
        buf = hb.Buffer()
        buf.add_str(text)
        buf.guess_segment_properties()
        hb.shape(self.hb, buf)
        scale = size / self.upem
        pen = SVGPathPen(self.glyph_set)
        for info, pos in zip(buf.glyph_infos, buf.glyph_positions, strict=True):
            gx, gy = x + pos.x_offset * scale, y - pos.y_offset * scale
            self.glyph_set[self.names[info.codepoint]].draw(TransformPen(pen, (scale, 0, 0, -scale, gx, gy)))
            x += pos.x_advance * scale + letter_spacing
        return pen.getCommands(), x


def build(sora_path):
    source = (HERE / 'logo_source.svg').read_text()
    source = re.sub(r'<metadata>.*?</metadata>', '', source, flags=re.S)
    source = source.replace(' xmlns:c2pa="http://c2pa.org/manifest"', '')

    name_font = Font(sora_path, NAME_WEIGHT)
    regular = Font(DEJAVU / 'DejaVuSans.ttf')
    italic = Font(DEJAVU / 'DejaVuSans-Oblique.ttf')

    name, _ = name_font.outline('spoligotyper', 500, 265, 138, letter_spacing=-2)
    tagline, x = [], 506.0
    for text, font in (('In silico', italic), (' spoligotyping of the ', regular), ('M. tuberculosis', italic),
                       (' complex', regular)):
        d, x = font.outline(text, x, 338, 29)
        tagline.append(d)
    paths = ('<path fill="#1C2541" d="{}"/>\n<path fill="#56627A" d="{}"/>'
             .format(name, ' '.join(tagline)))

    light = re.sub(r'<text.*</text>', lambda m: paths, source, flags=re.S)
    light = light.replace('<svg ', '<svg role="img" aria-label="spoligotyper logo" ', 1)
    dark = light
    for old, new in DARK_COLOURS.items():
        dark = dark.replace('fill="{}"'.format(old), 'fill="{}"'.format(new))

    (HERE / 'logo.svg').write_text(light)
    (HERE / 'logo_dark.svg').write_text(dark)
    cairosvg.svg2png(bytestring=light.encode(), write_to=str(HERE.parent / 'spoligotyper' / 'data' / 'logo.png'),
                     output_width=1200)


if __name__ == '__main__':
    build(sys.argv[1])
