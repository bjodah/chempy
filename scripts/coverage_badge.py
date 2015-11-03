#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script generates a "coverage" badge as a svg file from
the html report from coverage.py

Usage:

   $ ./coverage_badge.py htmlcov/ coverage.svg

"""

from __future__ import (absolute_import, division, print_function)
import os

# this template was generated from shields.io on 2015-10-11
template = """
<svg xmlns="http://www.w3.org/2000/svg" width="92" height="20">
<linearGradient id="b" x2="0" y2="100%">
<stop offset="0" stop-color="#bbb" stop-opacity=".1"/>
<stop offset="1" stop-opacity=".1"/>
</linearGradient>
<mask id="a">
<rect width="92" height="20" rx="3" fill="#fff"/>
</mask>
<g mask="url(#a)">
<path fill="#555" d="M0 0h63v20H0z"/>
<path fill="{0:s}" d="M63 0h29v20H63z"/>
<path fill="url(#b)" d="M0 0h92v20H0z"/>
</g>
<g fill="#fff" text-anchor="middle" font-family="DejaVu Sans,Verdana,Geneva,
    sans-serif" font-size="11">
<text x="31.5" y="15" fill="#010101" fill-opacity=".3">coverage</text>
<text x="31.5" y="14">coverage</text>
<text x="76.5" y="15" fill="#010101" fill-opacity=".3">{1:s}%</text>
<text x="76.5" y="14">{1:s}%</text>
</g>
</svg>
"""


def get_coverage(htmldir):
    for line in open(os.path.join(htmldir, 'index.html'), 'rt'):
        if 'pc_cov' in line:
            return int(line.split('pc_cov')[1].split(
                '>')[1].split('<')[0].rstrip('%'))
    raise ValueError("Could not find pc_cov in index.html")


def write_cov_badge_svg(path, percent):
    colors = '#e05d44 #fe7d37 #dfb317 #a4a61d #97CA00 #4c1'.split()
    limits_le = 50, 60, 70, 80, 90, 100
    c = next(clr for lim, clr in zip(limits_le, colors) if percent <= lim)
    with open(path, 'wt') as f:
        f.write(template.format(c, str(percent)))

if __name__ == '__main__':
    import sys
    assert len(sys.argv) == 3
    cov_percent = get_coverage(sys.argv[1])
    write_cov_badge_svg(sys.argv[2], cov_percent)
