# -*- coding: utf-8 -*-

import shutil
import subprocess
import tempfile

import pytest


from chempy.util.table import (
    rsys2tablines, rsys2table, rsys2pdf_table
)
from .test_graph import _get_rsys


try:
    pdflatex_missing = subprocess.call(['pdflatex', '--version']) != 0
except OSError:
    pdflatex_missing = True


def test_rsys2tablines():
    assert rsys2tablines(_get_rsys(), tex=False) == [
        '1 & 2 A & -> & B & 3 & - & None'
    ]


def test_rsys2table():
    assert rsys2table(_get_rsys()) == (
        r"""
\begin{table}
\centering
\label{tab:none}
\caption[None]{None}
\begin{tabular}{lllllll}
\toprule
Id. & Reactants &  & Products & {Rate constant} & Unit & Ref \\
\midrule
1 & \ensuremath{2 \boldsymbol{A}} & \ensuremath{\rightarrow} &""" +
        r""" \ensuremath{\boldsymbol{B}} & \ensuremath{3} & \ensuremath{-} & None \\
\bottomrule
\end{tabular}
\end{table}""")


@pytest.mark.slow
@pytest.mark.parametrize('longtable', (True, False))
@pytest.mark.skipif(pdflatex_missing, reason='latex not installed? (pdflatex command missing)')
def test_rsys2pdf_table(longtable):
    rsys = _get_rsys()
    tempdir = tempfile.mkdtemp()
    try:
        rsys2pdf_table(rsys, tempdir, longtable=longtable)
    finally:
        shutil.rmtree(tempdir)


@pytest.mark.slow
@pytest.mark.skipif(pdflatex_missing, reason='latex not installed? (pdflatex command missing)')
def test_rsys2pdf_table_no_output_dir():
    rsys = _get_rsys()
    rsys2pdf_table(rsys, save=False)
