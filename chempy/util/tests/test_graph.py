# -*- coding: utf-8 -*-

import os
import shutil
import tempfile
from chempy.chemistry import Reaction, ReactionSystem, Substance
from chempy.util.graph import rsys2dot, rsys2graph


def _get_rsys():
    r1 = Reaction({'A': 2}, {'B': 1}, param=3.0)
    A = Substance('A', latex_name='\\boldsymbol{A}')
    B = Substance('B', latex_name='\\boldsymbol{B}')
    rsys = ReactionSystem([r1], [A, B])
    return rsys


def test_rsys2dot():
    rsys = _get_rsys()
    assert list(map(str.strip, rsys2dot(rsys))) == [
        'digraph None{',
        '{',
        'node [label=r1 shape=diamond]',
        'r1',
        '}',
        '"A" -> "r1" [label ="2"];',
        '"r1" -> "B" [label =""];',
        '}'
    ]


def test_rsys2graph():
    rsys = _get_rsys()
    tempdir = tempfile.mkdtemp()
    try:
        rsys2graph(rsys, os.path.join(tempdir, 'out.png'))
        rsys2graph(rsys, os.path.join(tempdir, 'out.ps'))
        # rsys2graph(rsys, os.path.join(tempdir, 'out.tex'))
        # https://github.com/kjellmf/dot2tex/issues/48
    finally:
        shutil.rmtree(tempdir)
