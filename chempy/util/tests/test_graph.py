# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


import os
import subprocess
import shutil
import tempfile

import pytest

from chempy import Reaction, ReactionSystem, Substance
from ..graph import rsys2dot, rsys2graph
from ..testing import requires

try:
    dot_missing = subprocess.call(['dot', '-?']) != 0
except OSError:
    dot_missing = True


def _get_rsys():
    r1 = Reaction({'A': 2}, {'B': 1}, param=3.0)
    A = Substance('A', latex_name='\\boldsymbol{A}')
    B = Substance('B', latex_name='\\boldsymbol{B}')
    rsys = ReactionSystem([r1], [A, B])
    return rsys


@requires('numpy')
@pytest.mark.skipif(dot_missing, reason='graphviz not installed? (dot command missing)')
def test_rsys2dot():
    rsys = _get_rsys()
    assert list(map(str.strip, rsys2dot(rsys))) == [
        'digraph "None" {',
        '"A" [fontcolor=maroon label="A"];',
        '"B" [fontcolor=darkgreen label="B"];',
        '{',
        'node [label="r1",shape=diamond]',
        'r1',
        '}',
        '"A" -> "r1" [color=maroon,fontcolor=maroon,label="2"];',
        '"r1" -> "B" [color=darkgreen,fontcolor=darkgreen,label=""];',
        '}'
    ]


@pytest.mark.slow
@requires('numpy')
@pytest.mark.skipif(dot_missing, reason='graphviz not installed? (dot command missing)')
def test_rsys2graph():
    rsys = _get_rsys()
    tempdir = tempfile.mkdtemp()
    try:
        rsys2graph(rsys, os.path.join(tempdir, 'out.png'))
        rsys2graph(rsys, os.path.join(tempdir, 'out.ps'))
        try:
            subprocess.call(['dot2tex', '-v'])
        except Exception:
            pass
        else:
            rsys2graph(rsys, os.path.join(tempdir, 'out.tex'))
    finally:
        shutil.rmtree(tempdir)
