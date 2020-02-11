# -*- coding: utf-8 -*-
"""pytest lint tests for the tests."""

import pytest
import pycodestyle
import pydocstyle

import os

# Let's muck around and get the path to theses files.
dir = os.path.dirname(os.path.realpath(__file__))

files = [
    '_aqueous.py',
    ]


def test_pep8_conformance():
    """Test conformance to PEP 8."""
    style = pycodestyle.StyleGuide(quiet=True)
    for file in files:
        result = style.check_files([os.path.join(dir, file)])
        assert result.total_errors == 0


def test_pep257_conformance():
    """Test conformance to PEP 257."""
    for file in files:
        errors = pydocstyle.check([os.path.join(dir, file)])
        num = 0
        if errors:
            for error in errors:
                print(error)
                num += 1
        assert num == 0
