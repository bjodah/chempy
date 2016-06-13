# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from ..numbers import roman, number_to_scientific_html, number_to_scientific_latex, number_to_scientific_unicode


def test_roman():
    assert roman(4) == 'IV'
    assert roman(20) == 'XX'
    assert roman(94) == 'XCIV'
    assert roman(501) == 'DI'


def test_number_to_scientific_html():
    assert number_to_scientific_html(2e-17) == '2&sdot;10<sup>-17</sup>'
    assert number_to_scientific_html(1e-17) == '10<sup>-17</sup>'


def test_number_to_scientific_latex():
    assert number_to_scientific_latex(2e-17) == r'2\cdot 10^{-17}'
    assert number_to_scientific_latex(1e-17) == '10^{-17}'


def test_number_to_scientific_unicode():
    assert number_to_scientific_unicode(2e-17) == u'2·10⁻¹⁷'
    assert number_to_scientific_unicode(1e-17) == u'10⁻¹⁷'
