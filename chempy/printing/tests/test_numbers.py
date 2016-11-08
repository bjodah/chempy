# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from chempy.util.testing import requires
from chempy.units import units_library, default_units as u

from ..numbers import (
    roman, _float_str_w_uncert, number_to_scientific_html, number_to_scientific_latex,
    number_to_scientific_unicode
)


def test__float_str_w_uncert():
    assert _float_str_w_uncert(-5739, 16.34, 2) == '-5739(16)'
    assert _float_str_w_uncert(-5739, 16.9, 2) == '-5739(17)'
    assert _float_str_w_uncert(0.0123, 0.00169, 2) == '0.0123(17)'
    assert _float_str_w_uncert(0.01234, 0.00169, 2) == '0.0123(17)'
    assert _float_str_w_uncert(0.01234, 0.0016501, 2) == '0.0123(17)'
    assert _float_str_w_uncert(0.01234, 0.0016501, 1) == '0.012(2)'
    assert _float_str_w_uncert(-9.99752e15, 349e10, 2) == '-9.9975(35)e15'
    assert _float_str_w_uncert(-9.99752e5, 349, 2) == '-999750(350)'
    assert _float_str_w_uncert(-9.99752e5, 349, 3) == '-999752(349)'
    assert _float_str_w_uncert(315, 17.9e-4, 2) == '315.0000(18)'


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
    assert number_to_scientific_latex(315, 17.9e-4, fmt=2) == '315.0000(18)'


@requires(units_library)
def test_number_to_scientific_latex__units():
    assert number_to_scientific_latex(315*u.km, 17.9*u.dm, fmt=2) == r'315.0000(18)\,\mathrm{km}'
    assert number_to_scientific_latex(315*u.km, 17.9*u.dm, u.m, fmt=2) == r'315000.0(18)\,\mathrm{m}'
    assert number_to_scientific_latex(1319*u.km, 41207*u.m, u.m, fmt=1) == r'1.32(4)\cdot 10^{6}\,\mathrm{m}'


def test_number_to_scientific_unicode():
    assert number_to_scientific_unicode(2e-17) == u'2·10⁻¹⁷'
    assert number_to_scientific_unicode(1e-17) == u'10⁻¹⁷'
