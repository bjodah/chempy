# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from math import log10, floor

from ..units import html_of_unit, latex_of_unit, unicode_of_unit, to_unitless, unit_of
from ..util.parsing import _unicode_sup


def roman(num):
    """
    Examples
    --------
    >>> roman(4)
    'IV'
    >>> roman(17)
    'XVII'

    """
    tokens = 'M CM D CD C XC L XL X IX V IV I'.split()
    values = 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1
    result = ''
    for t, v in zip(tokens, values):
        cnt = num//v
        result += t*cnt
        num -= v*cnt
    return result


def _mag(num):
    return int(floor(log10(abs(num))))


def _float_str_w_uncert(x, xe, precision=2):
    """ Prints uncertain number with parenthesis

    Parameters
    ----------
    x : nominal value
    xe : uncertainty
    precision : number of significant digits in uncertainty

    Examples
    --------
    >>> _float_str_w_uncert(-9.99752e5, 349, 3)
    '-999752(349)'
    >>> _float_str_w_uncert(-9.99752e15, 349e10, 2)
    '-9.9975(35)e15'
    >>> _float_str_w_uncert(3.1416, 0.029, 1)
    '3.14(3)'
    >>> _float_str_w_uncert(3.1416e9, 2.9e6, 1)
    '3.142(3)e9'

    Returns
    -------
    shortest string representation of "x +- xe" either as
    ``x.xx(ee)e+xx`` or ``xxx.xx(ee)``

    Notes
    -----
    The code in this function is from a question on StackOverflow:
        http://stackoverflow.com/questions/6671053
        written by:
            Lemming, http://stackoverflow.com/users/841562/lemming
        the code is licensed under 'CC-WIKI'.
        (see: http://blog.stackoverflow.com/2009/06/attribution-required/)

    """
    # base 10 exponents
    x_exp = int(floor(log10(abs(x))))
    xe_exp = int(floor(log10(abs(xe))))

    # uncertainty
    un_exp = xe_exp-precision+1
    un_int = round(xe*10**(-un_exp))

    # nominal value
    no_exp = un_exp
    no_int = round(x*10**(-no_exp))

    # format - nom(unc)exp
    fieldw = x_exp - no_exp
    fmt = '%%.%df' % fieldw
    result1 = (fmt + '(%.0f)e%d') % (no_int*10**(-fieldw), un_int, x_exp)

    # format - nom(unc)
    fieldw = max(0, -no_exp)
    fmt = '%%.%df' % fieldw
    result2 = (fmt + '(%.0f)') % (no_int*10**no_exp, un_int*10**max(0, un_exp))

    # return shortest representation
    if len(result2) <= len(result1):
        return result2
    else:
        return result1


def _number_to_X(number, uncertainty, unit, fmt, unit_fmt, fmt_pow_10, space=' '):
    uncertainty = uncertainty or getattr(number, 'uncertainty', None)
    unit = unit or unit_of(number)
    if unit is 1:
        unit_str = ''
        mag = number
    else:
        unit_str = space + unit_fmt(unit)
        mag = to_unitless(number, unit)
        if uncertainty is not None:
            uncertainty = to_unitless(uncertainty, unit)

    if uncertainty is None:
        if fmt is None:
            fmt = 5
        if isinstance(fmt, int):
            flt = ('%%.%dg' % fmt) % mag
        else:
            flt = fmt(mag)
    else:
        if fmt is None:
            fmt = 2
        if isinstance(fmt, int):
            flt = _float_str_w_uncert(mag, uncertainty, fmt)
        else:
            flt = fmt(mag, uncertainty)
    if 'e' in flt:
        significand, mantissa = flt.split('e')
        return fmt_pow_10(significand, mantissa) + unit_str
    else:
        return flt + unit_str


def _latex_pow_10(significand, mantissa):
    if significand in ('1', '1.0'):
        fmt = '10^{%s}'
    else:
        fmt = significand + r'\cdot 10^{%s}'
    return fmt % str(int(mantissa))


def number_to_scientific_latex(number, uncertainty=None, unit=None, fmt=None):
    r""" Formats a number as LaTeX (optionally with unit/uncertainty)

    Parameters
    ----------
    number : float (w or w/o unit)
    uncertainty : same as number
    unit : unit
    fmt : int or callable

    Examples
    --------
    >>> number_to_scientific_latex(3.14) == '3.14'
    True
    >>> number_to_scientific_latex(3.14159265e-7)
    '3.1416\\cdot 10^{-7}'
    >>> import quantities as pq
    >>> number_to_scientific_latex(2**0.5 * pq.m / pq.s)
    '1.4142\\,\\mathrm{\\frac{m}{s}}'
    >>> number_to_scientific_latex(1.23456, .789, fmt=2)
    '1.23(79)'

    """
    return _number_to_X(number, uncertainty, unit, fmt, latex_of_unit, _latex_pow_10, r'\,')


def _unicode_pow_10(significand, mantissa):
    if significand in ('1', '1.0'):
        result = u'10'
    else:
        result = significand + u'·10'
    return result + u''.join(map(_unicode_sup.get, str(int(mantissa))))


def number_to_scientific_unicode(number, uncertainty=None, unit=None, fmt=None):
    u""" Formats a number as unicode (optionally with unit/uncertainty)

    Parameters
    ----------
    number : float (w or w/o unit)
    uncertainty : same as number
    unit : unit
    fmt : int or callable

    Examples
    --------
    >>> number_to_scientific_unicode(3.14) == u'3.14'
    True
    >>> number_to_scientific_unicode(3.14159265e-7) == u'3.1416·10⁻⁷'
    True
    >>> import quantities as pq
    >>> number_to_scientific_unicode(2**0.5 * pq.m / pq.s)
    '1.4142 m/s'

    """
    return _number_to_X(number, uncertainty, unit, fmt, unicode_of_unit, _unicode_pow_10)


def _html_pow_10(significand, mantissa):
    if significand in ('1', '1.0'):
        result = '10<sup>'
    else:
        result = significand + '&sdot;10<sup>'
    return result + str(int(mantissa)) + '</sup>'


def number_to_scientific_html(number, uncertainty=None, unit=None, fmt=None):
    r""" Formats a number as HTML (optionally with unit/uncertainty)

    Parameters
    ----------
    number : float (w or w/o unit)
    uncertainty : same as number
    unit : unit
    fmt : int or callable

    Examples
    --------
    >>> number_to_scientific_html(3.14) == '3.14'
    True
    >>> number_to_scientific_html(3.14159265e-7)
    '3.1416&sdot;10<sup>-7</sup>'
    >>> number_to_scientific_html(1e13)
    '10<sup>13</sup>'
    >>> import quantities as pq
    >>> number_to_scientific_html(2**0.5 * pq.m / pq.s)
    '1.4142 m/s'

    """
    return _number_to_X(number, uncertainty, unit, fmt, html_of_unit, _html_pow_10)
