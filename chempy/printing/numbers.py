# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


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


def number_to_scientific_latex(number, fmt='%.3g'):
    r"""
    Examples
    --------
    >>> number_to_scientific_latex(3.14) == '3.14'
    True
    >>> number_to_scientific_latex(3.14159265e-7)
    '3.14\\cdot 10^{-7}'
    >>> import quantities as pq
    >>> number_to_scientific_latex(2**0.5 * pq.m / pq.s)
    '1.41 \\mathrm{\\frac{m}{s}}'

    """
    try:
        unit = ' ' + number.dimensionality.latex.strip('$')
        number = number.magnitude
    except AttributeError:
        unit = ''
    s = fmt % number
    if 'e' in s:
        prefix, suffix = s.split('e')
        if prefix in ('1', '1.0'):
            result = '10^{%s}'
        else:
            result = prefix + r'\cdot 10^{%s}'
        return result % str(int(suffix)) + unit
    else:
        return s + unit


def number_to_scientific_unicode(number, fmt='%.3g'):
    u"""
    Examples
    --------
    >>> number_to_scientific_unicode(3.14) == u'3.14'
    True
    >>> number_to_scientific_unicode(3.14159265e-7) == u'3.14·10⁻⁷'
    True
    >>> import quantities as pq
    >>> number_to_scientific_unicode(2**0.5 * pq.m / pq.s)
    '1.41 m/s'

    """
    try:
        unit = ' ' + number.dimensionality.unicode
        number = number.magnitude
    except AttributeError:
        unit = ''
    s = fmt % number
    if 'e' in s:
        prefix, suffix = s.split('e')
        if prefix in ('1', '1.0'):
            result = u'10'
        else:
            result = prefix + u'·10'
        return result + u''.join(map(_unicode_sup.get, str(int(suffix)))) + unit
    else:
        return s + unit


def number_to_scientific_html(number, fmt='%.3g'):
    """
    Examples
    --------
    >>> number_to_scientific_html(3.14) == '3.14'
    True
    >>> number_to_scientific_html(3.14159265e-7)
    '3.14&sdot;10<sup>-7</sup>'
    >>> number_to_scientific_html(1e13)
    '10<sup>13</sup>'
    >>> import quantities as pq
    >>> number_to_scientific_html(2**0.5 * pq.m / pq.s)
    '1.41 m/s'

    """
    try:
        unit = ' ' + str(number.dimensionality)
        number = number.magnitude
    except AttributeError:
        unit = ''
    s = fmt % number
    if 'e' in s:
        prefix, suffix = s.split('e')
        if prefix in ('1', '1.0'):
            result = '10<sup>'
        else:
            result = prefix + '&sdot;10<sup>'
        return result + str(int(suffix)) + '</sup>' + unit
    else:
        return s + unit
