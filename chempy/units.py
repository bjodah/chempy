try:
    import quantities as pq
except ImportError:
    default_units = None
else:
    default_units = pq


def allclose(a, b, rtol=1e-8, atol=None):
    d = abs(a - b)
    lim = abs(a)*rtol
    if atol is not None:
        lim += atol
    try:
        n = len(d)
    except TypeError:
        n = 1

    if n == 1:
        return d < lim
    else:
        return all(_d < _lim for _d, _lim in zip(d, lim))
