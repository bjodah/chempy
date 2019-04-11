# -*- coding: utf-8 -*-
"""
Contains rudimentary tools for regression: (iteratively) (weighted) least squares
and functions for plotting the fit from the regression analysis.
"""

from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None

from ..units import latex_of_unit, is_unitless, to_unitless, unit_of
from ..printing import number_to_scientific_latex


def plot_fit(x, y, beta, yerr=None, vcv_beta=None, r2=None, kw_data=None,
             kw_fit=None, fit_label_cb=None, ax=True,
             x_unit=1, y_unit=1, nsigma=1):
    """ Plot the result of a fit

    Parameters
    ----------
    x : array_like
    y : array_like
    beta : array_like
    yerr : array_like
    vcv_beta : array_like
    kw_data : dict
        Keyword arguments to ``plot`` for x, y data
    kw_fit : dict
        Keyword arguments to ``plot`` for fitted data
    fit_label_cb: callable:
        signature (beta, variance_beta, r2) -> str
    ax : matplotlib.axes.Axes
        Alternatively ``True`` or ``None``
    x_unit : unit
    y_unit : unit
    nsigma : int
        Multiplier for errorbars when plotting.

    """
    x_ul = to_unitless(x, x_unit)
    y_ul = to_unitless(y, y_unit)
    if ax is True:
        import matplotlib.pyplot as plt
        ax = plt.subplot(1, 1, 1)
    kw_data, kw_fit = kw_data or {}, kw_fit or {}
    if fit_label_cb is not None and 'label' not in kw_fit:
        kw_fit['label'] = fit_label_cb(beta, vcv_beta, r2)

    if yerr is None:
        ax.plot(x_ul, y_ul, **kw_data)
    else:
        ax.errorbar(x_ul, y_ul, yerr=to_unitless(yerr*nsigma, y_unit), **kw_data)

    xlim = [np.min(x_ul), np.max(x_ul)]
    if 'marker' not in kw_fit:
        kw_fit['marker'] = 'None'

    beta_ul = [to_unitless(elem, y_unit*x_unit**-i) for i, elem in enumerate(beta)]
    yfit_ul = [sum([b*x_elem**i for i, b in enumerate(beta_ul)]) for x_elem in xlim]

    ax.plot(xlim, yfit_ul, **kw_fit)
    if 'label' in kw_fit:
        ax.legend(loc='best')

    if is_unitless(x_unit):
        ax.set_xlabel('$x$')
    else:
        ax.set_xlabel('$x / %s$' % latex_of_unit(x_unit))

    if is_unitless(y_unit):
        ax.set_ylabel('$y$')
    else:
        ax.set_ylabel('$y / %s$' % latex_of_unit(y_unit))

    return ax


def _beta_tup(beta, x_unit, y_unit):
    return tuple(coeff*y_unit/x_unit**i for i, coeff in enumerate(beta))


def plot_least_squares_fit(x, y, beta_vcv_r2, yerr=None, plot_cb=None,
                           plot_cb_kwargs=None, x_unit=1, y_unit=1):
    """ Performs Least-squares fit and plots data and fitted line

    Parameters
    ----------
    x : array_like
    y : array_like
    beta_vcv_r2 : tuple
        Result from :func:`least_squares_fit`.
    plot_cb : callable
        When ``None``: uses :func:`plot_fit`, when callable:
        signature ``(x, y, beta, yerr=None, fit_label_cb=lambda beta, vcv, r2: 'None') -> str``.
    plot_cb_kwargs: dict, optional
        Keyword arguments passed on to ``plot_cb`` (see :func:`plot_fit` for list of
        expected kwargs). If ``plot_cb`` is ``True`` it will be populated with defaults
        (kw_data, fit_label_cb, x_unit, y_unit).

    """
    plot_cb_kwargs = plot_cb_kwargs or {}
    if plot_cb is None:
        kw_data = plot_cb_kwargs.get('kw_data', {})
        if 'marker' not in kw_data and len(x) < 40:
            kw_data['marker'] = 'd'
        if 'ls' not in kw_data and 'linestyle' not in kw_data and len(x) < 40:
            kw_data['ls'] = 'None'
        plot_cb_kwargs['kw_data'] = kw_data
        if 'fit_label_cb' not in plot_cb_kwargs:
            plot_cb_kwargs['fit_label_cb'] = lambda b, v, r2: (
                '$y(x) = %s + %s \\cdot x$' % tuple(map(number_to_scientific_latex, b))
            )
        plot_cb = plot_fit
    if 'x_unit' not in plot_cb_kwargs:
        plot_cb_kwargs['x_unit'] = x_unit
    if 'y_unit' not in plot_cb_kwargs:
        plot_cb_kwargs['y_unit'] = y_unit

    plot_cb(x, y, beta_vcv_r2[0], yerr, **plot_cb_kwargs)


def least_squares_units(x, y, w=1):
    """ Units-aware least-squares (w or w/o weights) fit to data series.

    Parameters
    ----------
    x : array_like
    y : array_like
    w : array_like, optional

    See also
    --------
    - :func:`least_squares`

    """
    x_unit, y_unit = unit_of(x), unit_of(y)
    explicit_errors = w is not 1
    if explicit_errors:
        if unit_of(w) == y_unit**-2:
            _w = to_unitless(w, y_unit**-2)
        elif unit_of(w) == unit_of(1):
            _w = w
        else:
            raise ValueError("Incompatible units in y and w")
    else:
        _w = 1
    _x = to_unitless(x, x_unit)
    _y = to_unitless(y, y_unit)
    beta, vcv, r2 = least_squares(_x, _y, _w)
    beta_tup = _beta_tup(beta, x_unit, y_unit)
    return beta_tup, vcv, float(r2)


def least_squares(x, y, w=1):  # w == 1 => OLS, w != 1 => WLS
    """ Least-squares (w or w/o weights) fit to data series.

    Linear regression (unweighted or weighted).

    Parameters
    ----------
    x : array_like
    y : array_like
    w : array_like, optional

    Returns
    -------
    length 2 tuple : pair of parameter estimates (intercept and slope)
    2x2 array : variance-covariance matrix
    float : R-squared (goodness of fit)

    Examples
    --------
    >>> import numpy as np
    >>> beta, vcv, R2 = least_squares([0, 1, 2], [1, 3, 5])
    >>> all(abs(beta - np.array([1, 2])) < 1e-14), R2 == 1, (abs(vcv) < 1e-14).all()
    (True, True, True)
    >>> b1, v1, r2_1 = least_squares([1, 2, 3], [0, 1, 4], [1, 1, 1])
    >>> b2, v2, r2_2 = least_squares([1, 2, 3], [0, 1, 4], [1, 1, .2])
    >>> abs(b2[1] - 1) < abs(b1[1] - 1)
    True

    References
    ----------
    Wikipedia & standard texts on least sqaures method.
    Comment regarding R2 in WLS:
        Willett, John B., and Judith D. Singer. "Another cautionary note about R 2:
        Its use in weighted least-squares regression analysis."
        The American Statistician 42.3 (1988): 236-238.

    """
    sqrtw = np.sqrt(w)
    Y = np.asarray(y, dtype=np.float64) * sqrtw
    _x = np.asarray(x)
    X = np.ones((_x.size, 2))
    X[:, 1] = x
    if hasattr(sqrtw, 'ndim') and sqrtw.ndim == 1:
        sqrtw = sqrtw.reshape((sqrtw.size, 1))
    X *= sqrtw

    beta = np.linalg.lstsq(X, Y, rcond=2e-16*_x.size)[0]
    eps = X.dot(beta) - Y
    SSR = eps.T.dot(eps)  # sum of squared residuals
    vcv = SSR/(_x.size - 2)*np.linalg.inv(X.T.dot(X))
    TSS = np.sum(np.square(Y - np.mean(Y)))  # total sum of squares
    R2 = 1 - SSR/TSS
    return beta, vcv, R2


def irls_units(x, y, **kwargs):
    """ Units aware version of :func:`irls`

    Parameters
    ----------
    x : array_like
    y : array_like
    \\*\\*kwargs
        Keyword arguments passed on to :func:`irls`

    See also
    --------
    - :func:`irls`

    """
    x_unit, y_unit = unit_of(x), unit_of(y)
    x_ul, y_ul = to_unitless(x, x_unit), to_unitless(y, y_unit)
    beta, vcv, info = irls(x_ul, y_ul, **kwargs)
    beta_tup = _beta_tup(beta, x_unit, y_unit)
    return beta_tup, vcv, info


def irls(x, y, w_cb=lambda x, y, b, c: x**0, itermax=16, rmsdwtol=1e-8):
    """ Iteratively reweighted least squares

    Parameters
    ----------
    x : array_like
    y : array_like
    w_cb : callbable
        Weight callback, signature ``(x, y, beta, cov) -> weight``.
        Predefined:
            - ``irls.ones``: unit weights (default)
            - ``irls.exp``: :math:`\\mathrm{e}^{-\\beta_2 x}`
            - ``irls.gaussian``: :math:`\\mathrm{e}^{-\\beta_2 x^2}`
            - ``irls.abs_residuals``: :math:`\\lvert \\beta_1 + \\beta_2 x - y \\rvert`
    itermax : int
    rmsdwtol : float
    plot_cb : callble
        See :func:`least_squares`
    plot_cb_kwargs : dict
        See :func:`least_squares`

    Returns
    -------
    beta : length-2 array
        parameters
    cov : 2x2 array
        variance-covariance matrix
    info : dict
        Contains
           - success : bool
           - niter : int
           - weights : list of weighting arrays

    # Examples
    # --------
    # beta, cov, info = irls([1, 2, 3], [3, 2.5, 2.1], irls.abs_residuals)

    """
    if itermax < 1:
        raise ValueError("itermax must be >= 1")
    weights = []
    x, y = np.asarray(x), np.asarray(y)
    w = np.ones_like(x)
    rmsdw = np.inf
    ii = 0
    while rmsdw > rmsdwtol and ii < itermax:
        weights.append(w)
        beta, cov, r2 = least_squares(x, y, w)
        old_w = w.copy()
        w = w_cb(x, y, beta, cov)
        rmsdw = np.sqrt(np.mean(np.square(w - old_w)))
        ii += 1

    return beta, cov, {'weights': weights, 'niter': ii, 'success': ii < itermax}


irls.ones = lambda x, y, b, c: 1

if np is not None:
    irls.exp = lambda x, y, b, c: np.exp(b[1]*x)
    irls.gaussian = lambda x, y, b, c: np.exp(-(b[1]*x)**2)  # guassian weighting
    irls.abs_residuals = lambda x, y, b, c: np.abs(b[0] + b[1]*x - y)


def plot_avg_params(opt_params, cov_params, avg_params_result, label_cb=None, ax=None,
                    title=False, xlabel=False, ylabel=False, flip=False, nsigma=1):
    """ Calculates the average parameters from a set of regression parameters

    Parameters
    ----------
    opt_params : array_like
        Of shape ``(nfits, nparams)``.
    cov_params : array_like
        of shape (nfits, nparams, nparams)
    avg_params_result : length-2 tuple
       Result from :func:`avg_parrams`.
    label_cb : callable
        signature (beta, variance_beta) -> str
    ax : matplotlib.axes.Axes
    title : bool or str
    xlabel : bool or str
    ylabel : bool or str
    flip : bool
        for plotting: (x, y) -> beta1, beta0
    nsigma : int
        Multiplier for error bars

    Returns
    -------
    avg_beta: weighted average of parameters
    var_avg_beta: variance-covariance matrix

    """
    avg_beta, var_avg_beta = avg_params_result
    import matplotlib.pyplot as plt
    if label_cb is not None:
        lbl = label_cb(avg_beta, var_avg_beta)
    else:
        lbl = None
    if ax is None:
        ax = plt.subplot(1, 1, 1)
    xidx, yidx = (1, 0) if flip else (0, 1)
    opt_params = np.asarray(opt_params)
    cov_params = np.asarray(cov_params)
    var_beta = np.vstack((cov_params[:, 0, 0], cov_params[:, 1, 1])).T
    ax.errorbar(opt_params[:, xidx], opt_params[:, yidx], marker='s', ls='None',
                xerr=nsigma*var_beta[:, xidx]**0.5,
                yerr=nsigma*var_beta[:, yidx]**0.5)
    if xlabel:
        if xlabel is True:
            xlabel = r'$\beta_%d$' % xidx
        ax.set_xlabel(xlabel)
    if ylabel:
        if ylabel is True:
            xlabel = r'$\beta_%d$' % yidx
        ax.set_ylabel(ylabel)
    if title:
        if title is True:
            title = r'$y(x) = \beta_0 + \beta_1 \cdot x$'
        ax.set_title(title)
    ax.errorbar(avg_beta[xidx],
                avg_beta[yidx],
                xerr=nsigma*var_avg_beta[xidx]**0.5,
                yerr=nsigma*var_avg_beta[yidx]**0.5, marker='o', c='r',
                linewidth=2, markersize=10, label=lbl)
    ax.legend(numpoints=1)


def avg_params(opt_params, cov_params):
    """ Calculates the average parameters from a set of regression parameters.

    Parameters
    ----------
    opt_params : array_like
        of shape (nfits, nparams)
    cov_params : array_like
        of shape (nfits, nparams, nparams)

    Returns
    -------
    avg_beta: weighted average of parameters
    var_avg_beta: variance-covariance matrix

    """
    opt_params = np.asarray(opt_params)
    cov_params = np.asarray(cov_params)
    var_beta = np.vstack((cov_params[:, 0, 0], cov_params[:, 1, 1])).T
    avg_beta, sum_of_weights = np.average(opt_params, axis=0, weights=1/var_beta, returned=True)
    var_avg_beta = np.sum(np.square(opt_params - avg_beta)/var_beta, axis=0)/((avg_beta.shape[0] - 1) * sum_of_weights)
    return avg_beta, var_avg_beta
