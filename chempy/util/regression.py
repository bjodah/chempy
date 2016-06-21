# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import warnings


try:
    import numpy as np
except ImportError:
    np = None

from ..units import latex_of_unit, is_unitless, to_unitless, unit_of
from ..printing import number_to_scientific_latex


def plot_fit(x, y, beta, yerr=None, vcv_beta=None, r2=None, kw_data=None, kw_fit=None, fit_label_cb=None, ax=True,
             x_unit=None, y_unit=None):
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
        ax.errorbar(x, y, yerr=yerr, **kw_data)
    xlim = [np.min(x_ul), np.max(x_ul)]
    if 'marker' not in kw_fit:
        kw_fit['marker'] = 'None'
    yfit_ul = to_unitless(beta[0] + beta[1]*xlim*x_unit, y_unit)
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


def least_squares(x, y, w=1, plot_cb=None, plot_cb_kwargs=None):  # w == 1 => OLS, w != 1 => WLS
    """ Least-squares (w or w/o weights) fit to data series.

    Linear regression (unweighted or weighted).

    Parameters
    ----------
    x : array_like
    y : array_like
    w : array_like, optional
    plot_cb : callable or True
        When ``True``: uses :func:`plot_fit`, when callable:
        signature ``(x, y, beta, yerr=None, fit_label_cb=lambda beta, vcv, r2: 'None') -> str``,

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
    x_unit, y_unit = unit_of(x), unit_of(y)
    explicit_errors = w is not 1
    if explicit_errors:
        if unit_of(w) != y_unit**-2:
            raise ValueError("Incompatible units in y and w")
        _w = to_unitless(w, y_unit**-2)
    else:
        _w = 1
    sqrtw = np.sqrt(_w)
    _y = to_unitless(y, y_unit)
    Y = _y * sqrtw
    _x = to_unitless(x, x_unit)
    X = np.ones((_x.size, 2))
    X[:, 1] = _x
    if hasattr(sqrtw, 'ndim') and sqrtw.ndim == 1:
        sqrtw = sqrtw.reshape((sqrtw.size, 1))
    X *= sqrtw

    beta = np.linalg.lstsq(X, Y)[0]
    eps = X.dot(beta) - Y
    SSR = eps.T.dot(eps)  # sum of squared residuals
    vcv = SSR/(_x.size - 2)*np.linalg.inv(X.T.dot(X))
    TSS = np.sum(np.square(Y - np.mean(Y)))  # total sum of squares
    R2 = 1 - SSR/TSS
    plot_cb_kwargs = plot_cb_kwargs or {}
    if plot_cb is True:
        # plot_cb == True give some convenient defaults
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

    beta_tup = tuple(coeff*y_unit*x_unit**-i for i, coeff in enumerate(beta))
    if plot_cb is not None:
        plot_cb(x, y, beta_tup, w**-0.5 if explicit_errors else None, **plot_cb_kwargs)

    return beta_tup, vcv, float(R2)


def irls(x, y, w_cb=lambda x, y, b, c: x**0, itermax=16, rmsdwtol=1e-8, full_output=False, ax=None):
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
    full_output : bool
    ax : matplotlib.axes.Axes

    Returns
    -------
    beta : length-2 array
        parameters
    cov : 2x2 array
        variance-covariance matrix
    info : dict
        if ``full_output == True``, keys:
           - weights
           - niter

    # Examples
    # --------
    # beta, cov, info = irls([1, 2, 3], [3, 2.5, 2.1], irls.abs_residuals)


    """

    if itermax < 1:
        raise ValueError("itermax must be >= 1")
    weights = []
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
    if ii == itermax:
        warnings.warn('itermax reached')
    if ax is not None:
        plot_fit(x, y, beta, ax=ax)
    if full_output:
        return beta, cov, {'weights': weights, 'niter': ii, 'success': ii < itermax}
    else:
        return beta, cov

irls.ones = lambda x, y, b, c: 1

if np is not None:
    irls.exp = lambda x, y, b, c: np.exp(b[1]*x)
    irls.gaussian = lambda x, y, b, c: np.exp(-(b[1]*x)**2)  # guassian weighting
    irls.abs_residuals = lambda x, y, b, c: np.abs(b[0] + b[1]*x - y)


def avg_params(opt_params, cov_params, label_cb=None, ax=None,
               title=r'$y(x) = \beta_0 + \beta_1 \cdot x$',
               xlabel=r'$\beta_0$', ylabel=r'$\beta_1$'):
    """ Calculates the average parameters from a set of regression parameters

    Parameters
    ----------
    opt_params : array_like
        of shape (nfits, nparams)
    cov_params : array_like
        of shape (nfits, nparams, nparams)
    label_cb : callable
        signature (beta, variance_beta) -> str
    ax : matplotlib.axes.Axes
    xlabel : str
    ylabel : str

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
    if ax is not None:
        import matplotlib.pyplot as plt
        if label_cb is not None:
            lbl = label_cb(avg_beta, var_avg_beta)
        else:
            lbl = None
        if ax is True:
            ax = plt.subplot(1, 1, 1)
        ax.errorbar(opt_params[:, 0], opt_params[:, 1], marker='s', ls='None',
                    xerr=var_beta[:, 0]**0.5, yerr=var_beta[:, 1]**0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.errorbar(avg_beta[0], avg_beta[1], xerr=var_avg_beta[0]**0.5, yerr=var_avg_beta[1]**0.5, marker='o', c='r',
                    linewidth=2, markersize=10, label=lbl)
        ax.legend(numpoints=1)
    return avg_beta, var_avg_beta
