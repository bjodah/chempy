# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import warnings


try:
    import numpy as np
except ImportError:
    np = None


from ..printing import number_to_scientific_latex as ntsl


def plot_fit(x, y, beta, err=None, kw_data=None, kw_fit=None, ax=True):
    """ Plot the result of a fit

    Parameters
    ----------
    x : array_like
    y : array_like
    beta : array_like
    kw_data : dict
        Keyword arguments to ``plot`` for x, y data
    kw_fit : dict
        Keyword arguments to ``plot`` for fitted data
    ax : matplotlib.axes.Axes
        Alternatively ``True`` or ``None``
    label_cb : callable
        signature (beta, variance_beta) -> str

    """
    if ax is True:
        import matplotlib.pyplot as plt
        ax = plt.subplot(1, 1, 1)
    kw_data, kw_fit = kw_data or {}, kw_fit or {}
    if err is None:
        ax.plot(x, y, **kw_data)
    else:
        ax.errorbar(x, y, yerr=err, **kw_data)
    xlim = np.asarray([np.min(x), np.max(x)])
    ax.plot(xlim, beta[0] + beta[1]*xlim, marker='None', **kw_fit)


def least_squares(x, y, w=1, plot_cb=None):  # w == 1 => OLS, w != 1 => WLS
    """ Least-squares (w or w/o weights) fit to data series.

    Linear regression (unweighted or weighted).

    Parameters
    ----------
    x : array_like
    y : array_like
    w : array_like, optional
    plot_cb : callable or True
        Signature (x, y, beta), e.g. ``plot_fit`` (used with label when ``True``)

    Returns
    -------
    length 2 tuple : pair of parameter estimates (intercept and slope)
    2x2 array : variance-covariance matrix
    float : R-squared (goodness of fit)


    Examples
    --------
    >>> beta, vcv, R2 = least_squares([0, 1, 2], [1, 3, 5])
    >>> all(abs(beta - [1, 2]) < 1e-14), R2 == 1, (abs(vcv) < 1e-14).all()
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
    plot_errors = w is not 1
    sqrtw = np.sqrt(w)
    y = np.asarray(y)
    Y = y * sqrtw
    x = np.asarray(x)
    X = np.ones((x.size, 2))
    X[:, 1] = x
    if hasattr(sqrtw, 'ndim') and sqrtw.ndim == 1:
        sqrtw = sqrtw.reshape((sqrtw.size, 1))
    X *= sqrtw

    beta = np.linalg.lstsq(X, Y)[0]
    eps = X.dot(beta) - Y
    SSR = eps.T.dot(eps)  # sum of squared residuals
    vcv = SSR/(x.size - 2)*np.linalg.inv(X.T.dot(X))
    TSS = np.sum(np.square(Y - np.mean(Y)))  # total sum of squares
    R2 = 1 - SSR/TSS
    if plot_cb is not None:
        plot_cb(x, y, beta, w**-0.5 if plot_errors else None)
    if plot_cb is True:
        plot_fit(x, y, beta, w**-0.5 if plot_errors else None,
                 kw_fit={'label': '$y(x) = %s + %s \\cdot x$' % tuple(map(ntsl, beta))})
    return beta, vcv, float(R2)


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

    Examples
    --------
    >>> beta, cov, info = irls([1, 2, 3], [3, 2.5, 2.1], irls.abs_residuals)


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
