# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

try:
    import numpy as np
except ImportError:
    np = None

import warnings


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
    sqrtw = np.sqrt(w)
    y = np.asarray(y) * sqrtw
    x = np.asarray(x)
    X = np.ones((x.size, 2))
    X[:, 1] = x
    if hasattr(sqrtw, 'ndim') and sqrtw.ndim == 1:
        sqrtw = sqrtw.reshape((sqrtw.size, 1))
    X *= sqrtw

    beta = np.linalg.lstsq(X, y)[0]
    eps = X.dot(beta) - y
    SSR = eps.T.dot(eps)  # sum of squared residuals
    vcv = SSR/(x.size - 2)*np.linalg.inv(X.T.dot(X))
    TSS = np.sum(np.square(y - np.mean(y)))  # total sum of squares
    R2 = 1 - SSR/TSS
    return beta, vcv, float(R2)


def irls(x, y, w_cb=lambda x, y, b, c: x**0, itermax=16, rmsdwtol=1e-8, full_output=False):
    """ Iteratively reweighted least squares

    Parameters
    ----------
    x : array_like
    y : array_like
    w_cb : callbable
        signature (x, y, beta, cov) -> weight
    itermax : int
    rmsdwtol : float
    full_output : bool

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

    """

    if itermax < 1:
        raise ValueError("intermax must be >= 1")
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
    if full_output:
        return beta, cov, {'weights': weights, 'niter': ii, 'success': ii < itermax}
    else:
        return beta, cov

irls.ones = lambda x, y, b, c: 1

if np is not None:
    irls.exp = lambda x, y, b, c: np.exp(b[1]*x)
    irls.gaussian = lambda x, y, b, c: np.exp(-(b[1]*x)**2)  # guassian weighting
    irls.abs_residuals = lambda x, y, b, c: np.abs(np.exp(b[0] + b[1]*x) - np.exp(y))


def plot_fit(x, y, beta, kw_data=None, kw_fit=None):
    import matplotlib.pyplot as plt
    kw_data, kw_fit = kw_data or {}, kw_fit or {}
    plt.plot(x, y, linestyle='None', **kw_data)
    plt.plot(x[[0, -1]], beta[0] + beta[1]*x[[0, -1]], marker='None', **kw_fit)


def avg_params(opt_params, cov_params, label_cb=None, plot=False):
    var_beta = np.vstack((cov_params[:, 0, 0], cov_params[:, 1, 1])).T
    avg_beta, sum_of_weights = np.average(opt_params, axis=0, weights=1/var_beta, returned=True)
    var_avg_beta = np.sum(np.square(opt_params - avg_beta)/var_beta, axis=0)/((avg_beta.shape[0] - 1) * sum_of_weights)
    if plot:
        import matplotlib.pyplot as plt
        if label_cb is not None:
            lbl = label_cb(avg_beta, var_avg_beta)
        else:
            lbl = None
        plt.errorbar(opt_params[:, 0], opt_params[:, 1], marker='s', ls='None',
                     xerr=var_beta[:, 0]**0.5, yerr=var_beta[:, 1]**0.5)
        plt.xlabel(r'$\beta_0$')
        plt.ylabel(r'$\beta_1$')
        plt.title(r'$y(x) = \beta_0 + \beta_1 \cdot x$')
        plt.errorbar(avg_beta[0], avg_beta[1], xerr=var_avg_beta[0]**0.5, yerr=var_avg_beta[1]**0.5, marker='o', c='r',
                     linewidth=2, markersize=10, label=lbl)
        plt.legend(numpoints=1)
    return avg_beta, var_avg_beta
