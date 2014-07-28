# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from scipy import stats
import numpy as np
import scipy.special as special
from scipy import stats as Sstats
import sys


def get_significance(correlation, n, pthres=1.01):
    """
    calculate significance of correlation

    Parameters
    ----------
    correlation : float
        pearson correlation coefficient
    n : int
        number of samples
    pthres : float
        (optional) specifies threshold for p-value. Everything
        above this threshold will be masked

    Returns
    -------
    ndarray : returns the p-value
    """

    nf = n - 2  # degree of freedom

    if np.isscalar(correlation):
        if abs(correlation) < 1.:
            # abs() is important
            t_value = np.abs(correlation) * np.sqrt(nf / (1. - correlation ** 2.))
            # calculate two-sided p-value
            p = 2. * stats.t.sf(t_value, nf)
        else:
            p = 0.
    else:
        msk = correlation < 1.
        p = np.zeros_like(correlation)
        t_value = np.abs(correlation[msk]) * np.sqrt(nf / (1. - correlation[msk] ** 2.))
        p[msk] = 2. * stats.t.sf(t_value, nf)

    return np.ma.array(p, mask=p > pthres)

#-------------------------------------------------------------------------
# The following routines are derived from scipy.mstats.mstats_basic.py
#
# the reason to include them here is that there was a bug in the code
# which I reported to scipy developers. Once this is fixed,
# we can remove it again here
#
# http://projects.scipy.org/scipy/ticket/1777


def _chk2_asarray(a, b, axis):
    if axis is None:
        a = np.ma.ravel(a)
        b = np.ma.ravel(b)
        outaxis = 0
    else:
        a = np.ma.asanyarray(a)
        b = np.ma.asanyarray(b)
        outaxis = axis
    return a, b, outaxis


def betai(a, b, x):
    x = np.asanyarray(x)
    x = np.ma.where(x < 1.0, x, 1.0)  # if x > 1 then return 1.0
    return special.betainc(a, b, x)


def ttest_ind(a, b, axis=0):
    a, b, axis = _chk2_asarray(a, b, axis)
    (x1, x2) = (a.mean(axis), b.mean(axis))
    (v1, v2) = (a.var(axis=axis, ddof=1), b.var(axis=axis, ddof=1))
    (n1, n2) = (a.count(axis), b.count(axis))
    df = n1 + n2 - 2
    #AL <<<<<<<<<<<<<< fix, as float() functions from mstats_basic.py
    #does not work for multidimensional arrays!
    svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / (df * 1.)
    #svar == 0
    # N-D COMPUTATION HERE!!!!!!
    t = (x1 - x2) / np.ma.sqrt(svar * (1.0 / n1 + 1.0 / n2))
    t = np.ma.filled(t, 1)           # replace NaN t-values with 1.0
    #AL <<<<<<<<<<<<<<
    probs = betai(0.5 * df, 0.5, (df * 1.) / (df + t * t)).reshape(t.shape)
    # .squeeze() #<<< AL removed the squeeze, so I get back an array!
    return t, probs


#----------------------- END OF SCIPY IMPORT

def welchs_approximate_ttest(n1, mean1, sem1, n2, mean2, sem2, alpha):
    """
    REFERENCES:
    http://econpy.googlecode.com/svn/trunk/pytrix/stat.py
    http://comments.gmane.org/gmane.comp.python.scientific.user/12907

    Welch''s approximate t-test for the difference of two means of
    heteroscedasctic populations.

    :see: Biometry, Sokal and Rohlf, 3rd ed., 1995, Box 13.4

    :Parameters:
        n1 : int
            number of variates in sample 1
        n2 : int
            number of variates in sample 2
        mean1 : float
            mean of sample 1
        mean2 : float
            mean of sample 2
        sem1 : float
            standard error of mean1
        sem2 : float
            standard error of mean2
        alpha : float
            desired level of significance of test

    :Returns:
        significant : bool
            True if means are significantly different, else False
        t_s_prime : float
            t_prime value for difference of means
        t_alpha_prime : float
            critical value of t_prime at given level of significance

    :author: Angus McMorland
    :license: BSD_

    .. BSD: http://www.opensource.org/licenses/bsd-license.php
    """
    svm1 = sem1 ** 2 * n1
    svm2 = sem2 ** 2 * n2
    t_s_prime = (mean1 - mean2) / np.sqrt(svm1 / n1 + svm2 / n2)

    t_alpha_df1 = Sstats.t.ppf(1 - alpha / 2, n1 - 1)
    t_alpha_df2 = Sstats.t.ppf(1 - alpha / 2, n2 - 1)
    t_alpha_prime = (t_alpha_df1 * sem1 ** 2 + t_alpha_df2 * sem2 ** 2) /\
                    (sem1 ** 2 + sem2 ** 2)
    return abs(t_s_prime) > t_alpha_prime, t_s_prime, t_alpha_prime




def lomb_scargle_periodogram(t, p, y):
    """
    calculate the Lomb-Scargle periodogram
    This corresponds to a method to perform spectral analyis
    based on unevenly sampled data and/or data with gaps

    The estimation is based on a linear regression of a cosine model
    for each of the input frequencies

    References
    ----------
    [1] Hocke & Kaempfer: http://www.atmos-chem-phys.net/9/4197/2009/

    Parameters
    ----------
    t : ndarray
        time array [e.g. units of days]

    p : ndarray
        array with desired periods [day]; in general any kind of timeunit is
        valid, but needs to be simply consistent with time
    y : ndarray
        observations to fit

    Returns
    -------
    returns
    a) Amplitude for each frequency
    b) phase for each frequency
    c) Periodogram for each frequency

    Example
    -------

    """
    from scipy import optimize

    def func(x, A, B):
        return A*np.cos(x + B)

    resA = np.ones(len(p))
    resB = np.ones(len(p))

    cnt = 0
    for period in p:
        f = 2.*np.pi * t / period
        popt, pcov = optimize.curve_fit(func, f, y, p0=[1., 0.])

        resA[cnt] = popt[0]
        resB[cnt] = popt[1]
        del popt, pcov
        cnt += 1
    return resA, resB










