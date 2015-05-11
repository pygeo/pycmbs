# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import numpy as np


def lomb_scargle_periodogram(t, p, y, corr=True):
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
    corr : bool
        calculate also correlation of model with data (quality of fit)

    Returns
    -------
    returns
    a) A: Amplitude for each frequency
    b) B: phase for each frequency

    optional
    c) R: peasson correlation coefficient of model with data for each frequency
    d) P: p-value for linear correlation

    Example
    -------
    see file lomb.py in scripts subdirectory

    """
    from scipy import optimize
    from scipy import stats

    def func(x, A, B):
        return A * np.cos(x + B)

    resA = np.ones(len(p))
    resB = np.ones(len(p))
    resR = np.ones(len(p))
    resP = np.ones(len(p))

    cnt = 0
    hlp = 2. * np.pi * t

    for period in p:
        popt, pcov = optimize.curve_fit(func, hlp / period, y, p0=[1., 0.])   # f=2*pi*t/period
        resA[cnt] = popt[0]
        resB[cnt] = popt[1]
        if corr:
            ymod = func(hlp / period, resA[cnt], resB[cnt])
            ymod = np.ma.array(ymod, mask=ymod != ymod)

            import pickle
            pickle.dump({'y': y, 'ymod': ymod}, open('y.pkl', 'w'))

            #~ print ''

            slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(ymod, y)
            #~ print cnt, r_value, p_value

            resR[cnt] = r_value
            resP[cnt] = p_value
#~
            #~ print type(ymod), type(y)
            #~ print np.corrcoef(ymod, y)
            #~ print 'A,B: ', resA[cnt], resB[cnt], popt

        del popt, pcov
        cnt += 1

    if corr:
        return resA, resB, resR, resP
    else:
        return resA, resB
