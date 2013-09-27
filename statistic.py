#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
'''



from scipy import stats
import numpy as np
import scipy.special as special
from scipy import stats as Sstats

def get_significance(correlation,n,pthres=1.01):
    """
    calculate significance of correlation

    @param correlation: pearson correlation coefficient
    @type  correlation: float

    @param n: number of samples
    @type  n: integer

    @param pthres: (optional) specifies threshold for p-value. Everything
        above this threshold will be masked
    @type  pthres: integer

    @return: returns the p-value
    @rtype : numpy masked array

    """

    nf = n - 2. #degree of freedom
    t_value = np.abs(correlation) * np.sqrt( nf / (1.-correlation**2)) # abs() is important

    # calculate two-sided p-value
    p =   2. * stats.t.sf(t_value,nf)

    return np.ma.array(p,mask=p>pthres)

#-------------------------------------------------------------------------
# The following routines are derived from scipy.mstats.mstats_basic.py
#
# the reason to include them here is that there was a bug in the code
# which I reported to scipy developers. Once this is fixed, we can remove it again here
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
    df = n1+n2-2
    svar = ((n1-1)*v1+(n2-1)*v2) / (df*1.)  #AL <<<<<<<<<<<<<< fix, as float() functions from mstats_basic.py does not work for multidimensional arrays!
    #svar == 0
    t = (x1-x2)/np.ma.sqrt(svar*(1.0/n1 + 1.0/n2))  # N-D COMPUTATION HERE!!!!!!
    t = np.ma.filled(t, 1)           # replace NaN t-values with 1.0
    probs = betai(0.5*df,0.5,(df*1.)/(df+t*t)).reshape(t.shape)   #AL <<<<<<<<<<<<<<
    return t, probs #.squeeze() #<<< AL removed the squeeze, so I get back an array!


#----------------------- END OF SCIPY IMPORT




def welchs_approximate_ttest(n1, mean1, sem1, n2, mean2, sem2, alpha):
    '''
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
    '''
    svm1 = sem1**2 * n1
    svm2 = sem2**2 * n2
    t_s_prime = (mean1 - mean2)/np.sqrt(svm1/n1+svm2/n2)

    t_alpha_df1 = Sstats.t.ppf(1-alpha/2, n1 - 1)
    t_alpha_df2 = Sstats.t.ppf(1-alpha/2, n2 - 1)
    t_alpha_prime = (t_alpha_df1 * sem1**2 + t_alpha_df2 * sem2**2) /\
                    (sem1**2 + sem2**2)
    return abs(t_s_prime) > t_alpha_prime, t_s_prime, t_alpha_prime