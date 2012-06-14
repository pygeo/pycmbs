#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

'''
statistical module for pyCMBS
'''


from scipy import stats
import numpy as np

def get_significance(correlation,n,pthres=1.01):
    """
    calculate significance of correlation

    @param correlation: pearson correlation coefficient
    @type  correlation: float

    @param n: number of samples
    @type  n: integer

    @param pthres: (optional) specifies threshold for p-value. Everything
        above this threshold will be masked
    @type  n: integer

    @return: returns the p-value
    @rtype : numpy masked array

    """

    nf = n - 2. #degree of freedom
    t_value = np.abs(correlation) * np.sqrt( nf / (1.-correlation**2)) # abs() is important

    # calculate two-sided p-value
    p =   2. * stats.t.sf(t_value,nf)

    return np.ma.array(p,mask=p>pthres)
