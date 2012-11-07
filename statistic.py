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
    @type  pthres: integer

    @return: returns the p-value
    @rtype : numpy masked array

    """

    nf = n - 2. #degree of freedom
    t_value = np.abs(correlation) * np.sqrt( nf / (1.-correlation**2)) # abs() is important

    # calculate two-sided p-value
    p =   2. * stats.t.sf(t_value,nf)

    return np.ma.array(p,mask=p>pthres)
