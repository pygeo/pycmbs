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
routines for testing functionality in pyCMBS
'''

from pylab import *

import numpy as np

from pyCMBS import *

from scipy import stats




#/// test correlation ///


close('all')

#/// read some data ///
fname='/media/Data/dev/svn/alex/sahel_networks/data/ndvi/NDVI_conv_03_1984_2005_mmax_JAS_detrend.nc'
D1=Data(fname,'ndvi',read=True)

# generate some random vector
z=rand(21)

# correlate data matrix with this vector
R =  D1.corr_single(z)

# plot
imshow(R); colorbar()

show()






