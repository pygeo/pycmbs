#!/usr/bin/python
# -*- coding: utf-8 -*-

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






