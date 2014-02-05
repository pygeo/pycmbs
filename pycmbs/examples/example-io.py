# implicit imports

#from pyCMBS import *
import os
import numpy as np
from matplotlib import pylab as pl

pl.close('all')

file ='../example_data/air.mon.mean.nc'
#--- read data ---
D = Data(file, 'air', lat_name='lat',lon_name='lon',read=True)
D.save('nix.nc',delete=True)



