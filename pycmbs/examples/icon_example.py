# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from pycmbs.icon import Icon
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt

plt.close('all')

##### PLEASE NOTE THAT THIS EXAMPLE IS CURRENTLY LIMITED USERS WITH ACCESS TO SPECIFIC DATA


#//// read data ////
gridfile ='../..//example_data/icon/r2b4_amip.nc'
datafile = '../../example_data/icon/rms0006_atm_phy_DOM01_ML_0001.nc'

#/// read data ///
IC = Icon(datafile,gridfile,'rsns')
IC.read()
IC.label='shortwave net flux at surface'
IC.unit='$W/^2$'

#/// Do plotting ///
print('Doing first plot ...')
map_plot(IC, use_basemap=True)
print('Doing second plot ...')
map_plot(IC, use_basemap=False)

plt.show()
