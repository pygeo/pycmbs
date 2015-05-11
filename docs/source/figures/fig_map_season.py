# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""
from pycmbs.data import Data
from pycmbs.plots import map_season
import matplotlib.pyplot as plt

file_name = '../../../pycmbs/examples/example_data/air.mon.mean.nc'
air = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')
c = air.get_climatology(return_object=True)

# a quick plot as well as a projection plot
f1 = map_season(c, show_stat=False, vmin=-30., vmax=30., cticks=[-30., 0., 30.])  # unprojected
plt.show()
