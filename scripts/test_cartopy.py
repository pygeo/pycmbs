"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

"""
Test map layouts
"""
from pycmbs.mapping import SingleMap
from pycmbs.data import Data
import matplotlib.pyplot as plt

plt.close('all')


file='testdata.nc'
d=Data(file,'tmp',read=True)


import cartopy.crs as ccrs

fig = plt.figure(figsize=(8,10))
ax1= fig.add_subplot(2,1,1)
ax2= fig.add_subplot(2,1,2)

# map only with colorbars
m = SingleMap(d, ax=ax1, backend='basemap')  # this is supposed to make a baemap plot with stripes
m.plot(colorbar_orientation='vertical', vmin=0., vmax=30., proj_prop={'projection':'robin', 'lon_0':0.}, title='with Basemap')

m1 = SingleMap(d, ax=ax2, backend='cartopy')  # this is supposed to make a baemap plot with stripes
m1.plot(colorbar_orientation='vertical', vmin=0., vmax=30., proj_prop={'projection':'robin', 'lon_0':0.}, title='with Cartopy')

plt.show()



