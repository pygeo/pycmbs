"""
Test map layouts
"""
from pycmbs.mapping import SingleMap
from pycmbs.data import Data
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as grd
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.close('all')


file='testdata.nc'
d=Data(file,'tmp',read=True)


import cartopy.crs as ccrs

fig = plt.figure()
ax1= fig.add_subplot(2,1,1)
ax2= fig.add_subplot(2,1,2, projection=ccrs.Robinson())



# map only with colorbars
m = SingleMap(d, ax=ax1, backend='basemap')  # this is supposed to make a baemap plot with stripes
m.plot(colorbar_orientation='vertical', vmin=10., vmax=20., proj_prop={'projection':'robin', 'lon_0':0.})






#~ box_top = 45
#~ x, y = [-44, -44, 45, 45, -44], [-45, box_top, box_top, -45, -45]

x, y = [-110., 0.], [0., 50.]

#~ rotated_pole = ccrs.RotatedPole(pole_latitude=45, pole_longitude=180)

#~ rotated_pole = ccrs.Robinson(central_longitude=0.)

#~ ax = plt.subplot(111, projection=ccrs.Robinson())
ax2.set_global()
#~ ax.stock_img()
ax2.coastlines()


ax2.pcolormesh(d.lon, d.lat, d.data[0,:,:], transform=ccrs.PlateCarree())

ax2.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
ax2.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.PlateCarree())
ax2.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic(), color='blue')  # draws the orthodrome

#~ ax.plot(x, y, marker='o', transform=rotated_pole)
#~ ax.fill(x, y, color='coral', transform=rotated_pole, alpha=0.4)
ax2.gridlines()


plt.show()
