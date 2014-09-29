"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""


"""
test script to test plotting of polygons on top of maps
"""
from pycmbs.mapping import SingleMap
from pycmbs.data import Data
import matplotlib.pyplot as plt
from pycmbs.polygon import Polygon
import numpy as np

plt.close('all')

poly1 = [(-105., 60.), (-168.022, 60.), (-168.022, 72.554), (-105., 72.554)]  # ALA 1
P1 = Polygon(1, poly1)

poly2 = [(-66.377, -20.), (-79.729, -1.239), (-68.8, 11.439), (-50., 11.439), (-50., -20.)]
P2 = Polygon(7, poly2)



#~ AMZ 7 (20.000S, 66.377W) (1.239S, 79.729W) (11.439N, 68.800W) (11.439N, 50.000W) (20.000S, 50.000W)
#~ CAM 6 (11.439N, 68.800W) (1.239S, 79.729W) (28.566N, 118.323W) (28.566N, 90.315W)


tmp = np.ones((180, 360))
d = Data(None, None)
d.data = np.ma.array(tmp, mask=tmp!=tmp)
d.cell_area = np.ones_like(tmp)


lon = np.arange(-180., 180.) + 0.5
lat = np.arange(-90., 90.) + 0.5
d.lon, d.lat = np.meshgrid(lon, lat)

# Basemap plots
m = SingleMap(d)  # this is supposed to make a baemap plot with stripes
m.backend = 'basemap'  # overwrite default
m._draw = m._draw_basemap
m.plot(polygons=[P1, P2], proj_prop={'projection':'robin', 'lon_0':0.})
plt.title('Basemap')

# cartopy plots
m1 = SingleMap(d, backend='cartopy')
m1.plot(polygons=[P1, P2], proj_prop={'projection':'robin', 'lon_0':0.})
plt.title('Cartopy')

plt.show()

