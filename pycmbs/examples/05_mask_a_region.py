# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from pycmbs.region import RegionBboxLatLon
from pycmbs.examples import download
from pycmbs.data import Data
from pycmbs.mapping import map_plot

import matplotlib.pyplot as plt

# specify some region using bounding box
# here: 20deg W ... 30 DEG E, 40 DEG S ... 5 DEG N
r = RegionBboxLatLon(777, -20.,30.,-40.,5., label='testregion')   #777 is just the ID value

# read some data as Data object
filename = download.get_sample_file(name='air', return_object=False)
air = Data(filename, 'air', read=True)

# generate some plot BEFORE the masking
map_plot(air, title='before', use_basemap=True)

# now mask the data ...
air.get_aoi_lat_lon(r)

# generate some plot AFTER the masking
map_plot(air, title='after', use_basemap=True)

#... o.k., so far so good, but the dataset "air" still contains data for the entire domain.
# even if it is masked it will eat some of your memory. You can see this by plotting the size of the matrix
print (air.shape)

# wouldn't it be nice to just remove everything which is not needed?
# Here you go ...
air.cut_bounding_box()
map_plot(air, title='after cutting', use_basemap=True)

# still looks the same, isn't  it? No it isn't, look at the size!
print (air.shape)

# To illustrate what happens, you can plot WITHOUT using a map projection
map_plot(air, title='after cutting (no projection)')

plt.show()

