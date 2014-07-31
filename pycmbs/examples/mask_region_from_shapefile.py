# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
This example shows how to mask a region which is specified from a shapefile
"""

from pycmbs.region import RegionShape
from pycmbs.examples import download
from pycmbs.mapping import map_plot
import tempfile

shp_file = './shp_data/tp_1'  # specify name of shapefile; note that it should be done WITHOUT the file extension

# read regions from shapefile
# This gives an object which contains all regions stored in the shapefile
RS = RegionShape(shp_file)

# just print the region keys for illustration
for k in RS.regions.keys():
    print k

# if you now want to generate a particular mask we can do that
# in the following example we mask the airt temperature for the
# Tibetean plateau

# load the air temperature data
air = download.get_sample_file(name='air')
air.label = 'air temperature'
org = air.copy()

# and then mask it
r_tibet = RS.regions[1]  # gives a Region object
map_plot(air, title='before')
air1 = air.copy()

# mask with region
air.mask_region(r_tibet)
map_plot(air, title='after')
#~ b = air1.mask_region(r_tibet, return_object=True) # this is an alternative which will keep the original data unchanged

# probably you want to save you rmask in a file. If you do this it is much more efficient if you want to apply it again.
file_to_store_mask = tempfile.mktemp(suffix='.nc')  # use temporary file here for this example

x = org.copy()
y = x.mask_region(r_tibet, maskfile=file_to_store_mask, return_object=True)  #this will store the mask in a netcdf file --> have a look at it!

#if you now do the same masking again, then it is much faster as it looks first if the maskfile is already there. If so its beeing used
y1 = x.mask_region(r_tibet, maskfile=file_to_store_mask, return_object=True)  #this will store the mask in a netcdf file --> have a look at it!
map_plot(air, title='quick way')

print 'Mask was stored in file: ', file_to_store_mask





