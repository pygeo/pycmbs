"""
This example shows how to mask a region which is specified from a shapefile
"""

from pycmbs.region import RegionShape
from pycmbs.examples import download
from pycmbs.mapping import map_plot

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

# and then mask it
r_tibet = RS.regions[1]  # gives a Region object
map_plot(air, title='before')
air1 = air.copy()

# mask with region
air.mask_region(r_tibet)
map_plot(air, title='after')
#~ b = air1.mask_region(r_tibet, return_object=True) # this is an alternative which will keep the original data unchanged



