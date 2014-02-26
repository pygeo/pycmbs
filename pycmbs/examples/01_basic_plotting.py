"""
Basic plotting in pyCMBS
"""

from pycmbs.data import Data
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt

file_name = '../../pycmbs/examples/example_data/air.mon.mean.nc'

# core ot pyCMBS is the Data object.
# Read some data ...
air = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')

# a quick plot as well as a projection plot
f1 = map_plot(air)  # unprojected
f2 = map_plot(air, use_basemap=True)  # projected
plt.show()




