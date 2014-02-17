"""
Basic plotting in pyCMBS
"""

from pycmbs.data import Data
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt

file_name = '../../pycmbs/examples/example_data/air.mon.mean.nc'

air = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')

map_plot(air, show_timeseries=True, use_basemap=True, title='show_timeseries=True')
map_plot(air, show_zonal=True, use_basemap=True, title='show_zonal=True')
map_plot(air, show_histogram=True, use_basemap=True, title='show_histogram=True')

plt.show()







