from pycmbs.data import Data
from pycmbs.plots import map_difference
import matplotlib.pyplot as plt

file_name = '../../../pycmbs/examples/example_data/air.mon.mean.nc'
A = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')
B = A.copy()
B.mulc(2.3, copy=False)
a = A.get_climatology(return_object=True)
b = B.get_climatology(return_object=True)

# a quick plot as well as a projection plot
f1 = map_difference(a, b, show_stat=False, vmin=-30., vmax=30., dmin=-60., dmax=60.)  # unprojected
plt.show()
