from pycmbs.data import Data
from pycmbs.plots import GlobalMeanPlot
import matplotlib.pyplot as plt

file_name = '../../../pycmbs/examples/example_data/air.mon.mean.nc'
A = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')
B = A.copy()
B.mulc(2.3, copy=False)

# a quick plot as well as a projection plot
GM = GlobalMeanPlot()
GM.plot(A, label='model A')
GM.plot(B, label='model V')
GM.ax.grid()
GM.ax1.grid()
plt.show()
