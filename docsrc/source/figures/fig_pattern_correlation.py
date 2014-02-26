from pycmbs.data import Data
from pycmbs.diagnostic import PatternCorrelation
import matplotlib.pyplot as plt
import numpy as np

file_name = '../../../pycmbs/examples/example_data/air.mon.mean.nc'
A = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')
B = A.copy()
B.mulc(2.3, copy=False)
B.data = B.data + np.random.random(B.shape)*100.

# calculate spatial correlation for all timesteps ...
P = PatternCorrelation(A,B)
# ... and vizalize it
P.plot()

plt.show()
