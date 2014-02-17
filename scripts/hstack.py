"""
Hstackplotting
"""

from pycmbs.plots import HstackTimeseries
import numpy as np
import matplotlib.pyplot as plt


plt.close('all')

ht = HstackTimeseries()

for i in xrange(5):
    x = np.random.random(100)*2.-1.
    ht.add_data(x, 'model' + str(i).zfill(3) )

ht.plot(cmap='RdBu_r', interpolation='nearest', vmin=-1., vmax=1., nclasses=15)

plt.show()
