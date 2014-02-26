"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the files
LICENSE.md and COPYRIGHT.md
"""

"""
development script for pattern correlation analysis
"""

from pycmbs.diagnostic import PatternCorrelation
from pycmbs.data import Data
import numpy as np

import matplotlib.pyplot as plt

plt.close('all')
fname = '../pycmbs/examples/example_data/air.mon.mean.nc'

# generate two datasets
x = Data(fname, 'air', read=True)
xc = x.get_climatology(return_object=True)
yc = xc.copy()
yc.data = yc.data * np.random.random(yc.shape)*10.

PC = PatternCorrelation(xc, yc)
PC.plot()

plt.show()
