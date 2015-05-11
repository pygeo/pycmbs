# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""


"""
working with the data object
"""


import os
import numpy as np
import matplotlib.pyplot as plt

from pycmbs.examples import download
from pycmbs.mapping import map_plot

plt.close('all')

air = download.get_sample_file(name='air')
air.label = 'air temperature'

# generate some artificial data by rescaling the original data
air2 = air.copy()  # copy data object
air2.mulc(1.3, copy=False)  # copy=False applies operation directly on data, otherwise a copy is returned

# difference
print 'Calculate differences ...'
diff = air.sub(air2)  # calculate difference
f1 = map_plot(diff, cmap_data = 'RdBu_r', vmin=-10., vmax=10.)

# temporal evolution of mean and stdv
print 'Temporal evolution of mean and stdv of a field; note that values are calculated as area weigted values'
f = plt.figure()
ax = f.add_subplot(111)
ax.plot(air.fldmean(), label='fldmean()')
ax.plot(air.fldstd(), label='fldstd()')
ax.legend()
ax.grid()

plt.show()


