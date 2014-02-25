"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the files
LICENSE.md and COPYRIGHT.md
"""

"""
Hstackplotting
"""

from pycmbs.plots import HstackTimeseries
import numpy as np
import matplotlib.pyplot as plt


plt.close('all')

ht = HstackTimeseries()

for i in xrange(50):
    x = np.random.random(12)*2.-1.
    ht.add_data(x, 'model' + str(i).zfill(3) )

ht.plot(cmap='RdBu_r', interpolation='nearest', vmin=-1., vmax=1., nclasses=15, title='Testtitle', maxheight=0.3, monthly_clim_ticks=True)



plt.show()


#~ set_size_inches(18.5,10.5)
