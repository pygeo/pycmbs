"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""
from pycmbs.data import Data
from pycmbs.examples import download
import matplotlib.pyplot as plt

plt.close('all')

# load some sample data
filename = download.get_sample_file(name='air', return_object=False)
x = Data(filename, 'air', read=True)
print 'Data dimensions: ', x.shape

# calculate global mean temperature timeseries
t = x.fldmean()

# plot results as a figure
f = plt.figure()
ax = f.add_subplot(111)
ax.plot(x.date, t, label='global mean')
ax.set_xlabel('Years')
ax.set_ylabel('Temperature [degC]')

# perhaps you also want to calculate some statistics like the temperature trend
from scipy import stats
import numpy as np
slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(x.time, t)
# note that the slope has the same units like the time variable of the Data object. Here it is hours!
# if we want to express the slope in [K/decade] we need to rescale
slope = slope * 24. * 365.25 * 10.
print 'Temperature trend [K/decade]: ', slope

ax.set_title('Global temperature trend: ' + str(np.round(slope,3)) + ' [K/decade]')

plt.show()
