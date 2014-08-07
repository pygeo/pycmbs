"""
This is an example that should illustrate how you can scale
a dataset by the length of the month
"""

from pycmbs.examples import download
from pycmbs.data import Data
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt

plt.close('all')

# read some data as Data object
filename = download.get_sample_file(name='air', return_object=False)
air = Data(filename, 'air', read=True)

# this dataset has the following times
print air.date

# obviously the different months have different numbers of days.
# Let's say you want now to perform a proper averaging of the data
# taking into account the different lengths of the months
#
# the way how you would do it is like
# y = sum(w[i] * x[i])
# whereas w is a weighting factor for each timestep and 'x' is the input data

# how can you easily do that with the Data object?

# 1) calculate the weights ...
#     these are dependent on the number of days  which you get as ...

ml = air._get_days_per_month()
print ml

w = ml / float(sum(ml))  # relative weights
print w
print 'The sum of the weights should be 1.: ', sum(w)

# 2) now we need to multiply the data with the weights and then
#  sum up over time
new = air.mul_tvec(w, copy=True)  # check also the other options of this routine!
res = new.timsum(return_object=True)  # gives a 2D array with weighted results
print res.shape

map_plot(air, title='classical unweighted temporal mean', show_stat=True)
map_plot(res, title='weighted temporal mean', show_stat=True)
map_plot(air.sub(res), title='difference', show_stat=True, cmap_data='RdBu_r', vmin=-0.1, vmax=0.1)

# you will see that the differences are marginal, but they are nevertheless there
# for this example.

plt.show()






