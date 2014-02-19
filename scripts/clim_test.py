
from pycmbs.data import Data
import numpy as np

fname = '../pycmbs/examples/example_data/air.mon.mean.nc'

d=Data(fname, 'air', read=True)
c=d.get_climatology(return_object=True)

print 'c raw: ', c.fldmean()
print c.date
print ''

# create some invalid data
d1=d.copy()
t = d1.time*1.
d1.time[20:] = t[0:-20]
d1.time[0:20] = t[-20:]

tmp = d1.data*1.
d1.data[20:,:,:] = tmp[0:-20,:,:]
d1.data[0:20,:,:] = tmp[-20:,:,:]

c1=d1.get_climatology(return_object=True, ensure_start_first=True)

print ''
print c1.date
print c1.fldmean()
