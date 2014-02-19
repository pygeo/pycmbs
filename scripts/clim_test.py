
from pycmbs.data import Data
import numpy as np

fname = '../pycmbs/examples/example_data/air.mon.mean.nc'

d=Data(fname, 'air', read=True)

c=d.get_climatology(return_object=True)

# now manipulate timeseries such that it is not starting in january any more
t = c.time*1.
c.time[2:] = t[0:-2]
c.time[0:2] = t[-2:]

di = np.diff(c.time)
m = di < 0.
if m.sum() == 0:
    pass
elif m.sum() == 1:  # a single breakpoint
    n = m.argmax()+1 # position where the break in timeseries occurs
else:
    raise ValueError('More than a single breakpoint found. Can not process this data as it is not in cyclic ascending order')



c1 = c.copy()

print c.date
c.timeshift(n, shift_time=True)

print ''
print c.date

print ''
print c1.date
