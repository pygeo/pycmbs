"""
working with the data object
"""
from pyCMBS import *
import os
import numpy as np
from matplotlib import pylab as pl

pl.close('all')


D = Data('air.mon.mean.nc', 'air',read=True)
P = Data('pr_wtr.eatm.mon.mean.nc','pr_wtr',read=True)


print 'Calculate differences ...'
map_plot(D.sub(P)) # = D-P

print 'Temporal evolution of mean and stdv of a field; note that values are calculated as area weigted values'
f=pl.figure();ax=f.add_subplot(111)
ax.plot(D.fldmean(),label='fldmean()')
ax.plot(D.fldstd(),label='fldstd()')
ax.legend(); ax.grid()


print 'Perform EOF and plot results ...'
C = D.get_climatology(return_object=True)
E = EOF(C) #calculate EOF based on climatology because of performance issues for this example.
E.plot_EOF(0,show_coef=True,use_basemap=True)
E.plot_EOF(1,show_coef=True,use_basemap=True)
print '... note that map_plot() arguments can be used here.'

pl.show()

r=raw_input("Press Enter to continue...")

pl.close('all')

#partial correlation analysis

#todo


#~ simple difference operations etc.
