"""
Basic working with the data object
"""

from pyCMBS import *
import os
import numpy as np
from matplotlib import pylab as pl

pl.close('all')

file='air.mon.mean.nc'
if not os.path.exists(file):
    raise ValueError, 'Sample file not existing: see example-01.py'

#--- read data ---
D = Data('./air.mon.mean.nc', 'air', lat_name='lat',lon_name='lon',read=True)

#--- some analysis ---
print 'Temporal stdv. ...'
t = D.timstd(return_object=True)
map_plot(t,use_basemap=True,title='Temporal stdv.',show_stat=True)

print 'Temporal trend ...'
f=pl.figure()
ax1=f.add_subplot(221)
ax2=f.add_subplot(222)
ax3=f.add_subplot(223)
ax4=f.add_subplot(224)
R,S,I,P = D.temporal_trend(return_object=True)
map_plot(R,use_basemap=True,ax=ax1)
map_plot(S,use_basemap=True,ax=ax2)
map_plot(I,use_basemap=True,ax=ax3)
map_plot(P,use_basemap=True,ax=ax4)
f.suptitle('Example of temporal correlation analysis results', size=20)


print 'Calculate climatology and plot ...'
map_season(D.get_climatology(return_object=True),use_basemap=True)


#map_difference

#LinePlot
#ScatterPlot
#...

#ZonalPlot







print 'Hovmoeller diagrams ...'
hm = HovmoellerPlot(D)
hm.plot(climits=[-20.,30.])




#EOF analysis

#todo

pl.show()
