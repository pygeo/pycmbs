"""
Basic working with the data object
"""

# XXX: implicit import
# from pyCMBS import *
import os
import numpy as np
from matplotlib import pylab as pl

pl.close('all')

file='../example_data/air.mon.mean.nc'
if not os.path.exists(file):
    raise ValueError, 'Sample file not existing: see example-01.py'

#--- read data ---
D = Data('../example_data/air.mon.mean.nc', 'air',read=True)
P = Data('../example_data/pr_wtr.eatm.mon.mean.nc','pr_wtr',read=True)

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
#get_climatology() returns 12 values which are then used for plotting
map_season(D.get_climatology(return_object=True),use_basemap=True,vmin=-20.,vmax=30.)

print 'Map difference between datasets ...'
map_difference(D,P)

print 'ZonalPlot ...'
Z=ZonalPlot(); Z.plot(D)

print 'Some LinePlot'
L=LinePlot(regress=True, title='This is a LinePlot with regression')
L.plot(D,label='2m air temperature')
L.plot(P,label='Precipitable water',ax=L.ax.twinx() ,color='green') #use secondary axis for plotting here
L.legend()

print 'Scatterplot between different variables ...'
S=ScatterPlot(D) #scatterplot is initialized with definition of X-axis object
S.plot(P)
S.legend()

print 'Hovmoeller diagrams ...'
hm = HovmoellerPlot(D)
hm.plot(climits=[-20.,30.])

print '... generate Hovmoeller plot from deseasonalized anomalies'
ha=HovmoellerPlot(D.get_deseasonalized_anomaly(base='all'))
ha.plot(climits=[-2.,2.],cmap='RdBu_r')


pl.show()
r=raw_input("Press Enter to continue...")

pl.close('all')


