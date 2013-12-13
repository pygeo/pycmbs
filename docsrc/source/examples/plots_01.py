"""
Generate some sample map plots
"""
import matplotlib.pyplot as plt
from pyCMBS import *

example_dir='../../../pyCMBS/examples/example_data'
file =example_dir+'/air.mon.mean.nc'

# This is how you read a sample datafile
# Simply specify the filename, the variable name (here: "air") and set read=True
D = Data(file, 'air', read=True)

# some simple plotting
f = map_plot(D, title='This is the simplest and fastest plot')
f = map_plot(D,use_basemap=True, title='The same as projected figure')

# some further options
f = map_plot(D,show_timeseries=True,use_basemap=True,title='show_timeseries=True')
f = map_plot(D,show_zonal=True,use_basemap=True,title='show_zonal=True')
f = map_plot(D,show_histogram=True,use_basemap=True,title='show_histogram=True')


# now we create a figure with several subplots and illustrate further different options
# note the usage of "ax" in map_plot()
f = pl.figure()
ax1 = f.add_subplot(221)
ax2 = f.add_subplot(222)
ax3 = f.add_subplot(223)
ax4 = f.add_subplot(224)

map_plot(D,use_basemap=True,title='vmin=-30.,vmax=30.,cmap_data=RdBu_r',vmin=-30.,vmax=30.,cmap_data='RdBu_r',ax=ax1)
map_plot(D,contours=True,use_basemap=True,levels=np.arange(-50.,50.,10.),title='contours=True',ax=ax2)
map_plot(D,show_stat=True,use_basemap=True,title='show_stat=True',ax=ax3)
map_plot(D,show_stat=True,stat_type='median',use_basemap=True,title='show_stat=True,stat_type="median"',ax=ax4)

plt.show()

