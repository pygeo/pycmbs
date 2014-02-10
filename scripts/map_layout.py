"""
Test map layouts
"""
from pycmbs.mapping import SingleMap
from pycmbs.data import Data
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as grd
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.close('all')


#~ fig = plt.figure()
#~
#~ gs = grd.GridSpec(2, 1, height_ratios=[95,5], wspace=0.05)
#~
#~ pax = fig.add_subplot(gs[0])
#~ cax = fig.add_subplot(gs[1])
#~
#~ cax.set_xticklabels([])
#~ cax.set_xticks([])
#~
#~ plt.show()
#~
#~ stop



#~ file='/home/m300028/shared/data/SEP/variables/land/Ta_2m/cru_ts_3_00.1901.2006.tmp_miss_t63.nc'

f=plt.figure()
ax1=f.add_subplot(2,1,1)

file='testdata.nc'
d=Data(file,'tmp',read=True)

map_plot(d, ax=ax1, colorbar_orientation='vertical', nclasses=34, cticks=[-20.,0., 10.], cticklabels=['A','B','C'], show_zonal=True, use_basemap=True)


plt.show()
stop

#~ d=Data(None,None)
#~ x = np.random.random((100,100))
#~ d.data = np.ma.array(x, mask = x != x)


# map only with colorbars
m = SingleMap(d, backend='basemap')
m.plot(colorbar_orientation='horizontal', vmin=10., vmax=20., proj_prop={'projection':'robin', 'lon_0':0.})

mx = SingleMap(d, backend='imshow')
mx.plot(colorbar_orientation='horizontal', vmin=10., vmax=20.)


m1 = SingleMap(d)
m1.plot(colorbar_orientation='vertical')

m1 = SingleMap(d, backend='basemap')
m1.plot(colorbar_orientation='vertical', proj_prop={'projection':'robin', 'lon_0':0.})




# map with colorbar and zonal plot
m2 = SingleMap(d)
m2.plot(colorbar_orientation='horizontal', show_zonal=True, nclasses=4)


m3 = SingleMap(d, savefile='mym3')
m3.plot(colorbar_orientation='vertical', show_zonal=True, cmap='RdBu', nclasses=7, title='Mytest', ctick_prop={'ticks':[-15, 0., 3.], 'labels':['A','B','C'] })
#~ m3.pax.imshow(x)
#~ m3.cax.imshow(x)

m3 = SingleMap(d, savefile='mym3xx', backend='basemap')
m3.plot(colorbar_orientation='vertical', show_zonal=True, cmap='RdBu', nclasses=7, title='Mytest', ctick_prop={'ticks':[-15, 0., 3.], 'labels':['A','B','C'] }, proj_prop={'projection':'robin', 'lon_0':0.})


plt.show()
