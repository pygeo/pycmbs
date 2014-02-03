"""
Test map layouts
"""
from pyCMBS import *
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
file='testdata.nc'
d=Data(file,'tmp',read=True)

#~ d=Data(None,None)
#~ x = np.random.random((100,100))
#~ d.data = np.ma.array(x, mask = x != x)

# map only with colorbars
m = SingleMap(d)
m.plot(colorbar_orientation='horizontal', vmin=10., vmax=20.)


divider = make_axes_locatable(m.pax)



m1 = SingleMap(d)
m1.plot(colorbar_orientation='vertical')


# map with colorbar and zonal plot
m2 = SingleMap(d)
m2.plot(colorbar_orientation='horizontal', show_zonal=True, nclasses=4)


m3 = SingleMap(d, savefile='mym3')
m3.plot(colorbar_orientation='vertical', show_zonal=True, cmap='RdBu', nclasses=7, title='Mytest', ctick_prop={'ticks':[-15, 0., 3.], 'labels':['A','B','C'] })
#~ m3.pax.imshow(x)
#~ m3.cax.imshow(x)



plt.show()
