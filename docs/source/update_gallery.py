"""
update the gallery plots
"""

from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt
from pycmbs.examples import download
import numpy as np

# Read some sample data ...
air = download.get_sample_file(name='air')
air.label = 'air temperature'

# a quick plot as well as a projection plot
f1 = map_plot(air, title='unprojected data')  # unprojected
f2 = map_plot(air, use_basemap=True, title='projected data')  # projected

f1.savefig('./fig/figure1.png')
f2.savefig('./fig/figure2.png')

f = map_plot(air,show_zonal=True, use_basemap=True,title='show_zonal=True')
f.savefig('./fig/figure3.png')
del f


f = plt.figure()
ax1 = f.add_subplot(221)
ax2 = f.add_subplot(222)
ax3 = f.add_subplot(223)
ax4 = f.add_subplot(224)

map_plot(air, use_basemap=True, title='vmin=-30.,vmax=30.,cmap_data=RdBu_r', vmin=-30., vmax=30., cmap_data='RdBu_r', ax=ax1)
ax2.set_xticks([])
ax2.set_yticks([])
map_plot(air, show_stat=True, use_basemap=True,title='show_stat=True',ax=ax3)
map_plot(air, show_stat=True, stat_type='median', use_basemap=True, title='show_stat=True,stat_type="median"', ax=ax4)
f.savefig('./fig/figure4.png')
