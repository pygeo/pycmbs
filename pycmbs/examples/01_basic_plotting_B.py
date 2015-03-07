# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

"""
Basic plotting in pyCMBS
"""
from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt
import numpy as np
from pycmbs.examples import download

# Read some data ...
air = download.get_sample_file(name='air')

f=plt.figure(figsize=(10,6))
ax1=f.add_subplot(221)
ax2=f.add_subplot(222)
ax3=f.add_subplot(223)
ax4=f.add_subplot(224)

f1 = map_plot(air, use_basemap=False, title='vmin=-30.,vmax=30.,cmap_data=RdBu_r', vmin=-30., vmax=30., cmap_data='RdBu_r', ax=ax1)
#~ f2 = map_plot(air, contours=True, use_basemap=True, levels=np.arange(-50.,50.,10.), title='contours=True', ax=ax2)
f3 = map_plot(air, show_stat=True, use_basemap=False,title='show_stat=True',ax=ax3)
f4 = map_plot(air, show_stat=True, stat_type='median', use_basemap=False, title='show_stat=True,stat_type="median"', ax=ax4)

ax2.set_frame_on(False)
ax2.set_xticks([])
ax2.set_yticks([])
plt.show()
