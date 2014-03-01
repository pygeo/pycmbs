# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Basic plotting in pyCMBS
"""

from pycmbs.mapping import map_plot
import matplotlib.pyplot as plt
from pycmbs.examples import download

air = download.get_sample_file(name='air')
air.label = 'air temperature'
f1 = map_plot(air, show_timeseries=False, use_basemap=True, title='show_timeseries=True')
f2 = map_plot(air, show_zonal=True, use_basemap=True, title='show_zonal=True')
f3 = map_plot(air, show_histogram=True, use_basemap=True, title='show_histogram=True')
plt.show()
