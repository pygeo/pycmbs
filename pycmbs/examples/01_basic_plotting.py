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

# Read some sample data ...
air = download.get_sample_file(name='air')
air.label = 'air temperature'

# a quick plot as well as a projection plot
f1 = map_plot(air)  # unprojected
f2 = map_plot(air, use_basemap=True)  # projected
plt.show()




