# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
EOF analysis
"""

from pycmbs.examples import download
import matplotlib.pyplot as plt
from pycmbs.diagnostic import EOF

air = download.get_sample_file(name='air')
air.label = 'air temperature'

# calculate climatological mean
clim = air.get_climatology(return_object=True)

# calculate EOF based on climatology because of performance issues for this example.
E = EOF(clim)
E.plot_EOF(0,show_coef=True, use_basemap=True)  # map_plot argument can be used here
E.plot_EOF(1,show_coef=True, use_basemap=True)

plt.show()
