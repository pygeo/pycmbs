# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""
from pycmbs.data import Data
from pycmbs.plots import HovmoellerPlot
import matplotlib.pyplot as plt
import numpy as np
import datetime

file_name = '../../../pycmbs/examples/example_data/air.mon.mean.nc'
A = Data(file_name, 'air', lat_name='lat', lon_name='lon', read=True, label='air temperature')

# a quick plot as well as a projection plot
H = HovmoellerPlot(A)
H.plot(climits=(-20., 20.))
H.ax.set_yticks([])
H.ax.set_ylabel('lat')
plt.show()
