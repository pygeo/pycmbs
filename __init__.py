#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# pyCMBS -- Pythonic Climate Model Suite
#
#
# Copyright (C) 2011-2012 Alexander Loew 
# Authors: Alexander Loew <alexander.loew@zmaw.de.de>,
#          Michael Borsche,
#          Mikhail Itkin,
#          Gernot Geppert
#
# Max-Planck-Institute for Meteorology, Hamburg, Germany
#

__name__ = "pyCMBS"
"""The project name."""

__author__ = "Alexander Loew <alexander.loew@zmaw.de>"
"""The primary author of pyCMBS."""

__copyright__ = "Copyright (c) 2011-2012 Alexander Loew"
"""The copyright holder of pyCMBS."""

__license__ = "MIT Open Source License"
"""The license governing the use and distribution of pygeonetwork."""

__url__ = "http://www.mpimet.mpg.de/en/science/the-land-in-the-earth-system/terrestrial-remote-sensing-hoaps.html"
"""The URL for pyCMBS's homepage."""

__version__ = 0.1
"""Version number of pyCMBS."""

__date__ = "2012-03-14"
"""The release date of this version of pyCMBS."""

__docformat__ = "epytext en"
"""The epydoc documentation format for this file."""

#------------------------------
#--- Import classes
#------------------------------

#- generic classes
from matplotlib import pylab as pl

from scipy import stats

import numpy as np

from mpl_toolkits.basemap import Basemap,shiftgrid

import os

import sys

#- classes specific to pyCMBS
from statistic import *
from data   import *
from report import *
from region import *
from diagnostic import *
from plots  import *
from pyCDO  import *





#------------------------------
#--- Global constants
#------------------------------
