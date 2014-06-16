# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os
import json


ppath = os.path.dirname(os.path.realpath(__file__))

__name__ = "pycmbs"
"""The project name."""

__author__ = "Alexander Loew"
"""The primary author of pyCMBS."""

__institute__ = "Max-Planck-Institute for Meteorology (MPI-M)"

__copyright__ = "Copyright (c) 2011-2014 Alexander Loew"
"""The copyright holder of pyCMBS."""

__license__ = "MIT License, see LICENSE.md for details"
"""The license governing the use and distribution of pyCMBS."""

__url__ = "https://github.com/pygeo/pycmbs"
"""The URL for pyCMBS's homepage."""

__date__ = "2014-06-05"
"""The release date of this version of pyCMBS."""

__version__ = json.load(open(ppath + os.sep + 'version.json'))  # note that a proper version number is required as PIP otherwise can not detect the file
"""Version number of pyCMBS."""

__email__ = "alexander.loew@mpimet.mpg.de"

# set globally plotting backend
try:
    import matplotlib
    matplotlib.use('Agg')
except:
    print('Some problem with import of matplotlib happened')

try:
    from mapping import MultipleMap, SingleMap, MapPlotGeneric
except:
    pass

# set automatically directory where pycmbs is located
os.environ.update({'PYCMBSPATH': os.path.dirname(os.path.realpath(__file__))})
