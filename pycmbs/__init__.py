# -*- coding: utf-8 -*-

"""
Copyright (C) 2012-2014 Alexander Loew, alexander.loew@mpimet.mpg.de
See COPYING file for copying and redistribution conditions.
"""

__name__ = "pycmbs"
"""The project name."""

__author__ = "Alexander Loew"
"""The primary author of pyCMBS."""

__institute__ = "Max-Planck-Institute for Meteorology (MPI-M)"

__copyright__ = "Copyright (c) 2011-2014 Alexander Loew"
"""The copyright holder of pyCMBS."""

__license__ = "GNU General Public License"
"""The license governing the use and distribution of pyCMBS."""

__url__ = "https://code.zmaw.de/projects/pycmbs"
"""The URL for pyCMBS's homepage."""

__version__ = "1.0.0-dev"
"""Version number of pyCMBS."""

__date__ = "2013-06-14"
"""The release date of this version of pyCMBS."""

__email__ = "alexander.loew@mpimet.mpg.de"

# set globally plotting backend
import matplotlib
matplotlib.use('Agg')
from mapping import MultipleMap, SingleMap, MapPlotGeneric
