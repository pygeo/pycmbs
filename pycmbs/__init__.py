# -*- coding: utf-8 -*-

#details how to set global SVN variables can be found in http://blogchuck.com/2009/09/adding-svn-headers-revisited/

"""
# Copyright (C) 2012-2013 Alexander Loew, alexander.loew@mpimet.mpg.de
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
"""

__name__ = "pycmbs"
"""The project name."""

__author__ = "$Author: $"
"""The primary author of pyCMBS."""

__institute__ = "Max-Planck-Institute for Meteorology (MPI-M)"

__copyright__ = "Copyright (c) 2011-2013 Alexander Loew"
"""The copyright holder of pyCMBS."""

__license__ = "GNU General Public License"
"""The license governing the use and distribution of pyCMBS."""

__url__ = "https://code.zmaw.de/projects/pycmbs"
"""The URL for pyCMBS's homepage."""

__version__ = "0.2.0"
"""Version number of pyCMBS."""

__revision__ = filter(str.isdigit, "$Revision$")
"""SVN repository number, 
set SVN using svn propset
(see here: http://stackoverflow.com/questions/1449935/getting-svn-revision-number-into-a-program-automatically)"""

__date__ = "2013-06-14"
"""The release date of this version of pyCMBS."""

__docformat__ = "epytext en"
"""The epydoc documentation format for this file."""

__email__ = "alexander.loew@mpimet.mpg.de"

# set globally plotting backend
import matplotlib
matplotlib.use('Agg') 
