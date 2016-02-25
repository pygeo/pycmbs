# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""


"""
module to compile the required python extensions
This is for development purposes only! Later on
it might be integrated into the standard setup.py
"""

# http://docs.cython.org/src/tutorial/cython_tutorial.htmlfrom
# distutils.core import setup
from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(
        ["./pycmbs/geostatistic/variogram_base.pyx"]),
    # this is needed to get proper information on numpy headers
    include_dirs=[numpy.get_include()]
)

# run as ... to build extension
# $ python setup_extensions.py build_ext --inplace
