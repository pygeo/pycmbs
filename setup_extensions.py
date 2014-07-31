"""
module to compile the required python extensions
This is for development purposes only! Later on
it might be integrated into the standard setup.py
"""

# http://docs.cython.org/src/tutorial/cython_tutorial.htmlfrom distutils.core import setup
from distutils.core import setup
from Cython.Build import cythonize

setup(
  ext_modules = cythonize("./pycmbs/polygon_utils.pyx"),
)

# run as ...
# $ python setup.py build_ext --inplace
