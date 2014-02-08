from distutils.core import setup
import os

import pycmbs

install_requires = ["numpy>0.1", "cdo>1.2", ]

setup(name='pycmbs',
      version=pycmbs.__version__,
      packages=['pycmbs', 'pycmbs/benchmarking', 'pycmbs/tests',
                'pycmbs/benchmarking/logo', 'pycmbs/examples'],
      package_dir={'pycmbs': 'pycmbs'},
      package_data={'pycmbs': ['benchmarking/configuration/*',
                               'benchmarking/logo/*', 'LICENSE']},
      author="Alexander Loew",
      author_email='alexander.loew@mpimet.mpg.de',
      maintainer='Alexander Loew',
      maintainer_email='alexander.loew@mpimet.mpg.de',
      url='https://github.com/pygeo/pycmbs',
      description='pyCMBS - python Climate Model Benchmarking Suite',
      long_description='The pyCMBS project is a suite of tools to \
                        process, analyze, visualize and benchmark \
                        scientific model output against each other or \
                        against observational data. It is in particular \
                        useful for analyzing in an efficient way output \
                        from climate model simulations.',
      install_requires=install_requires,
      keywords=["data", "science", "climate", "meteorology",
                "model evaluation", "benchmarking", "metrics"],
      scripts=["pycmbs-benchmarking.py"],
      license="XXXX")

########################################################################
# Some useful information on shipping packages
########################################################################

#PIP
#~ python setup.py register
#~ python setup.py sdist
#~ python setup.py upload

#http://docs.python.org/2/distutils/setupscript.html#installing-package-data
#http://docs.python.org/2/distutils/setupscript.html#installing-additional-files
