# -*- coding: UTF-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

# good introduction into packing can be found in
# https://python-packaging-user-guide.readthedocs.org/en/latest/index.html

# the setuptools are supposed to be used as a standard. Thats why we ommit
# usage of distutils here

# example of setup.py can be found here:
# https://github.com/pypa/sampleproject/blob/master/setup.py


from setuptools import setup

import pycmbs
from Cython.Build import cythonize
import os

install_requires = ["numpy>0.1", "cdo>1.2", "netCDF4", "pytz", "matplotlib"]

setup(name='pycmbs',

    version=pycmbs.__version__,

    description='pyCMBS - python Climate Model Benchmarking Suite',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    # packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    packages=['pycmbs', 'pycmbs/benchmarking', 'pycmbs/tests',
            'pycmbs/benchmarking/logo', 'pycmbs/examples', 'pycmbs/diagnostic', 'pycmbs/colormaps', 'pycmbs/plots'],
    package_dir={'pycmbs': 'pycmbs'},
    package_data={'pycmbs': ['benchmarking/configuration/*',
                           'benchmarking/logo/*', 'version.json']},


    author="Alexander Loew",
    author_email='alex@geo2-consult.de',
    maintainer='Alexander Loew',
    maintainer_email='alex@geo2-consult.de',

    license='MIT',

    url='https://github.com/pygeo/pycmbs',

    long_description='The pyCMBS project is a suite of tools to \
                    process, analyze, visualize and benchmark \
                    scientific model output against each other or \
                    against observational data. It is in particular \
                    useful for analyzing in an efficient way output \
                    from climate model simulations.',

    # List run-time dependencies here. These will be installed by pip when your
    # project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/technical.html#install-requires-vs-requirements-files
    install_requires=install_requires,

    keywords=["data", "science", "climate", "meteorology",
            "model evaluation", "benchmarking", "metrics"],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.

    entry_points={
        'console_scripts': [
            'pycmbs-benchmarking = pycmbs-benchmarking:main'
        ]},

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
    # How mature is this project? Common values are
    # 3 - Alpha
    # 4 - Beta
    # 5 - Production/Stable
    'Development Status :: 4 - beta',
    # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Topic :: Scientific/Engineering :: GIS',
    'Topic :: Scientific/Engineering :: Visualization',

    # Pick your license as you wish (should match "license" above)
    'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 2.7'
    ]
)






    #~ entry_points={
        #~ 'console_scripts': [
            #~ 'foo = my_package.some_module:main_func',
            #~ 'bar = other_module:some_func',
        #~ ],
        #~ 'gui_scripts': [
            #~ 'baz = my_package_gui:start_func',
        #~ ]
    #~ }












#~ setup_dist(
  #~ ext_modules = cythonize("./pycmbs/polygon_utils.pyx"),
#~ )


# compile extensions
#~ os.system('sh compile_extensions.sh')


########################################################################
# Some useful information on shipping packages
########################################################################

# PIP
#~ python setup.py register
#~ python setup.py sdist
#~ python setup.py upload
