Installation checklist
======================

Here all obligatory and recommended dependencies are summarized.

This looks first of all like a huge list, but if you use a package manager, most of the dependencies are resolved automatically anyway.
The installation scripts provided in the installation instructions show that a fully installation can be done with a few commands only!

Core
----

python 2.7.x
    all development was done with python 2.7.x. Lower version numbers might work, but have not been tested so far.

python-dev
    an installation of the python development headers is also required

numpy
    python numpy is required for data handling

scipy >0.11
    python scipy is used for data analysis

netCDF4 library
    The netCDF4 library is the major backend for access to data files.

cdo >1.2
    climate data operators are required in the benchmarking environment for efficient data preprocessing and are also used for proper calculation of grid cell areas. The latter is important in many analysis functions to have a proper area weighting.

hdf5 library
    To access nc4 and hdf5 files, an installation of the HDF5 library is required. This is mandatory for usage of the CDO's. See cdo installation instructions. (libhdf5-openmpi-dev)

netcdf library
    also the installation of the netcdf library is required (libnetcdf-dev)

openmpi library
    the developer version of the openmpi library is required as well for the cdo's (libopenmpi-dev)

cdo python interface
    The python interface is used to make the cdo's and pyCMBS talk to each other

cython
    cython allows for very fast data processing from python code. It is used by some of the analysis functions and in particular for the netCDF4 library

latex
    a working installation of latex is required for generating automatic reports in the benchmarking framework. In particular your latex installation needs to support the *pdflatex* command.

nosetests
    nosetests are mandatory to run all the unittests in the development enviroment

Recommended
-----------

matplotlib >1.3
    some plotting features in pyCMBS will not work with older versions of matplotlib. Installation of a version > 1.3 is therefore highly recommended

freetype6 library
    library required by matplotlib. The development version (libfreetype6-dev) is required. See matplotlib installation instructions.

png library
    matplotlib also required the png library (libpng-dev). See matplotlib installation insstructions.

matplotlib basemap
    The Basemap module is required for plotting data on a map.

matplotlib basemap data
    Some additional data for matplotlib is required (e.g. background images)

cartopy
    cartopy is an efficient package for plotting data on a map

cartopy dependencies
    todo

pyshp
    library to read shapefiles using python





