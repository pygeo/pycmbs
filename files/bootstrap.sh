#!/usr/bin/env bash

#
# This file provides an installing procedure for pyCMBS WITHOUT Cartopy support
#    it was tested for ubuntu32
#

# update package database
apt-get update

#####################################################################
# DEPENDENCIES
#####################################################################

# the -qq option installs silent using defaults
apt-get -qq install texlive-latex-base texlive-latex-extra texlive-latex-recommended
apt-get -qq install python-pip python-dev
apt-get -qq install cdo libhdf5-openmpi-dev libnetcdf-dev libopenmpi-dev
apt-get -qq install python-numpy
apt-get -qq install cython
C_INCLUDE_PATH=/usr/include/mpi pip install netCDF4

# apt-get -qq install python-matplotlib  # this gives the system default package, which is currently v1.1.1 therefore it is not used here
# it is highly recommended to use matplotlib > 1.3
apt-get -qq install libfreetype6-dev libpng-dev  # required for matplotlib
sudo easy_install -U distribute
sudo pip install https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.3.1/matplotlib-1.3.1.tar.gz
apt-get -qq install python-mpltoolkits.basemap
apt-get -qq install python-mpltoolkits.basemap-data

apt-get -qq install python-scipy

#####################################################################
# pycmbs
#####################################################################
pip install --upgrade pycmbs

#####################################################################
# test environment
#####################################################################
pip install nose

echo "Now you can run the unittests as follows:"
echo "    cd /usr/local/lib/python2.7/dist-packages/pycmbs/tests"
echo "    nosetests"

