Installation
============

Below we describe how you properly install pyCMBS in your environment.
Currently, there are three different ways to install pyCMBS

 1. NOT RECOMMENDED AT THE MOMENT Easy installation using *pip* (note: not checked for a while be careful!)
 2. Source code installation from code repository (currently recommended)
 3. Source code installation from tarball

All approaches are detailed below. Special informations for users working
at the Max-Planck-Institute for Meteorology are provided in a
:doc:`installation_mpi`

Operating systems
-----------------

pyCMBS is developed purely in python and is expected to be in general independent from a special operating system.
The current version was however only developed on Linux and no tests at all were made on other operating systems.

Dependencies
------------

pyCMBS was built to have only a small number of external dependencies.
However, there is a minimum number of dependencies existing which are a
required for using pyCMBS sucessfully.


*Core python packages [obligatory]*

- python 2.7.x
- `matplotlib v1.3.x <http://matplotlib.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_

Note that you should really ensure that your matplotlib installation is greater than version 1.3, as otherwise some features will not work. Standard packages like e.g. shipped as standard with Ubuntu have smaller version numbers (e.g. v1.1.1).

*Obligatory additional dependencies*

For file I/O of netCDF files, a single library is supported at the moment:

- `netCDF4 library <http://code.google.com/p/netcdf4-python/>`_ is used as default for file I/O [recommended]

For an efficient data pre-processing the climate data operators are used. The core CDO's and the corresponding python wrapper is required.

- `climate data operators (cdo) <https://code.zmaw.de/projects/cdo>`_ [obligatory]
- `cdo python interface <https://code.zmaw.de/projects/cdo/wiki/Cdo%7Brbpy%7D>`_ [obligatory, can be easily installed using pip]

*Recommended dependencies*

For plotting projected data (map plots), pyCMBS currently supports two different plotting backends. These are

- `cartopy <http://scitools.org.uk/cartopy/>`_ [recommended]
- `matplotlib basemap <http://matplotlib.org/basemap/index.html>`_ [optional]


Installation of dependencies
----------------------------

For convenience, we have put together an installation procedure to install all required dependencies (except the ones needed for cartopy).
Please find it here__ for your convenience.

__ installation_dependencies_

An :doc:`installation_checklist` summarizing the required packages is available as well.

Quick installation from scratch for experts
-------------------------------------------

If you have not yet any of the above dependencies installed and are
working on a Debian like operating system (e.g. Ubunutu), the easiest way to
install pyCMBS is by executing the installation sequence which can be found in
the file *.travis.yml* in the root directory of the source code.

The developers are using automatic code checking and builds and the
installation sequence in the file *.travis.yml* is used to setup for each build
a working installation.


Detailed installation instructions for pyCMBS
---------------------------------------------

In the following, we will summarize the different approaches to install pyCMBS.

Installation using *pip* (the easiest way ... theoretically)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NOT RECOMMENDED AT THE MOMENT AS NOT TESTED!!!
USE INSTALLATION FROM GITHIB REPOSITORY INSTEAD!!!

Using `pip <https://pypi.python.org/pypi/pip>`_ is in general the easiest way to install pyCMBS if you dont want to develop by your own.

NOTICE however that this installation procedure was not tested thoroughly with present version. Please be therefore careful.

The installation with *pip* was tested without the cartopy plotting backend so far.

If you have *pip* not yet installed on your machine, then the first step
is to `install pip <https://pypi.python.org/pypi/pip>`_ .

Once *pip* is installed, the installation of pyCMBS is as simple as a two-liner::

    pip install numpy
    pip install pycmbs

Note that the *numpy* module needs to be installed separately before as otherwise the installation fails. The command *pip install pycmbs* will then install all remaining dependencies.

Check if everything is working, like described below_.


Installation of stable version from a *tarball*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can install specific version of pycmbs from a *tarball* from the `download site <https://github.com/pygeo/pycmbs/releases>`_ . All dependencies must have been installed before.

After you have obtained the pyCMBS code from a *tarball*, then first extract the archive into some new directory::

    mkdir temp_dir
    cd temp_dir
    tar -xvf pycmbs-vx.x.x.tar.gz

This expands the tar file and you obtain a subdirectory,  called *pycmbs-vx.x.x*

To install the package first change to the directory and then install
the package using the standard python setup tools as::

    cd pycmbs-vx.x.x
    python setup.py install

This will install the package in your python environment.
Check successful installation, like described below_.


github repository (for developers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another alternative to work with pyCMBS is to used the master branch of
the `code repository <https://github.com/pygeo/pycmbs>`_. The code is hosted on a github repository. All
dependencies  must have been installed before.

It is recommended that you first fork the entire project on github to make your own sub-project.
This will allow you also to make changes in the code and to contribute to the further development.
The description below applies to both, forked projects as well as the main project branch.

To retrieve the code do the following::

    # generate some directory
    mkdir pycmbs
    cd pycmbs

    # retrieve the code
    git clone

*Compilation*

Some sub-modules are written in cython to speed up processing. These modules need to be compiled prior to the final installation. This is done by just executing the following command::

    # compile cython code
    sh compile_extensions.py

*Final installation*

Now you have in principle two options. You either decide that the code should be installed in the python dist-packages directory, then you do::

    # do installation
    python setup.py install

or if you want to hack the code, it is highly recommended to simply set the right PYTHONPATH environment variable::

    # or as an alternative for developers, just set the PYTHONPATH
    # environment variable to the pycmbs root directory and also adapt
    # you systempath (PATH) such that includes the pycmbs rootdirectory

Check successful installation, like described below_.


.. _installation_dependencies:

Installation of dependencies
----------------------------

Please find here a working installation procedure, which was tested under Ubuntu 32-bit. It installs all pyCMBS dependencies, except the ones needed for cartopy and installs pyCMBS itself.::

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
    pip install pyshp

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


.. _below:

Final check of installation
---------------------------

Check that installation worked properly by going through the following
checklist. In case of problems, please refer to the troublesolver_ .

*Is the pyCMBS python module loaded properly?*::

    python -c "from pycmbs import *; print('Welcome to pyCMBS')"

This should give you a short welcome message, but no error messages.

*Is the benchmarking script working properly?*::

    pycmbs-benchmarking.py

This will you give a short message like::

   *******************************************
   * WELCOME to pycmbs.py                    *
   * Happy benchmarking ...                  *
   *******************************************

and will end with an error message that the configuration file is
not found (this is o.k.)

**If you see the above, the installation has worked! Congratulations!**

3. Check also the proper installation of the cdo's and the cdo.py
interface, as this is a prerequesite of beeing able to properly work
with pyCMBS::

     python -c "from cdo import *; cdo=Cdo(); print 'If you see this, everything went right ... have fun with pyCMBS and CDOs'"

Again, this should give you a short welcome message. Any error message
is a bad sign. In that case, please check your installation again.
Have a look at the troublesolver_.

pycmbs init


Running tests
-------------

pyCMBS code comes with a rich suite of test routines. We follow the concept of unittests using the nosetests tools. Tests should be always executed in the following cases:

* after installation
* before comitting code to the repository
* before merging a branch in the master branch

Tests can be simply executed using the *Makefile* in the main installation directory as::

    make tests

As an alternative you can also check the coverage of tests in the code using::

    make coverage

which gives you a report on test coverage in */coverage/index.html*.


.. _installation_details:

Further information and trouble solving
---------------------------------------

pyCMBS makes use of a standard directory to look for observations. This
directory is the Standard Evaluation Pool (SEP). The path to the SEP directory
needs to be specified in the $SEP environment variable. In you .bashrc write::

    export SEP=/path/to/directory

For users at MPI-M, the SEP variable needs to point to */pool/SEP*. It is
however possible to specify also for each observation an individual path where
the observation is located. Then the SEP evnironment variable is not required.
To check whether SEP is set, type::

    echo $SEP

.. _troublesolver:

Some hints for trouble solving
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your pyCMBS installation seems not to work properly, here are a few
recommendations where to start searching.

*Is python working properly?*::

    python -c "print 'Hello world'"

*Does your PYTHONPATH environment variable contain the path to pyCMBS?*::

    echo $PYTHONPATH

This should give you the path where python is searching for modules.
If it is empty you are most likely in trouble. Check if you have a
valid python installation.

*Is the script pycmb-benchmarking.py found in the system path?*::

    pycmbs-benchmarking.py

should give you a short Welcome Screen like described above. If this is not the
case then either the overall pyCMBS installation is incomplete or Your
systempath is not set appropriately. Type::

    echo $PATH

and verify if the directory where pycmbs-benchmarking.py is located is listed
in your PATH variable. If not, then you can try to change your PATH variable to
make it working.

*Further problems?*

In case that these recommendations did not solve your problem, please
feel free to ask a question or raise an issue on the pyCMBS `development site
<https://github.com/pygeo/pycmbs>`_.










