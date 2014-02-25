Installation
============

Below we describe how you properly install pyCMBS in your environment.
Currently, there are three different ways to install pyCMBS

 1. Easy installation using *pip* (recommended)
 2. Source code installation from code repository
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
- `matplotlib <http://matplotlib.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <http://www.scipy.org/>`_

*Obligatory additional dependencies*

For file I/O of netCDF files, there are in general two options supported at the moment. At least one interface library needs to be installed:

- `netCDF4 library <http://code.google.com/p/netcdf4-python/>`_ is used as default for file I/O [recommended]
- `Nio library <https://www.pyngl.ucar.edu/Nio.shtml>`_ could be installed as an alternative [optional, might be depreciated in the future]

For an efficient data pre-processing the climate data operators are used. The core CDO's and the corresponding python wrapper is required.

- `climate data operators (cdo) <https://code.zmaw.de/projects/cdo>`_ [obligatory]
- `cdo python interface <https://code.zmaw.de/projects/cdo/wiki/Cdo%7Brbpy%7D>`_ [obligatory, can be easily installed using pip]

*Recommended dependencies*

For plotting projected data (map plots), pyCMBS currently supports two different plotting backends. These are

- `cartopy <http://scitools.org.uk/cartopy/>`_ [recommended]
- `matplotlib basemap <http://matplotlib.org/basemap/index.html>`_ [optional]


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

Installation using *pip* (the easiest way)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using `pip <https://pypi.python.org/pypi/pip>`_ is the easiest (and recommended) way to install pyCMBS.
If you have *pip* not yet installed on your machine, then the first step
is to `install pip <https://pypi.python.org/pypi/pip>`_ .

Once *pip* is installed, the installation of pyCMBS is as simple as::

    pip install pycmbs

Check if everything is working, like described below_.


Installation of stable version from a *tarball*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can install specific version of pycmbs from a *tarball* from the `download site <https://code.zmaw.de/projects/pycmbs/files>`_ . All dependencies must have been installed before.

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

To retrieve the code do the following::

    # generate some directory
    mkdir pycmbs
    cd pycmbs

    # retrieve the code
    git clone

    # do installation
    python setup.py install

    # or as an alternative for developers, just set the PYTHONPATH
    # environment variable to the pycmbs root directory and also adapt
    # you systempath (PATH) such that includes the pycmbs rootdirectory

Check successful installation, like described below_.


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










