Installation
============

Below we describe how you properly install pyCMBS in your environment. Currently, there are three different ways to install pyCMBS

 1. *ZMAW users* at KlimaCampus, Hamburg can use pre-installed setup
 2. Easy installation using *pip*
 3. Installation from *source code*

All three approaches are detailed below.


Dependencies
------------

pyCMBS was built to have only a small number of external dependencies. However, there are certainly dependencies existing which are a requirement for using pyCMBS sucessfully.

Mandatory dependencies
~~~~~~~~~~~~~~~~~~~~~~

- netCDF4 library is used as default for file I/O.
- (Nio library could be installed as an alternative)

Recommended dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

- cartopy is highly recommended for the generation of high quality map plots.


ZMAW users
----------

For users working at the ZMAW (KlimaCampus, Hamburg), the usage of pyCMBS is straightforward. Nothing has to be installed! You just need to login into your shell and load the pyCMBS module as follows ::

    module load pyCMBS/default

That's it. This sets also automatically the settings for the *PYTHONPATH* and *SEP* environment variables which are needed.

Please check if anything is working like described below_.



For the rest of the world
-------------------------

Installation using PIP
~~~~~~~~~~~~~~~~~~~~~~

A very easy way to install python packages is by using PIP. If you have installed PIP on your environment,  pyCMBS installation is as simple as::

    pip install pycmbs

Check if everything is working,  like described below_.


Installation from a *tarball*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have obtained the pyCMBS code from a *tarball*, then do the following::

    tar -xvf pycmbs-vx.x.x.tar.gz

This expands the tar file and you obtain a subdirectory,  called *pycmbs-vx.x.x*

To install the package,  use::

    cd pycmbs-vx.x.x
    python setup.py install

This will install the package in your python environment. Check successful installation,  like described below_.



Direct usage of code repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another alternative to work with pyCMBS is to used the master branch of the code repository. Currently, pyCMBS code repository version control is based on ```svn```. The repository access is granted by the author to anybody interested to contribute to pyCMBS.

The procedure to work with pyCMBS from development branch is as follows:

1. Create a directory *pyCMBS*::
   mkdir pyCMBS
   cd pyCMBS


2. Check out the code::
   svn checkout https://svn.zmaw.de/svn/pycmbs/trunk .

Be aware of the '.'

3. Set the environment
   pyCMBS is a python module that needs to be found by the python interpreter during runtime. Normally, modules are installed in the a subdirectory of your python installation. Typically, this is something like */usr/local/lib/python2.7/dist-packages*.
   As you have now installed pyCMBS in a different location, you need to tell python where to find it. That's what the *PYTHONPATH* variable is used for. The following steps will lead you to a working environment (examples are given for *bash* shell)::

    # get name of current directory
    # (this should be the directory where the pyCMBS subfolder is located)

    m300028@alexnotebook:pycmbs$ pwd
    /home/m300028/shared/dev/svn/pyCMBS-GIT-Repo/pycmbs

    # now, you specify a new environment variable pointing to the root directory
    export PYCMBSROOT=/home/m300028/shared/dev/svn/pyCMBS-GIT-Repo/pycmbs

    # and another, pointing to the actual installation
    export PYCMBSPATH=$PYCMBSROOT/pyCMBS

    # extending your PYTHONPATH will ensure that python finds the module
    export PYTHONPATH=$PYCMBSROOT:$PYTHONPATH

    # and finally, you extend your path, so that a script is automatically found by the OS
    export PATH=$PYCMBSPATH/framework:$PATH

It is recommended to put the above commands in your *.bashrc*, so they are automatically loaded.


4. make script *pycmbs.py* executable::

   cd $PYCMBSPATH/framework
   chmod u+w pycmbs.py




.. _below:

Final check of installation
---------------------------

Check that installation works as follows:

1. check if pyCMBS python module is loaded properly::
   python -c "from pyCMBS import *; print('Welcome to pyCMBS')"

This should give you a short Welcome message, but no error messages.

2. Check if the benchmarking script 'pycmbs.py' is found, by typing::
   pycmbs.py

This will you give a message like::

   *******************************************
   * WELCOME to pycmbs.py                    *
   * Happy benchmarking ...                  *
   *******************************************

and will end with an error message that the configuration file is not found (this is o.k.)

If you see the above, the installation has worked! Congratulations!

3. Check also the proper installation of the cdo's and the cdo.py interface, as this is a prerequesite of beeing able to properly work with pyCMBS::

     python -c "from cdo import *; cdo=Cdo(); print 'If you see this, everything went right ... have fun with pyCMBS and CDOs'"

Again, this should give you a short welcome message. Any error message is a bad sign. In that case, please check your installation again. Have a look at the troublesolver_.


.. _installation_details:

Some more details on what the installation does
------------------------------------------------

setting SEP environment variable
    xxxxxxxxxxx

.. _troublesolver:

Some hints for trouble solving
------------------------------

1. Is python working properly?::

    python -c "print 'Hello world'"

2. Can pyCMBS be found in your *PYTHONPATH*?::

    echo $PYTHONPATH

This should give you the path where python is searching for modules. If it is empty you are most likely in trouble. Check if you have a valid python installation.













