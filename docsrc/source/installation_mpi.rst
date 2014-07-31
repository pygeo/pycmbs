Detailed installation guide for MPI-M users
===========================================

MPI-M users
-----------

For users working at the Max-Planck-Institute for Meteorology, the usage of pyCMBS is straightforward.
Nothing has to be installed! You just need to login into your shell and load the appropriate python module as follows::

    module avail
    # check which python module is available with latest version
    module load python/python-ve03   # or something similar

That's it. This sets also automatically the settings for the *PYTHONPATH* and *SEP* environment variables which are needed.

Please check if anything is working like described below_.


For advanced MPI-M users
------------------------

The pyCMBS distribution provided by CIS is associated always to a particular version.

If you want to work more with the code from the directly from the github repository, then do the following::

    module load python/python-ve03   # or something similar (will ensure that all dependencies are there)

    # decide for some directory to install the code
    export pycmbs_install_dir=$HOME/dev/pycmbs
    cd $pycmbs_install_dir

    # clone data from github
    git clone git@github.com:pygeo/pycmbs.git .

    # compile the necessary modules
    cd pycmbs
    python setup_extensions.py build_ext --inplace

    # update your PYTHONPATH
    export PYTHONPATH=$pycmbs_install_dir/pycmbs:$PYTHONPATH  # probably put this into your .bashrc

That's it. Please check if anything is working like described below_.

If you now want to work with the repository just update your local repository using the normal git procedure::

    git pull origin master

After that you should always compile ensure to recompile the cython modules::

    python setup_extensions.py build_ext --inplace








