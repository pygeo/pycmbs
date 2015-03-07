.. _bench-testing:

Testing the benchmarking framework (experts only)
=================================================

To technically test the benchmarking framework an own framework for testing has been developed. This is based on nosetests and the code is located in the file `test_bench.py`. The framework requires test data to run. Both, model output and observations are required. Fore details see :ref:`testdata`

.. _testdata:

Getting testdata for testing
----------------------------

Getting model data
~~~~~~~~~~~~~~~~~~

Tests were implemented so far for the model classes CMIP5RAWSINGLE and CMIP5RAW. Both rely on cmorized output in CMIP5 formats. To get testdata, results for a single model are sufficient. To extract the test data, the following steps are needed:

1. extract test data from MiKlip server archive using the `miklip_preprocessor`. Retrieve the required fields as follows::

    python miklip_preprocessor.py -i /home/zmaw/m300028/cmip5_miklip -o /scratch/m300028 -e amip \\
                                  -v rsds rsus pr tas sfcWind evspsbl hfls prw huss -cfg dummycfg

This will extract the data as described in :ref:`benchmarking-cmip`.

2. Choose the MPI-M folder and copy it to the testdata directory which is specified in `test_becnh.py`. by default, this is currently the directory `pycmbs_root/testdata/miklip`. Your structure should then look something like::

    /testdata/miklip.cfg
             /MPI-M/
                   /MPI-ESM-LR/
                              /amip/...
                   /MPI-ESM-MR/...

3. get also observation data as specified in the INI configuration files (see :ref:`getting-obs-data`)

4. execute the tests using `nosetests <https://nose.readthedocs.org/en/latest/>_`::

    # change into test directory
    cd <pycmbsdir>/benchmarking/tests/test_benchmarking
    nosetests

This will execute tests, for both, CMIP5RAWSINGLE and CMIP5RAW classes. The variables tested are specified in `test_bench.py`

**Note: The tests are only of technical character. Thus they only indicate that data I/O and benchmarking environment is working in principle well. A proper testing is nevertheless required to ensure that the results are valid (e.g. data scaling, valid areas ...)**


.. _getting-obs-data:

Getting observation data
~~~~~~~~~~~~~~~~~~~~~~~~

Observational data can come from various sources which depend on the input datasets specified in the INI configuration files. In general, the user can store the observational data in two different ways:

1. store data at custom places: provide a custom filename in the INI file
2. on a central data pool: use routine `get_data_pool_directory()` to specify the file location

The latter has the advantage, that when working on different file systems, the filenames don't need to be adapted, but only the name of the root directory where the data is stored on different servers. For examples see the code in the INI files provided.
