Your first benchmarking
-----------------------

The general workflow for using pyCMBS for benchmarking is shown in Figure ????
(TBD).
The main input is model data and observational datasets as well as a user specified configuration. After running the benchmarking, you get a report (PDF) as well as a lot of figures and statistics which contain usefull information.

The next steps will guide you through a benchmarking session to get you started.

1. set up a working directory and change to it::

    # set up a working directory somewhere on your machine
    mkdir -p ~/temp/my_first_benchmarking
    cd ~/temp/my_first_benchmarking

2. set up an initial configuration::

    # simply run the pycmbs.py script
    # if properly installed, you should be able to just run it from the console
    pycmbs.py init

This gives you::

    $ ls
    configuration  pyCMBS_template.cfg
    
If you do this for the first time,  then it is recommended that you make yourself familiar with the content of the *configuration* directory. This contains

* INI files for each variable which specify the plot and processing configuration
* *.json files which specify interface routines for data conversion

3. adapt the configuration file (*.cfg) to your needs
   The configuration file is the center part where you specify
   
   a) which variables shall be diagnosed
   b) which models shall be analysed
   
   Details about the configuration file sre specified here_.
   
4. Do it!
   Run the benchmarking now by executing::
   
       pycmbs.py your_config_script.cfg
       
       
--------

.. _here::

Structure of configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The configuration file has the following structure:

*Global options*::

    xxxxxxxxxxx






*Model options*::

    xxxxxxxxxxxxx












