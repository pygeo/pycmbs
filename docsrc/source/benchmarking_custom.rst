===========================================
Customizing pyCMBS benchmarking environment
===========================================

The model benchmarking framework can be easily customized and adapted to the user needs. IN the following, we will cover the following topics:

# How to add new model data for an already existing model output format?
# How to add new observational datasets?
# How to integrate new variables in the analysis?
# How to use already existing external scripts together with the pyCMBS environment?
# How to add a new model format?

How to add new model data for an already existing model output format?
----------------------------------------------------------------------

TBD






How to add new observational datasets?
--------------------------------------
The integration of new observational datasets is very simple as long as the datasets you use follow some standard conventions:

* datasets are in netCDF format
* optional: datasets have metadata attributes in the netCDF file. pyCMBS is making automatically use of CF conventions like netCDF attributes like *scale_factor*, *add_offset*, *_FillValue*, *units*. In case these attributes are provided, they are automatically used
* lat/lon coordinates are provided either as vectors (for simple lat/lon projections) or as 2D fields (a (lat, lon) tuple for each grid cell).
* observations are stored in a single file (all timesteps included)

These are the very basic requirements that pyCMBS can make use of new observations. 

Steps to integrate a new observational dataset into pyCMBS are as follows:

# decide for the variable the observational dataset belongs to --> variable name; you can look in the configuration file (.cfg) to get the currently supported variable names
# modify the corresponding INI file

 * let's say, that you have chosen *sis* (surface solar irradiance) as the variable and you have a new surface radiation dataset. Then the corresponding INI file would be *sis.ini*. The INI files can be found in the *framework/configuration* folder. 
 * You can however also generate an own, new configuration folder, by simply typing *pycmbs.py init* in a fresh directory

 * The content of the INI file is self explanatory. You have a global section which specifies how the analysis for this particular variable shall be made (e.g. which diagnostics and plots shall be generated). Below, you have for each observational dataset a section which specifies the details for each observation. Such a section looks e.g. like the following

```
[CERES]
obs_file = #get_data_pool_directory() + 'data_sources/CERES/EBAF/ED_26r_SFC/DATA/CERES_EBAF-Surface__Ed2.6r__sfc_sw_down_all_mon__1x1__200003-201002.nc'#
obs_var  = sfc_sw_down_all_mon
scale_data = 1.
gleckler_position = 2
add_to_report = True
valid_mask = global
```

The different entries have the following meaning

obs_file : path to the netCDF file with the observations
    here you can either specify an *absolute* path or you can out a python code snipped between two hashes '#'. The latter approach is usefull, when you have some function that directs to a particular directory, like in the example given. Otherwise, the easiest way is to just put the *absolute path* to the netCDF file.

obs_var : specifies the name of the variable in the netCDF file

scale_data : a multiplicative scaling factor that is applied to the data when it is read

gleckler_position : [1...4]; This number specifies, where the observational dataset will be placed in the portraet diagram (???)

add_to_report : [True,False], specifies if the observational dataset should be included in the report or not. This allows to have a lot of configurations in the INI file, but use only a few of them.

valid_mask : [land,ocean,global], specifies which area(s) are supposed to contain valid data. Other regions are automatically masked. Thus if you specify e.g. *land*, then the ocean will be masked. It is important, that you use a 

** Adding a new observation is as simple as copy/paste an already existing section and modify the entries like you need it. That's it!**



Recepies for handling problems:

In 80% of the cases, pyCMBS will handle your new data smoothly. However, it might happen that your file(s) are different from the files pyCMBS was tested so far with. For these cases the following steps might help to solve your problem:

# Is the file o.k?
 
 * Have a look at the file with other tools like e.g. ncview or panoply
 * make also an "ncdump -h" to check the metadata of the file

# Can *cdo's* work with the file?: The preprocessing capabilities of pyCMBS largely rely on the usage of the climate data operators (cdo). If the *cdo's* can not work with your file, then pyCMBS will most likely have also problems.

 * check if *cdo's* can in general read the file: **cdo sinfo <filename>**
 * check if grid of the file is recognized by trying to remap the file manually using **cdo remapcon,t63grid <infile> nothing.nc**

If one of the two tests above fail, then your file is missing some essential metadata or has a strange grid or grid description that is not automatically recognized. In these cases, it would be best, if you try to figure out, why the *cdo's* are not capable to work with your dataset. Try to pose your question to the *cdo's* help forum (don't forget to provide details about your file; e.g. by sending the results of *ncdump -h*)




How to integrate new variables in the analysis?
-----------------------------------------------

To add new variables in pyCMBS implies the following steps:

1. **Define I/O routine:** Implement for each model class that shall support the new variable a routine
that allows to read the data. Let's say you have a variable *sis*, then you
would need e.g. to implement a routine *get_sis()* for the CMIP5 model class.
Note that there is already a routine which can be used for generic I/O.

2. **Register I/O routine**: After you have implemented the routine to read the
data, you need to let the program know about it. All data is read using a
routine called *get_data()*. This routine gets the information which
subroutines to call from details provided in a configuration file. The
configuration file is found in::
    ./configuration/model_data_routines.json

The file is a simple JSON dictionary. Make yourself a bit familar with the
structure and it should not be a problem to implement your new routine there.

3. **Analysis script:** Now you have the analysis script that can be used to
read the data. However, you still need to tell pyCMBS how to make use of this
new information. This you do by implementing an analysis routine in
*analysis.py*. For most variables supported so far, this analysis routine is
just a wrapper which is calling a very generic analysis routine that basically
does everything you tell it to do. What to do is specified in the INI files for
each variable. Note however, that you are free to do what you want and you can
implement a new analysis routine which is doing right the thing you want it to
do.

4. **Last step** is to tell pyCMBS that the analysis script you implemented is
existing. This is again done, by simply registering it in the following file::
    ./configuration/analysis_scripts.json


How to use external scripts?
----------------------------
TBD


How to add a new model format?
------------------------------

TBD



