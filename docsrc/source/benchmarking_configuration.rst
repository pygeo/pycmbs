Benchmarking configuration (o.k.)
=================================

Plot configuration details
--------------------------

For each variable, a configuration file with the extension *.ini* is used to specify. The INI files are expected to be located in a directory which is specified in the main configuration file (.cfg).
The plot configuration file specifies for each variable

 * which diagnostics should be applied to a certain variable
 * how plots for a particular diagnostic should look like (e.g. colorbars,
   limits ...)
 * which observational datasets should be used

An INI file has two major parts:

 1. Global plot options
 2. Observation specific plot options (for each used observational dataset)

Global plot options
~~~~~~~~~~~~~~~~~~~

The global plot options have the following structure (example below)::

    [OPTIONS]
    map_difference =  True
    map_seasons    =  True
    map_season_difference = False
    reichler_plot  =  True
    gleckler_plot   =  True
    hovmoeller_plot   =  False
    regional_analysis = True
    global_mean    = True
    vmin           =  0.
    vmax           =  8.
    dmin           =  -1.
    dmax           =  1.
    units          =  $mm/day$
    label          =  Daily evaporation
    cticks         = [0.,2.,4.,6.,8.,10.]
    nclasses       = 8
    preprocess     = True
    interpolation  = conservative
    targetgrid     = t63grid
    projection     = robin
    region_file    = /home/m300028/shared/data/CMIP5/evap/evspsbl/merged/dummy_mask2.nc
    region_file_varname = regmask


*map_difference* [True,False]
    use diagnostic to plot difference between models and observations

*map_seasons* [True,False]
    use diagnostic to plot climatological monthly mean or seasonal mean maps of models and observations

*map_season_difference* [True,False]
    same as map_seasons, but for difference between models and observations.

*reichler_plot* [True,False]
    Summarize error skill score for this variable at the end of the section for this variable.

*gleckler_plot* [True,False]
    Use this variable in the *Portraet Diagram* at the end of the report.

*hovmoeller_plot* [True,False]
    Generate a hovmoeller plot for the variable, for both observations and models.

*regional_analysis* [True,False]
    Perform regional analysis (statistics and correlation) of observations and models per variable.

*region_file*
    Name of netCDF file which contains the rasterized region IDs (user needs to ensure that the same geometry as the target grid is provided)

*region_file_varname*
    name of variable in *region_file*, which shall be read to identify regions; note that the data is interpreted as integer values.

*global_mean*
    generate a global mean plot for this variable (see XXXX)

*vmin*
    minimum plotting limit for data

*vmax*
    maximum plotting limit for data

*dmin*
    minimum plotting limit for difference plot

*dmax*
    maximum plotting limit for difference plot

*units*
    string to specify units of the variable. This is used for automatic labeling of plots. Note that all text can be used which can also be used for labelling in matplotlib. In particular the usage of $ is usefull to render text using latex (e.g. $\frac{a}{b}$ will plot you the a/b in a nice way).

*label*
    label text to be used for the variable

*cticks*
    tick labels for colormap


Observation specific plot options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below the global options, one can include an arbitrary number of observations. 

Each observation is specified by a block of configuration parameters, like in the following example.::

    [CLARASAL]
    obs_file =  #get_data_pool_directory() + 'data_sources/CMSAF/CLARA-SAL/DATA/SAL_all_t63.nc'#
    obs_var  =  sal
    scale_data = 0.01
    gleckler_position = 2
    add_to_report = True
    valid_mask = land

*[Observation_Identifier]* : str
    unique identified for the observation. Will be used e.g. in plots as labels

*obs_file* : str
    name of observation file. Here the user can either specify a full path name
    to a file or, like shown in the example above, execute a python command
    that is used to construct the filename. In the above example, the hash (#)
    is used to identify a python command. If the value of obs_file starts and
    ends with a hash, then the string in between is executed like you would
    execute a python command. Here, the routine get_data_pool_directory() is
    called, which returns a path name and then the remaining path to the
    observational data file is appended.

*obs_var* : str
    name of variable in observation file

*scale_data* : float
    scaling factor to be applied on data of the file. This is e.g. usefull if
    the netCDF file does not contain an own scale_factor attribute or if you
    want to apply simple conversions (e.g. from kg/m**2 s to mm/day for
    precipitation). The data is multiplied by the scaling factor.

*gleckler_position* : int
    [1,2,3,4] position of the observational dataset in the Portraet diagram. Up
    to four different datasets can be shown at once. The meaning if the numbers
    is as follows: 1=top, 2=bottom, 3=left, 4=right

*add_to_report* : str
    add this observation to the report [True,False]

*valid_mask* : str
    [land,ocean,global]; specifies if a mask shall be applied to the dataset
    and model. If 'land', then all ocean areas are masked if 'ocean', then all
    land areas are masked. For any other options, the whole globe is used.
