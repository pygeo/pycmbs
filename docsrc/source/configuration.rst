Plot configuration
------------------

For each variable, a configuration file with the extension *.ini* is used to specify. The INI files are expected to be located in a directory which is specified in the main configuration file (??? cross link here ???)

 * which diagnostics should be applied to a certain variable
 * how plots for a particular diagnostic should look like
 * which observational datasets should be used

An INI file has two major parts:

 1. Global plot options
 2. Observation specific plot options

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
    use diagnostic to plot difference between models and observations (see XXXXX).

*map_seasons* [True,False]
    use diagnostic to plot climatological monthly mean or seasonal mean maps of models and observations (see XXXX).

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
    xxxxxxxxxxxxxxxxx


Observation specific plot options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TBD





