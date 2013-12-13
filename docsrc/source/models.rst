Model interfaces in pyCMBS
--------------------------

Reading model data is implemented through dedicated classes which implement the logic how to read different model data.

The following model types are supported so far:

CMIP5:
    CMIP5 climate model output

JSBACH RAW:
    JSBACH output files. JSBACH is the land component of the MPI-ESM

Each class implements routines to read data for different variables. The central routine is *get_data()*. This is responsible to call all subroutines that read individual variables using functions the user can specify.

Details on how functions for reading new variables can be implemented are described in detail *here*. The following table provides an overview about the currently supported variables for the implemented model classes.

+--------------------------------+-------+------------+
|Variable                        | CMIP5 | JSBACH RAW |
+================================+=======+============+
|Surface radiation               |       |            |
+--------------------------------+-------+------------+
|Albedo                          |       |            |
+--------------------------------+-------+------------+
|Surface solar downwelling flux  |       |            |
+--------------------------------+-------+------------+
|Surface solar upwelling flux    |       |            |
+--------------------------------+-------+------------+

For a detailed description of the model individual interfaces,  please refer to the inline code documentation.
