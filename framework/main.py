#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
pyCMBS - climate model bechmarking framework
'''


#--- Analysis compliance matrix (which was already checked?)

#                    | CMIP5 | JSBACH_BOT | JSBACH RAW |
#albedo analysis     |    x  |            |            |
#SIS analysis        |    x  |            |            |
#rainfall analysis   |       |      X     |            |
#temperature         |       |            |            |
#veg. fraction       |       |            |            |
# phenology          |       |            |            |

#preprocessing of observations and/or data ???

#land sea mask for model in SIS analysis


#@todo: implement JSBACH raw data
#@todo: implement ATM/BOT files

#@todo: systematic validation of zonal mean statistics using som reference cases

#todo
#
# @todo: implement temperature analysis
# @todo: implement precipitation analysis
#                 implement more observational datasets!!!


#
# @todo: implement temporary directory for pyCMBS processing


# TODO CMIP5:
# - seldate appropriately
# - check timestamp!
# seasmean calculations automatically!!

#todo check datetime; something is wrong! as data starts in Dcember 1978!


# area weigting of zonal means and also area weighting of RMS errors etc.


# - regional analysis based on an input mask
#@todo: correlation and RMSE analysis and Taylor plotting
#
# other colorbar for vegetation fraction analysis
# implement grass cover fraction analysis

#
# todo: interpolation method for fractional coverage !!!
# area weighting for correlation analysis

# - significance tests of differences
# - pre-processing scripts from external configuration file
# - regional subsetting options ??



__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "0000/00/00"

#============ IMPORTS ==================================================

from pyCMBS import *

import matplotlib.pylab as pl

#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes

#--- framework specific modules ---
from models   import *
from config   import *
from analysis import *

#=======================================================================



#=======================================================================

def get_analysis_scripts():
    '''
    returns names of analysis scripts for all variables as a dictionary
    in general, these names can be also read from an ASCII file
    '''
    d={}
    d.update({'rain':'rainfall_analysis'})
    d.update({'albedo':'albedo_analysis'})
    d.update({'sis':'sis_analysis'})
    d.update({'tree':'tree_fraction_analysis'})
    d.update({'grass':'grass_fraction_analysis'})

    return d



#=======================================================================



#####################################################################
#####################################################################
#####################################################################

from pyCMBS import *
pl.close('all')

#-----------------------------------------------------------------------

def get_methods4variables(variables):
    '''
    for a given list of variables, return a dictionary
    with information on methods how to read the data
    '''

    hlp={};
    hlp.update({'rain' : 'get_rainfall_data()'})
    hlp.update({'albedo' : 'get_albedo_data()'})
    hlp.update({'sis' : 'get_surface_shortwave_radiation_down()'})
    hlp.update({'tree' : 'get_tree_fraction()'})
    hlp.update({'grass' : 'get_grass_fraction()'})

    res={}
    for k in hlp.keys(): #only use the variables that should be analyzed!
        if k in variables:
            res.update({k:hlp[k]})

    return res




#=======================================================================
#=======================================================================

'''
adding a new variable:
1) register variable in get_methods4variables()
2) implement for each data object a routine how to read the data
3) implement an analysis script that performs the actual analysis
4) regist this analysis script in get_analysis_scripts()
'''



#/// read configuration file ///
file='pyCMBS.cfg'
CF = ConfigFile(file)

if CF.options['basemap']:
    f_fast = False
else:
    f_fast=True
shift_lon = use_basemap = not f_fast

print 'Using Basemap: ', use_basemap

s_start_time = CF.start_date #'1983-01-01' #todo where is this used ?
s_stop_time  = CF.stop_date  #'2005-12-31'
start_time = pl.num2date(pl.datestr2num(s_start_time))
stop_time  = pl.num2date(pl.datestr2num(s_stop_time ))

#--- names of analysis scripts for all variables ---
scripts = get_analysis_scripts()

#/// get dictionary with methods how to read data for variables to be analyzed ///
variables = CF.variables
varmethods = get_methods4variables(variables)

#=======================================================================

#/// READ DATA ///
'''
create a Model instance for each model specified
in the configuration file

read the data for all variables and return a list
of Data objects for further processing
'''
model_cnt   = 1; proc_models = []
for i in range(len(CF.models)):
    #--- assign model information from configuration ---
    data_dir   = CF.dirs[i]
    model      = CF.models[i]
    experiment = CF.experiments[i]

    #--- create model object and read data ---
    # results are stored in individual variables namex modelXXXXX
    if CF.dtypes[i].upper() == 'CMIP5':
        themodel = CMIP5Data(data_dir,model,experiment,varmethods,lat_name='lat',lon_name='lon',label=model,start_time=start_time,stop_time=stop_time)
    elif CF.dtypes[i].upper() == 'JSBACH_BOT':
        themodel = JSBACH_BOT(data_dir,varmethods,experiment,start_time=start_time,stop_time=stop_time,name=model,shift_lon=shift_lon)
    else:
        print CF.dtypes[i]
        raise ValueError, 'Invalid model type!'

    #--- read data for current model ---
    themodel.get_data()
    cmd = 'model' + str(model_cnt).zfill(4) + ' = ' + 'themodel.copy(); del themodel'
    exec(cmd) #store copy of cmip5 model in separate variable

    #--- append model to list of models ---
    proc_models.append('model' + str(model_cnt).zfill(4))
    model_cnt += 1


#/// prepare global becnhmarking metrices
#generate a global variable for glecker plot!
global_glecker = GleckerPlot()


########################################################################
# MAIN ANALYSIS LOOP: perform analysis for each model and variable
########################################################################
skeys = scripts.keys()
for variable in variables:

    #/// register current variable in Glecker Plot
    global_glecker.add_variable(variable)

    #/// call analysis scripts for each variable
    for k in range(len(skeys)):
        if variable == skeys[k]:

            print 'Doing analysis for variable ... ', variable
            print '   ... ', scripts[variable]
            model_list = str(proc_models).replace("'","")  #model list is reformatted so it can be evaluated properly
            eval(scripts[variable]+'(' + model_list + ',GP=global_glecker,shift_lon=shift_lon,use_basemap=use_basemap)') #run analysis

#/// generate Glecker analysis plot for all variables and models analyzed ///
global_glecker.plot(vmin=-0.8,vmax=0.8,nclasses=15)


#~ if __name__ == '__main__':
    #~ main()

