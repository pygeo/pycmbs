#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
'''


#--- Analysis compliance matrix (which was already checked)
# a P indicates that the preprocessing to seasonal data is performed automatically
# and works

#                    | CMIP5 | JSBACH_BOT | JSBACH RAW | report | Remark
#albedo analysis     |   P   |            |      X     |    x   | preprocessing includes calculation of albedo from SW up and downward fluxes as well as regridding to T63 grid if needed
#SIS analysis        |   P   |            |      X     |    x   |
#rainfall analysis   |   P   |      P     |      X     |    x   |
#temperature         |   P   |            |            |    x   |
#veg. fraction       |       |            |            |        |
#phenology           |       |            |            |        |  external framework
#snow fraction       |       |            |            |        |



# - regional analysis based on an input mask
#@todo: correlation and RMSE analysis and Taylor plotting
#
# other colorbar for vegetation fraction analysis
# implement grass cover fraction analysis

#
# todo: interpolation method for fractional coverage !!!




__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "0000/00/00"

#============ IMPORTS ==================================================

from pyCMBS import *

#--- always use plot backend which is not interactive for benchmarking framework
import matplotlib
matplotlib.use('agg')

import matplotlib.pylab as pl
import sys
import os

#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html

#--- framework specific modules ---
from models   import *
from config   import *
from analysis import *

#=======================================================================



#=======================================================================

def get_analysis_scripts():
    """
    returns names of analysis scripts for all variables as a dictionary
    in general, these names can be also read from an ASCII file
    """
    d={}
    d.update({'rain':'rainfall_analysis'})
    d.update({'albedo':'albedo_analysis'})
    d.update({'sis':'sis_analysis'})
    d.update({'surface_upward_flux':'surface_upward_flux_analysis'})
    d.update({'tree':'tree_fraction_analysis'})
    d.update({'grass':'grass_fraction_analysis'})
    d.update({'phenology_faPAR':'phenology_faPAR_analysis'})
    d.update({'temperature':'temperature_analysis'})
    d.update({'evap':'evaporation_analysis'})
    d.update({'wind':'wind_analysis'})
    d.update({'twpa':'twpa_analysis'})
    d.update({'wvpa':'wvpa_analysis'})
    d.update({'hair':'hair_analysis'})
    d.update({'late':'late_analysis'})
    d.update({'budg':'budg_analysis'})

    return d



#=======================================================================



#####################################################################
#####################################################################
#####################################################################

pl.close('all')

#-----------------------------------------------------------------------

def get_methods4variables(variables, model_dict):
    """
    for a given list of variables, return a dictionary
    with information on methods how to read the data

    IMPORTANT: all options provided to the routines need to be
    specified here and arguments must be set in calling routine get_data()
    """

    hlp={}
    hlp.update({'rain' : 'get_rainfall_data(interval=interval)'})
    hlp.update({'albedo' : 'get_albedo_data(interval=interval)'})
    hlp.update({'sis' : 'get_surface_shortwave_radiation_down(interval=interval)'})
    hlp.update({'surface_upward_flux' : 'get_surface_shortwave_radiation_up(interval=interval)'})
    hlp.update({'tree' : 'get_tree_fraction(interval=interval)'})
    hlp.update({'grass' : 'get_grass_fraction(interval=interval)'})
    hlp.update({'phenology_faPAR' : 'get_faPAR(interval=interval)'})
    hlp.update({'temperature' : 'get_temperature_2m(interval=interval)'})
    hlp.update({'snow' : 'get_snow_fraction(interval=interval)'})
    hlp.update({'evap': 'get_model_data_generic(interval=interval, **%s)' % model_dict['evap']})
    hlp.update({'wind': 'get_model_data_generic(interval=interval, **%s)' % model_dict['wind']})
    hlp.update({'twpa': 'get_model_data_generic(interval=interval, **%s)' % model_dict['twpa']})
    hlp.update({'wvpa': 'get_model_data_generic(interval=interval, **%s)' % model_dict['wvpa']})
    hlp.update({'late': 'get_model_data_generic(interval=interval, **%s)' % model_dict['late']})
    hlp.update({'budg': 'get_model_data_generic(interval=interval, **%s)' % model_dict['budg']})
    hlp.update({'hair': 'get_model_data_generic(interval=interval, **%s)' % model_dict['hair']})


    res={}
    for k in hlp.keys(): #only use the variables that should be analyzed!
        if k in variables:
            res.update({k:hlp[k]})

    #--- implement here also dependencies between variables for anylssi
    #e.g. phenology needs faPAR and snow cover fraction. Ensure here that
    # snow cover is also read, even if only phenology option is set
    if ('phenology_faPAR' in variables) and not ('snow' in variables):
        res.update({'snow' : 'get_snow_fraction()'})

    return res




#=======================================================================
#=======================================================================

"""
HOWTO

Add a new variable:
1) register variable in get_methods4variables()
2) implement for each data object a routine how to read the data
3) implement an analysis script that performs the actual analysis
4) register this analysis script in get_analysis_scripts()
"""



#/// check commandline options ///
if len(sys.argv) > 1:
    if len(sys.argv) == 2:
        file = sys.argv[1] #name of config file
        if not os.path.exists(file):
            raise ValueError, 'Configuration file can not be found: ' + file
    else:
        raise ValueError, 'Currently not more than one command line parameter supported!'
else: #default
    file='pyCMBS.cfg'

#/// read configuration file ///
CF = ConfigFile(file)

#/// read plotting options ///
PCFG = PlotOptions()
PCFG.read(CF)

plot_options=PCFG





for thevar in plot_options.options.keys():


    if thevar in plot_options.options.keys():
        print 'Variable: ', thevar
        for k in plot_options.options[thevar].keys():

            print '    Observation: ', k




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



# get observations?
model_dict = {'rain':  {'variable': 'pr',
                        'unit': 'mm/day',
                        #'interval': 'monthly',
                        'lat_name': 'lat',
                        'lon_name': 'lon',
                        'model_suffix': 'ensmean',
                        'model_prefix': 'Amon',
                        'file_format' : 'nc',
                        'scale_factor': 86400.,
                        'mask_area': 'ocean'},

              'evap':   {'variable': 'evspsbl',
                        'unit': 'mm/day',
                        #'interval': 'monthly',
                        'lat_name': 'lat',
                        'lon_name': 'lon',
                        'model_suffix': 'ensmean',
                        'file_format' : 'nc',
                        'model_prefix': 'Amon',
                        'scale_factor': 86400.,
                        'mask_area': 'ocean'},

              'twpa':   { 'variable': 'clwvi',
                         'unit': 'kg/m^2',
                         #'interval': 'monthly',
                         'lat_name': 'lat',
                         'lon_name': 'lon',
                         'model_suffix': 'ensmean',
                         'file_format' : 'nc',
                         'model_prefix': 'Amon',
                         'scale_factor': 1.,
                         'mask_area': 'ocean'},

             'wind':    {'variable': 'sfcWind',
                         'unit': 'm/s',
                         #'interval': 'monthly',
                         'lat_name': 'lat',
                         'lon_name': 'lon',
                         'model_suffix': 'ensmean',
                         'file_format' : 'nc',
                         'model_prefix': 'Amon',
                         'scale_factor': 1.,
                         'mask_area': 'ocean'},

              'wvpa':   {'variable': 'prw',
                        'unit': 'kg m^2',
                        #'interval': 'monthly',
                        'lat_name': 'lat',
                        'lon_name': 'lon',
                        'model_suffix': 'ensmean',
                        'file_format' : 'nc',
                        'model_prefix': 'Amon',
                        'scale_factor': 1,
                        'mask_area': 'ocean'},

              'late':   {'variable': 'hfls',
                        'unit': 'W/m^2',
                        #'interval': 'monthly',
                        'lat_name': 'lat',
                        'lon_name': 'lon',
                        'model_suffix': 'ensmean',
                        'file_format' : 'nc',
                        'model_prefix': 'Amon',
                        'scale_factor': 1,
                        'mask_area': 'ocean'},

              'hair':   {'variable': 'huss',
                         'unit': '$kg/kg^2$',
                         #'interval': 'monthly',
                         'lat_name': 'lat',
                         'lon_name': 'lon',
                         'model_suffix': 'ensmean',
                         'file_format' : 'nc',
                         'model_prefix': 'Amon',
                         'scale_factor': 1,
                         'mask_area': 'ocean'},

              'budg':   {'variable': 'budg',
                        'unit': 'mm/d',
                        #'interval': 'monthly',
                        'lat_name': 'lat',
                        'lon_name': 'lon',
                        'model_suffix': 'ensmean',
                        'file_format' : 'nc',
                        'model_prefix': 'Amon',
                        'scale_factor': 86400.,
                        'mask_area': 'ocean',
                        'custom_path' : '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data'}}


#--- names of analysis scripts for all variables ---
scripts = get_analysis_scripts()

#/// get dictionary with methods how to read data for variables to be analyzed ///
variables = CF.variables
varmethods = get_methods4variables(variables, model_dict)


#=======================================================================
#global_settings_dict = {'landsea_mask': {'filename': ''}}
#default_analysis = {}



#/// READ DATA ///
'''
create a Model instance for each model specified
in the configuration file

read the data for all variables and return a list
of Data objects for further processing
'''

#mean_model = Data(None,None)
model_cnt   = 1; proc_models = []

for i in range(len(CF.models)):
    #--- assign model information from configuration ---
    data_dir   = CF.dirs[i]
    model      = CF.models[i]
    experiment = CF.experiments[i]

    #--- create model object and read data ---
    # results are stored in individual variables namex modelXXXXX
    if CF.dtypes[i].upper() == 'CMIP5':
        themodel = CMIP5Data(data_dir,model,experiment,varmethods,intervals=CF.intervals,lat_name='lat',lon_name='lon',label=model,start_time=start_time,stop_time=stop_time,shift_lon=shift_lon)
    elif CF.dtypes[i].upper() == 'JSBACH_BOT':
        themodel = JSBACH_BOT(data_dir,varmethods,experiment,intervals=CF.intervals,start_time=start_time,stop_time=stop_time,name=model,shift_lon=shift_lon)
    elif CF.dtypes[i].upper() == 'JSBACH_RAW':
        themodel = JSBACH_RAW(data_dir,varmethods,experiment,intervals=CF.intervals,start_time=start_time,stop_time=stop_time,name=model,shift_lon=shift_lon)
    else:
        print CF.dtypes[i]
        raise ValueError, 'Invalid model type!'

    #--- read data for current model ---
    themodel.get_data()

    #--- copy current model to a variable named modelXXXX ---
    cmd = 'model' + str(model_cnt).zfill(4) + ' = ' + 'themodel.copy(); del themodel'
    exec(cmd) #store copy of cmip5 model in separate variable

    #--- append model to list of models ---
    proc_models.append('model' + str(model_cnt).zfill(4))
    model_cnt += 1


#########################################################
# MULTIMODEL MEAN
#########################################################
#--- here we have now all the model and variables read. The list of all models is contained in the variable proc_models.
# CALCULATE MULTIMODEL MEAN

if True:

    #raise ValueError, 'This mean model appraoch can not work, as the timesteps are not the same!!! We need to use the mean climatologigy! --> preprocessing (generic analysis ???)'

    for i in range(len(proc_models)):
        exec('actmodel = ' + proc_models[i] + '.copy()')

        if i == 0:
            model_mean = actmodel.copy()
            model_mean.name = 'mean-model'
            model_mean._unique_name = 'model_mean'
        else:
            for k in model_mean.variables.keys():

                print 'Processing ... ', k, proc_models[i]

                #the variables[] list contains Data objects!
                hlp1 = model_mean.variables[k] #is a Data object or a tuple!
                hlp2 = actmodel.variables[k]

                if isinstance(hlp1,tuple):
                    #the mean makes only sense for climatological mean values. Anything else should be not supported due to possbily different timestamps
                    #theD = hlp1[2].copy()

                    #theD[0][:] = None; theD[1][:] = None; theD.add(hlp2[2],copy=False)
                    #theD.label='Mean-model'
                    model_mean.variables.update( { k : (None,None,None) } )
                else:
                    theD = hlp1.copy()
                    theD.add(hlp2,copy=False) #SUM: by using masked arrays, the resulting field is automatically only valid, when both datasets contain valid information!
                    theD.label = 'Mean-model'
                    model_mean.variables.update( { k : theD } )
                del hlp1,hlp2 #, theD
        del actmodel

    #now we have the sum and can calculate the average
    for k in model_mean.variables.keys():
        hlp1 = model_mean.variables[k]
        if isinstance(hlp1,tuple):
            pass
            #model_mean.variables.update( { k : (hlp1[0],hlp1[1],hlp1[2].mulc(1./float(len(proc_models)),copy=False )  ) } )
        else:
            model_mean.variables[k].mulc(1./float(len(proc_models)),copy=False) #weight with number of models

    #add to list of models to process
    proc_models.append('model_mean')


#########################################################
# END MULTIMODEL MEAN
#########################################################




#/// prepare global benchmarking metrices
#generate a global variable for Gleckler plot!
global_gleckler = GlecklerPlot()

########################################################################
# REPORT
########################################################################
rep = Report(CF.options['report'],'pyCMBS report - ' + CF.options['report'],CF.options['author'],outdir='./report_' + CF.options['report'] + '/',dpi=300,format='pdf')
cmd = 'cp ../logo/Phytonlogo5.pdf ' + rep.outdir
os.system(cmd )



########################################################################
# MAIN ANALYSIS LOOP: perform analysis for each model and variable
########################################################################
skeys = scripts.keys()
for variable in variables:

    #/// register current variable in Gleckler Plot
    global_gleckler.add_variable(variable)

    #/// call analysis scripts for each variable
    for k in range(len(skeys)):
        if variable == skeys[k]:

            print 'Doing analysis for variable ... ', variable
            print '   ... ', scripts[variable]
            model_list = str(proc_models).replace("'","")  #... model list is reformatted so it can be evaluated properly
            cmd = scripts[variable]+'(' + model_list + ',GP=global_gleckler,shift_lon=shift_lon,use_basemap=use_basemap,report=rep,interval=CF.intervals[variable],plot_options=PCFG)'
            eval(cmd) #run analysis

#/// generate Gleckler analysis plot for all variables and models analyzed ///
global_gleckler.plot(vmin=-0.8,vmax=0.8,nclasses=25)


rep.section('Summary error statistics')
rep.figure(global_gleckler.fig,caption='Gleckler et al. (2008) model preformance index')

#legend for gleckler plot
for variable in variables:
    if variable not in PCFG.options.keys(): #check if variable existing
        continue
    varoptions = PCFG.options[variable]
    thelabels={}
    for k in varoptions.keys(): #keys of observational datasets
        if k == 'OPTIONS':
            continue
        else:
            thelabels.update({int(varoptions[k]['gleckler_position']) : k}) #generate dictionary for GlecklerPLot legend
    fl = global_gleckler._draw_legend(thelabels,title=variable.upper())
    rep.figure(fl,width='8cm',bbox_inches=None)
    del fl, thelabels

#---



#/// close report ///
rep.close()

pl.close('all')


#~ if __name__ == '__main__':
    #~ main()


