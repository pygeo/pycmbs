#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2011-2013 Alexander Loew, alexander.loew@zmaw.de
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
import matplotlib as mpl
mpl.use('agg')

#import mpl.pylab as pl
pl = mpl.pylab

import sys
import os

#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html

#--- framework specific modules ---
from models   import *
from config   import *
from analysis import *
import pickle

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
    d.update({'albedo_vis':'albedo_analysis_vis'})
    d.update({'albedo_nir':'albedo_analysis_nir'})
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
    d.update({'seaice_extent':'seaice_extent_analysis'})
    d.update({'seaice_concentration':'seaice_concentration_analysis'})
    d.update({'gpp':'gpp_analysis'})
    d.update({'cfc':'cfc_analysis'})


    return d


def get_methods4variables(variables, model_dict):
    """
    for a given list of variables, return a dictionary
    with information on methods how to read the data

    IMPORTANT: all options provided to the routines need to be
    specified here and arguments must be set in calling routine get_data()
    """

    hlp={}
    #hlp.update({'rain' : 'get_rainfall_data(interval=interval)'})
    hlp.update({'rain': 'get_model_data_generic(interval=interval, **%s)' % model_dict['rain']})
    hlp.update({'cfc': 'get_model_data_generic(interval=interval, **%s)' % model_dict['cfc']})
    hlp.update({'albedo' : 'get_albedo_data(interval=interval)'})
    hlp.update({'albedo_vis' : 'get_albedo_data_vis(interval=interval)'})
    hlp.update({'albedo_nir' : 'get_albedo_data_nir(interval=interval)'})
    #hlp.update({'sis' : 'get_surface_shortwave_radiation_down(interval=interval)'})
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
    hlp.update({'seaice_concentration': 'get_model_data_generic(interval=interval, **%s)' % model_dict['seaice_concentration']})
    hlp.update({'seaice_extent': 'get_model_data_generic(interval=interval, **%s)' % model_dict['seaice_extent']})
    hlp.update({'gpp': 'get_gpp_data(interval=interval)'})


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




#####################################################################

pl.close('all')

#-----------------------------------------------------------------------


"""
HOWTO

Add a new variable:
1) register variable in get_methods4variables()
2) implement for each data object a routine how to read the data
3) implement an analysis script that performs the actual analysis
4) register this analysis script in get_analysis_scripts()
"""


#######################################################################
# START
# read command line options ...
#######################################################################

if len(sys.argv) > 1:
    if len(sys.argv) == 2:
        file = sys.argv[1] #name of config file
        if not os.path.exists(file):
            raise ValueError, 'Configuration file can not be found: ' + file
    else:
        raise ValueError, 'Currently not more than one command line parameter supported!'

else: #default
    file='pyCMBS.cfg'


######################################################################
# CONFIGURATION and OPTIONS
######################################################################

#/// read configuration file ///
CF = ConfigFile(file)

#/// read plotting options ///
PCFG = PlotOptions(); PCFG.read(CF); plot_options=PCFG



#/// init regions ///
REGIONS = AnalysisRegions()


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


#####################################################################
# TIMES
#####################################################################
s_start_time = CF.start_date #'1983-01-01' #todo where is this used ?
s_stop_time  = CF.stop_date  #'2005-12-31'
start_time   = pl.num2date(pl.datestr2num(s_start_time))
stop_time    = pl.num2date(pl.datestr2num(s_stop_time ))



# get observations? #todo: this needs to be done for each model !!!
model_dict = {'rain':  {'CMIP5':
                            {
                                'variable': 'pr',
                                'unit': 'mm/day',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': 'ensmean',
                                'model_prefix': 'Amon',
                                'file_format' : 'nc',
                                'scale_factor': 86400.,
                                'valid_mask': 'ocean'
                            },


                        'JSBACH_RAW2':
                            {
                                'variable': 'precip_acc',
                                'unit': 'mm/day',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'file_format' : 'nc',
                                'scale_factor': 86400.,
                                'valid_mask': 'global'
                            }
                        },

               'cfc':  {'CMIP5':
                            {
                                'variable': 'clt',
                                'unit': 'per',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': 'ensmean',
                                'model_prefix': 'Amon',
                                'file_format' : 'nc',
                                'scale_factor': 1,
                                'valid_mask': 'global',
                                'custom_path': '/data/workspace/projects/evaclimod/data/CMIP5/clt/merged'
                            },
                        },

              'evap':   {'CMIP5':
                             {
                                 'variable': 'evspsbl',
                                 'unit': 'mm/day',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'model_suffix': 'ensmean',
                                 'file_format' : 'nc',
                                 'model_prefix': 'Amon',
                                 'scale_factor': 86400.,
                                 'valid_mask': 'ocean'
                             }
                        },

              'twpa':   {'CMIP5':
                             {
                                 'variable': 'clwvi',
                                 'unit': 'kg/m^2',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'model_suffix': 'ensmean',
                                 'file_format' : 'nc',
                                 'model_prefix': 'Amon',
                                 'scale_factor': 1.,
                                 'valid_mask': 'ocean'
                             }
              },

             'wind':    {'CMIP5':
                             {
                                 'variable': 'sfcWind',
                                 'unit': 'm/s',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'model_suffix': 'ensmean',
                                 'file_format' : 'nc',
                                 'model_prefix': 'Amon',
                                 'scale_factor': 1.,
                                 'valid_mask': 'ocean'
                             }
                        },

              'wvpa':   {'CMIP5':
                             {
                                 'variable': 'prw',
                                 'unit': 'kg m^2',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'model_suffix': 'ensmean',
                                 'file_format' : 'nc',
                                 'model_prefix': 'Amon',
                                 'scale_factor': 1,
                                 'valid_mask': 'ocean'
                             }
                        },

              'late':   {'CMIP5':
                             {
                                 'variable': 'hfls',
                                 'unit': 'W/m^2',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'model_suffix': 'ensmean',
                                 'file_format' : 'nc',
                                 'model_prefix': 'Amon',
                                 'scale_factor': 1,
                                 'valid_mask': 'ocean'
                             }
                        },

              'hair':   {'CMIP5':
                             {
                                 'variable': 'huss',
                                 'unit': '$kg/kg^2$',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'model_suffix': 'ensmean',
                                 'file_format' : 'nc',
                                 'model_prefix': 'Amon',
                                 'scale_factor': 1,
                                 'valid_mask': 'ocean'
                             }
                        },

              'seaice_concentration':   {'CMIP5':
                               {
                                'variable': 'sic',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': 'ens_mean_185001-200512',
                                'file_format' : 'nc',
                                'model_prefix': 'OImon',
                                'scale_factor': 1,
                                'valid_mask': 'ocean',
                                'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk'
                               },

                           'CMIP3':
                               {
                                'variable': 'SICOMO',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': '1860-2100.ext',
                                'file_format' : 'nc',
                                'model_prefix': '',
                                'scale_factor': 100.,
                                'valid_mask': 'ocean',
                                'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk',
                                'level':0
                               },
                          },

              'seaice_extent':   {'CMIP5':
                               {
                                'variable': 'sic',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': 'ens_mean_185001-200512',
                                'file_format' : 'nc',
                                'model_prefix': 'OImon',
                                'scale_factor': 1,
                                'valid_mask': 'ocean',
                                'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk'
                               },

                           'CMIP3':
                               {
                                'variable': 'SICOMO',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': '1860-2100.ext',
                                'file_format' : 'nc',
                                'model_prefix': '',
                                'scale_factor': 100.,
                                'valid_mask': 'ocean',
                                'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk',
                                'level':0
                               },
                          },

              'budg':   {'CMIP5':
                             {
                                'variable': 'budg',
                                'unit': 'mm/d',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': 'ensmean',
                                'file_format' : 'nc',
                                'model_prefix': 'Amon',
                                'scale_factor': 86400.,
                                'valid_mask': 'ocean',
                                'custom_path' : '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data'
                             }
                        },


              'sis':   {'JSBACH_RAW2':
                             {
                                'variable': 'swdown_acc',
                                'unit': '$W/m^2$',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'file_format' : 'nc',
                                'scale_factor' : 1.,
                                'valid_mask' : 'land'
                             }
                        },

              'surface_upward_flux':   {'JSBACH_RAW2':
                             {
                                'variable': 'swdown_reflect_acc',
                                'unit': '$W/m^2$',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'file_format' : 'nc',
                                'scale_factor' : 1.,
                                'valid_mask' : 'land'
                             }
                        },

              'albedo_vis':   {'JSBACH_RAW2':
                             {
                                'variable': 'var14',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'file_format' : 'nc',
                                'scale_factor' : 1.,
                                'valid_mask' : 'land'
                             }
                        },

              'albedo_nir':   {'JSBACH_RAW2':
                             {
                                'variable': 'var15',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'file_format' : 'nc',
                                'scale_factor' : 1.,
                                'valid_mask' : 'land'
                             }
                        },

              'temperature':   {'JSBACH_RAW2':
                             {
                                'variable': 'temp2',
                                'unit': 'K',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'file_format' : 'nc',
                                'scale_factor' : 1.,
                                'valid_mask' : 'global'
                             }
                        }

              }


########################################################################################################################
# INIT METHODS
########################################################################################################################
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
    elif CF.dtypes[i].upper() == 'JSBACH_RAW2':
        themodel = JSBACH_RAW2(data_dir,varmethods,experiment,intervals=CF.intervals,start_time=start_time,stop_time=stop_time,name=model,shift_lon=shift_lon,model_dict=model_dict)

    elif CF.dtypes[i].upper() == 'CMIP3':
        themodel = CMIP3Data(data_dir,model,experiment,varmethods,intervals=CF.intervals,lat_name='lat',lon_name='lon',label=model,start_time=start_time,stop_time=stop_time,shift_lon=shift_lon)
    else:
        #logger.error('Invalid model type!')
        print CF.dtypes[i]
        raise ValueError, 'Invalid model type!'

    #--- read data for current model ---
    themodel.plot_options = plot_options #options that specify regrid options etc.

    #print themodel.plot_options
    themodel.get_data()


    #--- copy current model to a variable named modelXXXX ---
    cmd = 'model' + str(model_cnt).zfill(4) + ' = ' + 'themodel.copy(); del themodel'
    exec(cmd) #store copy of cmip5 model in separate variable

    #--- append model to list of models ---
    proc_models.append('model' + str(model_cnt).zfill(4))
    model_cnt += 1


########################################################################################################################
# MULTIMODEL MEAN
########################################################################################################################
#--- here we have now all the model and variables read. The list of all models is contained in the variable proc_models.

f_mean_model = False #todo put this as an option!
if f_mean_model:
    #calculate climatological mean values: The models contain already climatological information in the variables[] list. Thus there is not need to take care for the different timesteps here. This should have been handled already in the preprocessing.
    #generate instance of MeanModel to store result
    MEANMODEL = MeanModel(varmethods,intervals=CF.intervals)

    #sum up all models
    for i in range(len(proc_models)):
        exec('actmodel = ' + proc_models[i] + '.copy()')
        MEANMODEL.add_member(actmodel); del actmodel

    #calculate ensemble mean
    MEANMODEL.ensmean()

    #add mean model to general list of models to process in analysis
    proc_models.append('MEANMODEL')


########################################################################################################################
# END MULTIMODEL MEAN
########################################################################################################################




########################################################################################################################
# INIT reporting and plotting and diagnostics
########################################################################################################################
# Gleckler Plot
global_gleckler = GlecklerPlot()

# Report
outdir='.' + os.sep + 'report_' + CF.options['report'] + os.sep
rep = Report(CF.options['report'],'pyCMBS report - ' + CF.options['report'],CF.options['author'],
             outdir=outdir,
             dpi=300,format=CF.options['report_format'])
cmd = 'cp ../logo/Phytonlogo5.pdf ' + rep.outdir
os.system(cmd )


########################################################################################################################
########################################################################################################################
########################################################################################################################
# MAIN ANALYSIS LOOP: perform analysis for each model and variable
########################################################################################################################
########################################################################################################################
########################################################################################################################
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
            cmd = scripts[variable]+'(' + model_list + ',GP=global_gleckler,shift_lon=shift_lon,use_basemap=use_basemap,report=rep,interval=CF.intervals[variable],plot_options=PCFG,regions=REGIONS.regions)'
            eval(cmd) #run analysis


########################################################################################################################
# GLECKLER PLOT finalization ...
########################################################################################################################

#/// generate Gleckler analysis plot for all variables and models analyzed ///
global_gleckler.plot(vmin=-0.1,vmax=0.1,nclasses=25,show_value=True)
oname = outdir + 'gleckler.pkl'
if os.path.exists(oname):
    os.remove(oname)
#pickle.dump(global_gleckler,open(oname,'w')) #store gleckler plot as separate file for further finetuning if necessary

rep.section('Summary error statistics')
rep.figure(global_gleckler.fig,caption='Gleckler et al. (2008) model preformance index')

#/// legend for gleckler plot ///
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



########################################################################################################################
# CLEAN up and finish
########################################################################################################################
rep.close()
pl.close('all')


print '##########################################'
print '# BENCHMARKING FINIHSED!                 #'
print '##########################################'

#~ if __name__ == '__main__': todo
    #~ main()


