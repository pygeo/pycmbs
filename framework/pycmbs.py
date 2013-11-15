# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@mpimet.mpg.de"

'''
# Copyright (C) 2011-2013 Alexander Loew, alexander.loew@mpimet.mpg.de
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




__author__ = "Alexander Loew"
__version__ = "0.1.4"
__date__ = "2013/11/12"

#============ IMPORTS ==================================================

from pyCMBS import *

#--- always use plot backend which is not interactive for benchmarking framework
import matplotlib as mpl
mpl.rcParams['backend'] = 'Agg'
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




#=======================================================================



#####################################################################
#####################################################################
#####################################################################

pl.close('all')

#-----------------------------------------------------------------------




#=======================================================================
#=======================================================================

"""
HOWTO

Add a new variable:
1) register variable in ./configuration/model_data_routines.json
2) implement for each Data object a routine how to read the data
3) implement an analysis script that performs the actual analysis
4) register this analysis script in the file ./configuration/analysis_scripts.json
"""


########################################################################################################################
# START
# read command line options ...
########################################################################################################################

if len(sys.argv) > 1:
    if len(sys.argv) == 2:
        file = sys.argv[1] #name of config file
        if not os.path.exists(file):
            raise ValueError, 'Configuration file can not be found: ' + file
    else:
        raise ValueError, 'Currently not more than one command line parameter supported!'
else: #default
    file='pyCMBS.cfg'




########################################################################################################################
# CONFIGURATION and OPTIONS
########################################################################################################################

#/// read configuration file ///
CF = ConfigFile(file)

#/// read plotting options ///
PCFG = PlotOptions(); PCFG.read(CF); plot_options=PCFG

########################################################################################################################
# REMOVE previous Data warnings
########################################################################################################################
outdir='.' + os.sep + 'report_' + CF.options['report'] + os.sep
os.environ['DATA_WARNING_FILE'] = outdir + 'data_warnings_' + CF.options['report'] + '.log'

if os.path.exists(os.environ['DATA_WARNING_FILE']):
    os.remove(os.environ['DATA_WARNING_FILE'])


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


########################################################################################################################
# TIMES
########################################################################################################################
s_start_time = CF.start_date #'1983-01-01' #todo where is this used ?
s_stop_time  = CF.stop_date  #'2005-12-31'
start_time   = pl.num2date(pl.datestr2num(s_start_time))
stop_time    = pl.num2date(pl.datestr2num(s_stop_time ))




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
scripts = CF.get_analysis_scripts()

#/// get dictionary with methods how to read data for variables to be analyzed ///
variables  = CF.variables
varmethods = CF.get_methods4variables(variables, model_dict)







#/// READ DATA ///
"""
create a Model instance for each model specified
in the configuration file

read the data for all variables and return a list
of Data objects for further processing
"""

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

f_mean_model = True #todo put this as an option!
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
pickle.dump(global_gleckler.models,open(outdir + 'gleckler_models.pkl','w'))
pickle.dump(global_gleckler.variables,open(outdir + 'gleckler_variables.pkl','w'))
pickle.dump(global_gleckler.data,open(outdir + 'gleckler_data.pkl','w'))
pickle.dump(global_gleckler._raw_data,open(outdir + 'gleckler_rawdata.pkl','w'))

rep.section('Summary error statistics')
rep.subsection('Gleckler metric')
rep.figure(global_gleckler.fig,caption='Gleckler et al. (2008) model performance index',width='10cm')


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
            if varoptions[k]['add_to_report']: #only add observation to legend, if option in INI file is set
                thelabels.update({int(varoptions[k]['gleckler_position']) : k}) #generate dictionary for GlecklerPLot legend
    fl = global_gleckler._draw_legend(thelabels,title=variable.upper())
    rep.figure(fl,width='8cm',bbox_inches=None)
    del fl, thelabels


#/// plot model ranking between different observational datasets ///
rep.subsection('Model ranking consistency')
for v in global_gleckler.variables:
    rep.subsubsection(v.upper())
    tmpfig = global_gleckler.plot_model_ranking(v,show_text=True)
    rep.figure(tmpfig,width='8cm',bbox_inches=None,caption='Model RANKING for different observational datasets: ' + v.upper())
    del tmpfig

    #/// plot absolut model error

    tmpfig = global_gleckler.plot_model_error(v)
    rep.figure(tmpfig,width='8cm',bbox_inches=None,caption='Model ERROR for different observational datasets: ' + v.upper())
    del tmpfig




#---


########################################################################################################################
# CLEAN up and finish
########################################################################################################################
pl.close('all')
rep.close()




print '##########################################'
print '# BENCHMARKING FINIHSED!                 #'
print '##########################################'

#~ if __name__ == '__main__': todo
    #~ main()


