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

from utils import *
from external_analysis import *
from datetime import *
from dateutil.rrule import *
from cdo import *



#///////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////


def preprocess_seasonal_data(raw_file,interval=None,themask = None,force=False,obs_var=None,label='',shift_lon=None):
    """
    This subroutine performs pre-processing of some raw data file. It does

    1) generate monthly timeseries
    2) generate seasonal timeseries
    3) remaps data to T63 grid
    4) calculates standard deviation and number of valid data
    5) returns seasonal/monthly climatology as well as monthly mean raw data

    ALL generated data will be stored in temporary data directory

    @param raw_file: name of original data file to be processed
    @type raw_file: str

    @param interval: [season,monthly]
    @type interval: str

    @param force: force caluclations
    @type force: bool

    @param themask: mask to be applied to the data (default=None)
    @type themask: numpy bool array or Data

    @param obs_var: variable to analyze
    @type obs_var: str
    """

    sys.stdout.write(' *** Preprocessing ' + raw_file)

    if obs_var == None:
        raise ValueError, 'Name of variable to be processed needs to be specified!'
    if shift_lon == None:
        raise ValueError, 'Lon shift parameter needs to be specified!'

    #--- PREPROCESSING of observational data  ---
    cdo = Cdo()

    #1) generate monthly mean file projected to T63
    obs_mon_file     = get_temporary_directory() + os.path.basename(raw_file)
    obs_mon_file = obs_mon_file[:-3] + '_monmean.nc'
    cdo.monmean(options='-f nc',output=obs_mon_file,input='-remapcon,t63grid ' + raw_file,force=force)

    #2) generate monthly mean or seasonal mean climatology as well as standard deviation
    if interval == 'monthly':
        obs_ymonmean_file     = obs_mon_file[:-3] + '_ymonmean.nc'
        obs_ymonstd_file = obs_mon_file[:-3] + '_ymonstd.nc'
        obs_ymonsum_file = obs_mon_file[:-3] + '_ymonsum.nc'
        obs_ymonN_file   = obs_mon_file[:-3] + '_ymonN.nc'
        cdo.ymonmean(options='-f nc -b 32',output=obs_ymonmean_file,     input=obs_mon_file,force=force)
        cdo.ymonsum (options='-f nc -b 32',output=obs_ymonsum_file, input=obs_mon_file,force=force)
        cdo.ymonstd (options='-f nc -b 32',output=obs_ymonstd_file, input=obs_mon_file,force=force)
        cdo.div(options='-f nc -b 32',output = obs_ymonN_file,input=obs_ymonsum_file + ' ' + obs_ymonmean_file, force=force) #number of samples

        time_cycle = 12
    elif interval == 'season':
        obs_ymonmean_file     = obs_mon_file[:-3] + '_yseasmean.nc'
        obs_ymonstd_file = obs_mon_file[:-3] + '_yseasstd.nc'
        obs_ymonsum_file = obs_mon_file[:-3] + '_yseassum.nc'
        obs_ymonN_file   = obs_mon_file[:-3] + '_yseasN.nc'
        cdo.yseasmean(options='-f nc -b 32',output=obs_ymonmean_file,     input=obs_mon_file,force=force)
        cdo.yseassum (options='-f nc -b 32',output=obs_ymonsum_file, input=obs_mon_file,force=force)
        cdo.yseasstd (options='-f nc -b 32',output=obs_ymonstd_file, input=obs_mon_file,force=force)
        cdo.div(options='-f nc -b 32',output = obs_ymonN_file,input=obs_ymonsum_file + ' ' + obs_ymonmean_file, force=force) #number of samples

        time_cycle = 4

    else:
        print interval
        raise ValueError, 'Unknown temporal interval. Can not perform preprocessing! '

    #--- READ DATA ---
    obs     = Data(obs_ymonmean_file,obs_var,read=True,label=label,unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,time_cycle=time_cycle)
    obs_std = Data(obs_ymonstd_file,obs_var,read=True,label=label + ' std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,time_cycle=time_cycle) #,mask=ls_mask.data.data)
    obs.std = obs_std.data.copy(); del obs_std
    obs_N = Data(obs_ymonN_file,obs_var,read=True,label=label + ' N',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,time_cycle=time_cycle) #,mask=ls_mask.data.data)
    obs.n = obs_N.data.copy(); del obs_N

    #sort climatology to be sure that it starts always with January
    obs.adjust_time(year=1700,day=15) #set arbitrary time for climatology
    obs.timsort()

    #read monthly data (needed for global means and hovmoeller plots)
    obs_monthly = Data(obs_mon_file,obs_var,read=True,label=label,unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,time_cycle=12) #,mask=ls_mask.data.data)

    #/// center dates of months
    obs_monthly.adjust_time(day=15)
    obs.adjust_time(day=15)

    if themask != None:
        obs._apply_mask(themask)
        obs_monthly._apply_mask(themask)

    return obs,obs_monthly

    #--- preprocessing END ---


#///////////////////////////////////////////////////////////////////////
# ANALYSIS SECTION
#///////////////////////////////////////////////////////////////////////

#-----------------------------------------------------------------------------------------------------------------------

def evaporation_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'

    report.section('Evaporation')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'evap', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------

def rainfall_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):

    if report == None:
        raise ValueError, 'You need to specify report option!'

    report.section('Precipitation')

#    report.subsection('HOAPS')
#    generic_analysis(obs_dict, model_list, 'rain', 'HOAPS', GP = GP, report = report, use_basemap = use_basemap, shift_lon = shift_lon)

    report.subsection('GPCP')
    generic_analysis(obs_dict, model_list, 'rain', 'GPCP', GP = GP, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    report.subsection('CRU')
    generic_analysis(obs_dict, model_list, 'rain', 'CRU', GP = GP, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    #report.subsection('GPCC') #todo activate
    #generic_analysis(obs_dict, model_list, 'rain', 'GPCC', GP = GP, report = report, use_basemap = use_basemap, shift_lon = shift_lon)

    #todo add TMPA data

#-----------------------------------------------------------------------------------------------------------------------

# WIND ANALYSIS
def wind_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'

    report.section('10 meters surface winds')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'wind', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------

#TWPA analysis
def twpa_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'
    report.section('Total column water content')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'twpa', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------

#WVPA analysis
def wvpa_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'
    report.section('Water Vapor Path')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'wvpa', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------

#HAIR analysis
def hair_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'
    report.section('Surface specific humidity')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'hair', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------

#LATE analysis
def late_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'
    report.section('Upward latent heat flux')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'late', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------

#BUDG analysis
def budg_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
    if report == None:
        raise ValueError, 'You need to specify report option!'
    report.section('Upward freshwater flux')
    report.subsection('HOAPS')
    generic_analysis(obs_dict, model_list, 'budg', 'HOAPS', GP=GP,report=report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

#-----------------------------------------------------------------------------------------------------------------------













#-------------------------------------------------------------------------------------------------------------



obs_dict = {
    'rain':
        {
            'CRU':
             {'obs_file':get_data_pool_directory() + 'variables/land/precipitation/CRU/cru_ts_3_00.1901.2006.pre_miss.nc',
              'start_time': '1979-01-01',
              'stop_time': '2006-12-31',
              'obs_var': 'pre',
              'scale_data': 12./365.,
              'glecker_position': 4,
              'map_difference': True,
              'map_seasons': True,
              'reichler_plot': False,
              'glecker_plot': True,
              'units': '$mm/day$',
              'mask_area': 'ocean',
              'add_to_report': True,
              #'interval': 'monthly',
              'label': 'CRU',
              'preprocess': True,
              'parameter': 'rain',
              'report': True
             }, #end CRU


            'GPCC':
                 {'obs_file':get_data_pool_directory() + 'variables/land/precipitation/GPCC/gpcc_full_vs4_1951-2007.nc',
                  'start_time': '1979-01-01',
                  'stop_time': '2007-12-31',
                  'obs_var': 'var260',
                  'scale_data': 12./365.,
                  'glecker_position': 3,
                  'map_difference': True,
                  'map_seasons': True,
                  'reichler_plot': False,
                  'glecker_plot': True,
                  'units': '$mm/day$',
                  'mask_area': 'ocean',
                  'add_to_report': True,
                  #'interval': 'monthly',
                  'label': 'GPCC',
                  'preprocess': True,
                  'parameter': 'rain',
                  'report': True
                 }, #end GPCC


            'HOAPS':
                  {'obs_file':'/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/rain/hoaps-g.t63.m01.rain.1987-2008.nc',
                   'start_time': '1989-01-01',
                   'stop_time': '2007-12-31',
                   'obs_var': 'rain',
                   'scale_data': 1.,
                   'glecker_position': 1,
                   'map_difference': True,
                   'map_seasons': True,
                   'reichler_plot': False,
                   'glecker_plot': True,
                   'units': '$mm/day$',
                   'mask_area': 'ocean',
                   'add_to_report': True,
                   #'interval': 'monthly',
                   'label': 'HOAPS',
                   'preprocess': True,
                   'parameter': 'rain',
                   'report': True
                  }, #end HOAPS

            'GPCP':
                   {'obs_file':get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasmean.nc',
                    #'obs_file_std':get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasstd.nc',
                    'start_time': '1989-01-01',  #todo clarify time
                    'stop_time': '2007-12-31',
                    'obs_var': 'precip',
                    'scale_data': 1.,
                    'glecker_position': 2,
                    'map_difference': True,
                    'map_seasons': True,
                    'reichler_plot': False,
                    'glecker_plot': True,
                    'units': '$mm/day$',
                    'mask_area': 'none',
                    'add_to_report': True,
                    #'interval': 'monthly',
                    'label': 'GPCP',
                    'preprocess': True,
                    'parameter': 'rain',
                    'report': True
                   } #end GPCP

            }, #end rain



    'evap':
        {'HOAPS':
             {'obs_file':'/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/evap/hoaps-g.t63.m01.evap.1987-2008.nc',
              'start_time': '1989-01-01',
              'stop_time': '2007-12-31',
              'obs_var': 'evap',
              'scale_data': 1.,
              'glecker_position': 1,
              'map_difference': True,
              'map_seasons': True,
              'reichler_plot': True,
              'glecker_plot': True,
              'units': '$mm/day$',
              'mask_area': 'ocean',
              'add_to_report': True,
              #'interval': 'monthly',
              'label': 'Daily evaporation',
              'preprocess': True,
              'parameter': 'evap',
              'report': True
             }
        }, #end EVAP


    'wind':
        {'HOAPS':
             {'obs_file': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/wind/hoaps-g.t63.m01.wind.1987-2008.nc',
              'start_time': '1989-01-01',
              'stop_time': '2007-12-31',
              'obs_var': 'wind',
              'scale_data': 1.,
              'glecker_position': 1,
              'map_difference': True,
              'map_seasons': True,
              'reichler_plot': True,
              'glecker_plot': True,
              'units': 'm/s',
              'mask_area': 'ocean',
              'add_to_report': True,
              #'interval': 'monthly',
              'label': '10m surface winds',
              'preprocess': True,
              'parameter': 'wind',
              'report': True
             }
        },

    'twpa':
        {'HOAPS':
             { 'obs_file': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/twpa/hoaps-g.t63.m01.twpa.1987-2008.nc',
               'start_time': '1989-01-01',
               'stop_time': '2007-12-31',
               'obs_var': 'twpa',
               'scale_data': 1.,
               'glecker_position': 1,
               'map_difference': True,
               'map_seasons': True,
               'reichler_plot': True,
               'glecker_plot': True,
               'units': '$kg m^2$',
               'mask_area': 'ocean',
               'add_to_report': True,
               #'interval': 'monthly',
               'label': 'Total water path',
               'preprocess': True,
               'parameter': 'twpa',
               'report': True
             }
        },

    'wvpa': {'HOAPS':
                 {'obs_file': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/wvpa/hoaps-g.t63.m01.wvpa.1987-2008.nc',
                  'start_time': '1989-01-01',
                  'stop_time': '2007-12-31',
                  'obs_var': 'wvpa',
                  'scale_data': 1.,
                  'glecker_position': 1,
                  'map_difference': True,
                  'map_seasons': True,
                  'reichler_plot': True,
                  'glecker_plot': True,
                  'units': 'kg/m^2',
                  'mask_area': 'ocean',
                  'add_to_report': True,
                  #'interval': 'monthly',
                  'label': 'Water Vapor Path',
                  'preprocess': True,
                  'parameter': 'wvpa',
                  'report': True
                 }
    },


    'hair': {'HOAPS':
                 {'obs_file': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/hair/hoaps-g.t63.m01.hair.1987-2008.nc',
                  'start_time': '1989-01-01',
                  'stop_time': '2007-12-31',
                  'obs_var': 'hair',
                  'scale_data': 1.,
                  'glecker_position': 1,
                  'map_difference': True,
                  'map_seasons': True,
                  'reichler_plot': True,
                  'glecker_plot': True,
                  'units': 'g/kg',
                  'mask_area': 'ocean',
                  'add_to_report': True,
                  #'interval': 'monthly',
                  'label': 'Surface specific humidity',
                  'preprocess': True,
                  'parameter': 'hair',
                  'report': True
                 }
    },

    'late': {'HOAPS':
                 {'obs_file': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/late/hoaps-g.t63.m01.late.1987-2008.nc',
                  'start_time': '1989-01-01',
                  'stop_time': '2007-12-31',
                  'obs_var': 'late',
                  'scale_data': 1.,
                  'glecker_position': 1,
                  'map_difference': True,
                  'map_seasons': True,
                  'reichler_plot': True,
                  'glecker_plot': True,
                  'units': 'W/m^2',
                  'mask_area': 'ocean',
                  'add_to_report': True,
                  #'interval': 'monthly',
                  'label': 'Upward latent heat flux',
                  'preprocess': True,
                  'parameter': 'late',
                  'report': True
                 }
             },



     'budg': {'HOAPS':
                  {'obs_file': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data/budg/hoaps-g.t63.m01.budg.1987-2008.nc',
                   'start_time': '1989-01-01',
                   'stop_time': '2007-12-31',
                   'obs_var': 'budg',
                   'scale_data': 1.,
                   'glecker_position': 1,
                   'map_difference': True,
                   'map_seasons': True,
                   'reichler_plot': True,
                   'glecker_plot': True,
                   'units': 'mm/d',
                   'mask_area': 'ocean',
                   'add_to_report': True,
                   #'interval': 'monthly',
                   'label': 'Upward freshwater flux',
                   'preprocess': True,
                   'parameter': 'budg',
                   'report': True}
     },


     'sis':{'ISCCP': {'obs_file': get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/isccp/T63_ISCCP__versD1__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc',
            'start_time': '1984-01-01',
            'stop_time': '2005-12-31',
            'obs_var': 'BfISC84',
            'scale_data': 1.,
            'glecker_position': 3,
            'map_difference': True,
            'map_seasons': True,
            'reichler_plot': True,
            'glecker_plot': True,
            'units': '$W/m^2$',
            'mask_area': 'none',
            'add_to_report': True,
            #'interval': 'monthly',
            'label': 'Shortwave downward radiation flux in air',
            'preprocess': True,
            'parameter': 'sis',
            'report': True},

            'SRB': {'obs_file': get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/srb/T63_SRB__vers28__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc',
                      'start_time': '1984-01-01',
                      'stop_time': '2005-12-31',
                      'obs_var': 'BfSRB84',
                      'scale_data': 1.,
                      'glecker_position': 4,
                      'map_difference': True,
                      'map_seasons': True,
                      'reichler_plot': True,
                      'glecker_plot': True,
                      'units': '$W/m^2$',
                      'mask_area': 'none',
                      'add_to_report': True,
                      #'interval': 'monthly',
                      'label': 'Shortwave downward radiation flux in air',
                      'preprocess': True,
                      'parameter': 'sis',
                      'report': True},

            'CMSAF': {'obs_file': get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/cmsaf_sis/SISmm_all_t63.nc',
                    'start_time': '1984-01-01',
                    'stop_time': '2005-12-31',
                    'obs_var': 'SIS',
                    'scale_data': 1.,
                    'glecker_position': 1,
                    'map_difference': True,
                    'map_seasons': True,
                    'reichler_plot': True,
                    'glecker_plot': True,
                    'units': '$W/m^2$',
                    'mask_area': 'none',
                    'add_to_report': True,
                    #'interval': 'monthly',
                    'label': 'Shortwave downward radiation flux in air',
                    'preprocess': True,
                    'parameter': 'sis',
                    'report': True},

            'CERES': {'obs_file': get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres_ebaf2.6/CERES_EBAF-Surface__Ed2.6r__sfc_sw_down_all_mon__1x1__200003-201002.nc',
                      'start_time': 'xxxx-01-01',
                      'stop_time': 'xxxx-12-31',
                      'obs_var': 'sfc_sw_down_all_mon',
                      'scale_data': 1.,
                      'glecker_position': 2,
                      'map_difference': True,
                      'map_seasons': True,
                      'reichler_plot': True,
                      'glecker_plot': True,
                      'units': '$W/m^2$',
                      'mask_area': 'none',
                      'add_to_report': True,
                      #'interval': 'monthly',
                      'label': 'Shortwave downward radiation flux in air',
                      'preprocess': True,
                      'parameter': 'sis',
                      'report': True}
         },

    'albedo' : {'MODIS':{
                'obs_file': get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged.nc',
                'stop_time': 'xxxx-12-31',
                'obs_var': 'surface_albedo_WSA',
                'scale_data': 1.,
                'glecker_position': 1,
                'map_difference': True,
                'map_seasons': True,
                'reichler_plot': True,
                'glecker_plot': True,
                'units': '-',
                'mask_area': 'none',
                'add_to_report': True,
                #'interval': 'monthly',
                'label': 'Surface albedo',
                'preprocess': True,
                'parameter': 'albedo',
                'report': True,
                'vmin':0.,
                'vmax':0.6,
                'dmin':-0.09,
                'dmax':0.09}

                }


    } #end of dict


global_settings_dict = {'landsea_mask':
                            {'filename': ''}}

#=======================================================================
# GENERIC - start
#=======================================================================

def generic_analysis(obs_dict, model_list, obs_type, obs_name, GP=None, GM = None, shift_lon=False, use_basemap=False, report=None,interval=None):
    """
    function for performing common analysis actions
    it is not a parameter specific function
    use it as a template for specific analysis

    @param      obs_dict dictionary that contains settings
                and parameters from the config file
    """

    if interval not in ['monthly','season']:
        raise ValueError, 'invalid interval in generic_analysis()'


    local_obs_dict = obs_dict[obs_type][obs_name]

    # get observations information from the global dictionary
    #s_start_time        = local_obs_dict['start_time'] #todo appropriate usage of start/stop times! (from config file ???)
    #s_stop_time         = local_obs_dict['stop_time']
    obs_raw             = local_obs_dict['obs_file']
    obs_var             = local_obs_dict['obs_var']
    #obs_type = local_obs_dict['obs_provider']
    scale_data          = local_obs_dict['scale_data'] #todo: how to best use ???
    #interval            = local_obs_dict['interval'] #does not make sense to take this from dictionary, as information needs to come from config file!!!
    #mask_area           = local_obs_dict['mask_area']
    glecker_pos         = local_obs_dict['glecker_position']
    param_name          = local_obs_dict['parameter'] # model variable name ##todo is this really needed ????; obs_name should do ????

    if 'vmin' in local_obs_dict.keys():
        vmin = local_obs_dict['vmin']
    else:
        vmin = None
    if 'vmax' in local_obs_dict.keys():
        vmax = local_obs_dict['vmax']
    else:
        vmax = None

    if 'dmin' in local_obs_dict.keys():
        dmin = local_obs_dict['dmin']
    else:
        dmin = None
    if 'dmax' in local_obs_dict.keys():
        dmax = local_obs_dict['dmax']
    else:
        dmax = None




    #model = model_list[0]
    #m_data = param_name
    m_data_org = param_name + '_org'

    if report == None:
        raise ReportError("Report option was not enabled")


    ls_mask = get_T63_landseamask(shift_lon, area = local_obs_dict['mask_area'])

    # preprocessing
    if local_obs_dict['preprocess'] == True:
        print interval
        obs_orig, obs_monthly = preprocess_seasonal_data(obs_raw, interval = interval,  themask = ls_mask, force = False, obs_var = obs_var, label = local_obs_dict['label'], shift_lon = shift_lon)


    #--- initialize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    if GM == None:
        fG = plt.figure(); axg = fG.add_subplot(211); axg1 = fG.add_subplot(212)
        GM = GlobalMeanPlot(ax=axg,ax1=axg1,climatology=False) #global mean plot

    if GM != None:
        GM.plot(obs_monthly, linestyle = '--')

    if local_obs_dict['map_seasons'] == True: #seasonal mean plot
        f_season = map_season(obs_orig,use_basemap=use_basemap,cmap_data='jet',show_zonal=True,zonal_timmean=True,nclasses=6,vmin=vmin,vmax=vmax)
        report.figure(f_season,caption='Seasonal mean ' + obs_name)


    for model in model_list:

        sys.stdout.write('\n *** %s analysis of model: ' % (param_name) + model.name + "\n")
        model_data = model.variables[param_name]
        GP.add_model(model.name) #register model name in GlecklerPlot

        if local_obs_dict['report'] == True:
            #/// report results
            sys.stdout.write('\n *** Making report figures. \n')
            report.subsubsection(model.name)

        if GM != None:
            if m_data_org in model.variables.keys():
                GM.plot(model.variables[m_data_org][2],label=model.name) #(time,meandata) replace rain_org with data_org

        if model_data == None:
            print 'Data not existing for model ', model.name; continue

        if model_data.data.shape != obs_orig.data.shape:
            print 'Inconsistent geometries' # add here parameter name
            print 'Model: ', model_data.data.shape
            print 'Observation: ', obs_orig.data.shape
            raise ValueError, "Invalid geometries"

        if local_obs_dict['map_difference'] == True:
            sys.stdout.write('\n *** Map difference plotting. \n')
            #--- generate difference map
            #pdb.set_trace()
            f_dif  = map_difference(model_data, obs_orig, nclasses=6,use_basemap=use_basemap,show_zonal=True,zonal_timmean=False,dmin=dmin,dmax=dmax,vmin=vmin,vmax=vmax)
            report.figure(f_dif,caption='Mean and relative differences')

        if local_obs_dict['map_seasons'] == True:
            sys.stdout.write('\n *** Seasonal maps plotting\n')

            #seasonal map
            f_season = map_season(model_data,use_basemap=use_basemap,cmap_data='jet',show_zonal=True,zonal_timmean=True,nclasses=6,vmin=vmin,vmax=vmax)
            report.figure(f_season,caption='Seasonal means model')

        if local_obs_dict['reichler_plot'] == True:
            #/// Reichler statistics ///
            sys.stdout.write('\n *** Computing diagnostics (Reichler index). \n')
            Diag = Diagnostic(obs_orig, model_data)
            e2   = Diag.calc_reichler_index()
            #print e2
            Rplot.add(e2,model_data.label,color='red')

        if local_obs_dict['glecker_plot'] == True:
            #/// Gleckler plot ///
            sys.stdout.write('\n *** Glecker plot. \n')
            e2a = GP.calc_index(obs_orig,model_data,model,param_name)
            #e2a = 0
            GP.add_data(param_name,model.name,e2a,pos=glecker_pos)

    del obs_monthly
    sys.stdout.write('\n *** Reichler plot.\n')
    f_reich = Rplot.bar(title='relative model error: %s' % param_name.upper())
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()

    sys.stdout.write('\n *** Processing finished. \n')

#=======================================================================
# GENERIC - end
#=======================================================================
























#=======================================================================
# VEGETATION COVER FRACTION -- begin
#=======================================================================

def phenology_faPAR_analysis(model_list,GP=None,shift_lon=None,use_basemap=False,report=None):
    """
    Analysis of vegetation phenology

    @param model_list: list of model C{Data} objects
    @param GP: Gleckler plot instance
    @param shift_lon: shift data in longitude direction
    @param use_basemap: use basemap for the analysis
    @param report: instance of Report

    Reference:
    Dahlke, C., Loew, A. & Reick, C., 2012. Robust identification of global greening phase patterns from remote sensing vegetation products. Journal of Climate, p.120726140719003. Available at: http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-11-00319.1 [Accessed December 6, 2012].

    """

    #/// OPTIONS
    f_execute = True #execute scripts

    #script names
    template1 = './external/phenology_benchmarking/Greening_Phase_Analysis_I.m'
    template2 = './external/phenology_benchmarking/Greening_Phase_Analysis_II.m'
    template3 = './external/phenology_benchmarking/Greening_Phase_Analysis_III.m'
    template4 = './external/phenology_benchmarking/Greening_Phase_Analysis_IV.m'

    if GP == None:
        GP = GlecklerPlot()

    for model in model_list:

        print '*** PHENOLOGY analysis for model: ' + model.name

        fapar_data = model.variables['phenology_faPAR'] #data object for phenology variable
        snow_data  = model.variables['snow']            #data object for snow variable

        #/// output directory for analysis results of current model
        outdir=os.getcwd() + '/output/GPA-01-' + model.name.replace(' ','')

        #/// directory where model simulations are performed and templates are copied to
        rundir = './output/tmp_' + model.name.replace(' ','') + '/'

        #/// Gleckler plot
        GP.add_model(model.name)

        #/// file and variable names ///
        data_file = fapar_data.filename
        varname = fapar_data.varname

        snowfilename = snow_data.filename
        snowvar = snow_data.varname


        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP1
        #//////////////////////////////////////
        tags = []
        tags.append({'tag':'<STARTYEAR>','value':'1995'})
        tags.append({'tag':'<STOPYEAR>','value':'2005'})
        tags.append({'tag':'<OUTDIR>','value':outdir })
        tags.append({'tag':'<FFTFIGNAME>','value':'FFT-Mask-'+model.name.replace(' ','')})
        tags.append({'tag':'<INPUTDATAFILE>','value':data_file})
        tags.append({'tag':'<DATAVARNAME>','value':varname})
        tags.append({'tag':'<FFTMASKFILE>','value':outdir+'/fft_mask.mat'})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>',template1,tags,options=',quit',output_directory=rundir ) #temporary scripts are stored in directories that have same name as model
        E.run(execute=f_execute,remove_extension=True)

        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP2
        #//////////////////////////////////////
        tags.append({'tag':'<INPUTSNOWFILE>','value': snowfilename })
        tags.append({'tag':'<SNOWVARNAME>','value':snowvar})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>',template2,tags,options=',quit',output_directory=rundir )
        E.run(execute=f_execute,remove_extension=True)

        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP3
        #//////////////////////////////////////
        tags.append({'tag':'<INPUTDIRECTORY>','value':outdir + '/'})
        tags.append({'tag':'<SENSORMASKFILE>','value':os.getcwd() + '/external/phenology_benchmarking/sensor_mask.mat'})
        tags.append({'tag':'<RESULTDIR_AVHRR>'  ,'value':'/net/nas2/export/eo/workspace/m300028/GPA/sensors/AVH_T63'}) #@todo this needs to compe from POOL SEP !!!!
        tags.append({'tag':'<RESULTDIR_SEAWIFS>','value':'/net/nas2/export/eo/workspace/m300028/GPA/sensors/SEA_T63'})
        tags.append({'tag':'<RESULTDIR_CYCLOPES>','value':'/net/nas2/export/eo/workspace/m300028/GPA/sensors/CYC_T63'})
        tags.append({'tag':'<RESULTDIR_MODIS>','value':'/net/nas2/export/eo/workspace/m300028/GPA/sensors/MCD_T63'})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>',template3,tags,options=',quit',output_directory=rundir )
        E.run(execute=f_execute,remove_extension=True)

        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP4
        #//////////////////////////////////////
        tags.append({'tag':'<FIG1CAPTION>','value':'Fig1_FFT_Mask_Evergreen_Overview'})
        tags.append({'tag':'<FIG2CAPTION>','value':'Fig2_GPAII_Results'})
        tags.append({'tag':'<FIG3CAPTION>','value':'Fig3_GPAIII_Shifts_between_Model_Sensors'})
        tags.append({'tag': '<FIG4CAPTION>', 'value': 'Fig4_Time_series_from_various_biomes'})
        tags.append({'tag':'<RUNDIR>','value':rundir})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>',template4,tags,options=',quit',output_directory=rundir )
        E.run(execute=f_execute,remove_extension=True)

        #/// Gleckler plot ///
        #~ e2a = GP.calc_index(gpcp,model_data,model,'pheno_faPAR')
        e2a = 0. #@todo
        GP.add_data('pheno_faPAR',model.name,e2a,pos=1)









def grass_fraction_analysis(model_list):
    #use same analysis script as for trees, but with different argument
    tree_fraction_analysis(model_list,pft='grass')



def tree_fraction_analysis(model_list,pft='tree'):

    def fraction2netcdf(pft):
        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/' + pft + '_05.dat'
        #save 0.5 degree data to netCDF file
        t=pl.loadtxt(filename)
        outname = filename[:-4] + '.nc'
        if os.path.exists(outname):
            os.remove(outname) #todo: really desired?
        F = Nio.open_file(outname,'w')
        F.create_dimension('lat',t.shape[0])
        F.create_dimension('lon',t.shape[1])
        var = F.create_variable(pft + '_fraction','d',('lat','lon'))
        lat = F.create_variable('lat','d',('lat',))
        lon = F.create_variable('lon','d',('lon',))
        lon.standard_name = "grid_longitude"
        lon.units = "degrees"
        lon.axis = "X"
        lat.standard_name = "grid_latitude"
        lat.units = "degrees"
        lat.axis = "Y"
        var._FillValue = -99.
        var.assign_value(t)
        lat.assign_value(np.linspace(90.-0.25,-90.+0.25,t.shape[0])) #todo: center coordinates or corner coordinates?
        lon.assign_value(np.linspace(-180.+0.25,180.-0.25,t.shape[1]))

        #todo: lon 0 ... 360 ****

        F.close()

        return outname


    #//// load tree fraction observations ////
    #1) convert to netCDF
    fraction_file = fraction2netcdf(pft)
    #2) remap to T63
    y1 = '2001-01-01'; y2='2001-12-31'
    cdo = pyCDO(fraction_file,y1,y2)
    t63_file = cdo.remap(method='remapcon')

    #3) load data
    ls_mask = get_T63_landseamask()
    hansen  = Data(t63_file,pft + '_fraction',read=True,label='MODIS VCF ' + pft + ' fraction',unit='-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)

    #~ todo now remap to T63 grid; which remapping is the most useful ???


    #//// PERFORM ANALYSIS ///
    vmin = 0.; vmax = 1.
    for model in model_list:

        model_data = model.variables[pft]

        if model_data.data.shape != hansen.data.shape:
            print 'WARNING Inconsistent geometries for GPCP'
            print model_data.data.shape; print hansen.data.shape

        dmin=-1; dmax=1.
        dif = map_difference(model_data,hansen,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0.,0.25,0.5,0.75,1.])

        #/// ZONAL STATISTICS
        #--- append zonal plot to difference map
        ax1 = dif.get_axes()[0]; ax2 = dif.get_axes()[1]; ax3 = dif.get_axes()[2]

        #create net axis for zonal plot
        #http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
        divider1 = make_axes_locatable(ax1); zax1 = divider1.new_horizontal("50%", pad=0.05, axes_class=maxes.Axes,pack_start=True)
        dif.add_axes(zax1)

        divider2 = make_axes_locatable(ax2); zax2 = divider2.new_horizontal("50%", pad=0.05, axes_class=maxes.Axes,pack_start=True)
        dif.add_axes(zax2)

        divider3 = make_axes_locatable(ax3); zax3 = divider3.new_horizontal("50%", pad=0.05, axes_class=maxes.Axes,pack_start=True)
        dif.add_axes(zax3)

        #--- calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,xlim=[vmin,vmax])
        zon2 = ZonalPlot(ax=zax2); zon2.plot(hansen,xlim=[vmin,vmax])
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(hansen))
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #--- calculate RMSE statistics etc.
        ES = CorrelationAnalysis(model_data,hansen)
        ES.do_analysis()


#=======================================================================
# VEGETATION COVER FRACTION -- end
#=======================================================================





#=======================================================================
# ALBEDO -- begin
#=======================================================================

def surface_upward_flux_analysis(model_list,GP=None,shift_lon=None,use_basemap=False,report=None,interval='season'):

    if shift_lon == None:
        raise ValueError, 'You need to specify shift_lon option!'
    if use_basemap == None:
        raise ValueError, 'You need to specify use_basemap option!'
    if report == None:
        raise ValueError, 'You need to specify report option!'

    print
    print '************************************************************'
    print '* BEGIN SURFACE UPWARD FLUX analysis ...'
    print '************************************************************'

    report.section('Surface upward flux')

    fG = plt.figure(); axg = fG.add_subplot(211); axg1 = fG.add_subplot(212)
    GM = GlobalMeanPlot(ax=axg,ax1=axg1) #global mean plot

    #- MODIS white sky albedo
    report.subsection('xxxxxxxxxxxxxx')
    report.write('The upward flux observations are calculated using the downward shortwave radiation flux from CERES and the surface albedo from MODIS white sky albedo (WSA).')
    surface_upward_flux_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,interval=interval,obs_type='MODIS_ceres',GM=GM)


    #DIRECT COMPARISON WITH CERES UPWARD FLUX
    report.subsection('CERES upward flux')
    surface_upward_flux_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,interval=interval,obs_type='CERES',GM=GM)



    #CERES ALBEDO scaled with MODEL downward


    #- CERES surface albedo from all sky fluxes
    #report.subsection('CERES albedo')
    #report.write('The CERES surface albedo is calculated as the ratio of the upward and downward surface all sky shortwave radiation fluxes based on CERES EBAF v2.6.' )
    #albedo_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,interval=interval,obs_type='CERES',GM=GM)

    report.figure(fG,caption='Global means for land surface upward flux')

    print '************************************************************'
    print '* END SURFACE UPWARD FLUX analysis ...'
    print '************************************************************'
    print



def surface_upward_flux_analysis_plots(model_list,GP=None,shift_lon=None,use_basemap=False,report=None,interval=None,obs_type=None,GM=None):
    """
    model_list = list which contains objects of data type MODEL
    """

    vmin = 0.; vmax = 150.

    if interval == None:
        raise ValueError, 'Interval period in surface_upward_flux_analyis not specified!'
    if obs_type == None:
        raise ValueError, 'Observation type for surface_upward_flux_analyis was not specified!'

    print 'Doing surface_upward_flux_analyis analysis ...'

    #--- GlecklerPlot
    if GP == None:
        GP = GlecklerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- loading observation data
    if obs_type == 'CERES':
        #CERES EBAF ...
        up_file   = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres_ebaf2.6/CERES_EBAF-Surface__Ed2.6r__sfc_sw_up_all_mon__1x1__200003-201002.nc'
        obs_raw_file = up_file
        obs_var = 'sfc_sw_up_all_mon'
        gleckler_pos = 1

    elif obs_type == 'MODIS_ceres':
        cdo = Cdo()

        down_file = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres_ebaf2.6/CERES_EBAF-Surface__Ed2.6r__sfc_sw_down_all_mon__1x1__200003-201002.nc'

        #remap CERES data first
        obs_down, obs_down_monthly = preprocess_seasonal_data(down_file,interval=interval,themask = ls_mask,force=False,obs_var='sfc_sw_down_all_mon',label='CERES down flux',shift_lon=shift_lon)

        #--- Albedo data
        alb_file_raw  = get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged.nc'
        alb_data, alb_monthly = preprocess_seasonal_data(alb_file_raw,interval=interval,themask = ls_mask,force=False,obs_var='surface_albedo_WSA',label='MODIS WSA',shift_lon=shift_lon)

        #--- get climatological mean upward flux (not for monthly data, as it is not ensured that same timeperiod is covered)
        # climatologies are already sorted with timsort()
        up_flux = alb_data.mul(obs_down)
        obs_raw_file = get_temporary_directory() + 'MODIS_CERES_surface_upward_flux.nc'
        obs_var = 'surface_shortwave_upward_flux'
        up_flux.unit = 'W/m**2'
        up_flux.save(obs_raw_file,varname=obs_var,delete=True)

        del obs_down, obs_down_monthly, alb_data, alb_monthly, up_flux

    else:
        raise ValueError, 'Invalid option for observation type in surface_upward_flux_analyis: ' + obs_type

    #/// do data preprocessing ///
    obs_up,obs_monthly = preprocess_seasonal_data(obs_raw_file,interval=interval,themask=ls_mask,force=False,obs_var=obs_var,label=obs_type,shift_lon=shift_lon)

    if GM != None:
        GM.plot(obs_monthly,linestyle='--')

    #--- initialize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        print '    SURFACE UPWARD FLUX analysis of model: ', model.name

        GP.add_model(model.name) #register model for Gleckler Plot

        if GM != None:


            if 'surface_upward_flux_org' in model.variables.keys():
                GM.plot(model.variables['surface_upward_flux_org'][2],label=model.name,mask=ls_mask) #(time,meandata)

        #--- get model data
        model_data = model.variables['surface_upward_flux']
        model_data._apply_mask(ls_mask)

        #--- use only valid albedo data (invalid values might be due to polar night effects)
        #model_data.data = np.ma.array(model_data.data,mask = ((model_data.data<0.) | (model_data.data > 1.)) )

        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != obs_up.data.shape:
            print 'Inconsistent geometries for surface_upward_flux'
            print model_data.data.shape; print obs_up.data.shape
            raise ValueError, "Invalid geometries"

        #--- generate difference map
        dmin = -20.; dmax = 20.
        f_dif  = map_difference(model_data ,obs_up,nclasses=6,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,show_zonal=True,zonal_timmean=False,vmin_zonal=0.,vmax_zonal=0.7,cticks=[0.,50.,100.,150.],cticks_diff=[-20.,-10.,0.,10.,20.])

        #seasonal map
        f_season = map_season(model_data.sub(obs_up),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu_r',show_zonal=True,zonal_timmean=True,cticks=[-20.,-10.,0.,10.,20.],nclasses=6)


        #todo hovmoeller plots !!!!


        #/// Reichler statistics ///
        Diag = Diagnostic(obs_up,model_data)
        e2   = Diag.calc_reichler_index()
        if e2 != None:
            Rplot.add(e2,model_data.label,color='red')

        #/// Gleckler plot ///
        e2a = GP.calc_index(obs_up,model_data,model,'surface_upward_flux')
        GP.add_data('surface_upward_flux',model.name,e2a,pos=1)

        #/// report results
        report.subsubsection(model.name)
        report.figure(f_season,caption='Seasonal differences')
        report.figure(f_dif,caption='Mean and relative differences')


    del obs_monthly

    f_reich = Rplot.bar(title='relative model error: surface_upward_flux')
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()















#=======================================================================
# ALBEDO -- begin
#=======================================================================

def albedo_analysis(model_list,GP=None,shift_lon=None,use_basemap=False,report=None,interval='season'):
    #todo: for a proper albedo analyis it would be useful to actually compare the all-sky albedo !

    if shift_lon == None:
        raise ValueError, 'You need to specify shift_lon option!'
    if use_basemap == None:
        raise ValueError, 'You need to specify use_basemap option!'
    if report == None:
        raise ValueError, 'You need to specify report option!'

    print
    print '************************************************************'
    print '* BEGIN ALBEDO analysis ...'
    print '************************************************************'

    report.section('Surface albedo')

    fG = plt.figure(); axg = fG.add_subplot(211); axg1 = fG.add_subplot(212)
    GM = GlobalMeanPlot(ax=axg,ax1=axg1) #global mean plot

    #- MODIS white sky albedo
    report.subsection('MODIS WSA')
    report.write('MODIS albedo is based on the MODIS white-sky albedo product. Snow covered areas remain in the data product, but all pixels flagged as invalid was discarded.')
    #albedo_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,interval=interval,obs_type='MODIS',GM=GM)
    generic_analysis(obs_dict, model_list, 'albedo', 'MODIS', GP = GP, GM = GM, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    #- CERES surface albedo from all sky fluxes
    report.subsection('CERES albedo')
    #CERES is not easily possible without pre-processing!!!!
    report.write('The CERES surface albedo is calculated as the ratio of the upward and downward surface all sky shortwave radiation fluxes based on CERES EBAF v2.6.' )
    albedo_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,interval=interval,obs_type='CERES',GM=GM)

    report.figure(fG,caption='Global means for land surface albedo')

    print '************************************************************'
    print '* END ALBEDO analysis ...'
    print '************************************************************'
    print


def albedo_analysis_plots(model_list,GP=None,shift_lon=None,use_basemap=False,report=None,interval=None,obs_type=None,GM=None):
    """
    model_list = list which contains objects of data type MODEL
    """

    vmin = 0.; vmax = 0.6

    if interval == None:
        raise ValueError, 'Interval period in albedo analyis not specified!'
    if obs_type == None:
        raise ValueError, 'Observation type for albedo was not specified!'

    print 'Doing ALBEDO analysis ...'

    #--- GlecklerPlot
    if GP == None:
        GP = GlecklerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- loading observation data
    if obs_type == 'MODIS':
        #--- load MODIS data
        #use file which was already mapped to T63; anything else takes too long!
        obs_raw_file     = get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged.nc'
        obs_var = 'surface_albedo_WSA'
        gleckler_pos = 1

    elif obs_type == 'CERES':
        #CERES EBAF ... calculate albedo from raw files
        up_file   = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres_ebaf2.6/CERES_EBAF-Surface__Ed2.6r__sfc_sw_up_all_mon__1x1__200003-201002.nc'
        down_file = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres_ebaf2.6/CERES_EBAF-Surface__Ed2.6r__sfc_sw_down_all_mon__1x1__200003-201002.nc'

        #- calculate albedo using cdos
        cdo = Cdo()
        obs_raw_file = get_temporary_directory() + os.path.basename(up_file)[:-3]+'_albedo.nc'
        cdo.div(options='-f nc', output=obs_raw_file,force=False,input=up_file + ' ' + down_file)

        obs_var = 'sfc_sw_up_all_mon'
        gleckler_pos = 2


    #todo integrate CLARA SAL DATA and CERES DATA
    else:
        raise ValueError, 'Invalid option for observation type in albedo analysis: ' + obs_type




    #/// do data preprocessing ///
    obs_alb,obs_monthly = preprocess_seasonal_data(obs_raw_file,interval=interval,themask=ls_mask,force=False,obs_var=obs_var,label=obs_type,shift_lon=shift_lon)


    if GM != None:
        GM.plot(obs_monthly,linestyle='--')

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        print '    ALBEDO analysis of model: ', model.name

        GP.add_model(model.name) #register model for Gleckler Plot

        if GM != None:
            if 'albedo_org' in model.variables.keys():
                GM.plot(model.variables['albedo_org'][2],label=model.name) #(time,meandata)

        #--- get model data
        model_data = model.variables['albedo']

        #--- use only valid albedo data (invalid values might be due to polar night effects)
        #model_data.data = np.ma.array(model_data.data,mask = ((model_data.data<0.) | (model_data.data > 1.)) )

        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != obs_alb.data.shape:
            print 'Inconsistent geometries for ALBEDO'
            print model_data.data.shape; print obs_alb.data.shape
            raise ValueError, "Invalid geometries"

        #--- generate difference map
        dmin = -0.09; dmax = 0.09
        f_dif  = map_difference(model_data ,obs_alb,nclasses=6,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,show_zonal=True,zonal_timmean=False,vmin_zonal=0.,vmax_zonal=0.7,cticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6],cticks_diff=[-0.09,-0.06,-0.03,0.,0.03,0.06,0.09])

        #seasonal map
        f_season = map_season(model_data.sub(obs_alb),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu_r',show_zonal=True,zonal_timmean=True,cticks=[-0.09,-0.06,-0.03,0.,0.03,0.06,0.09],nclasses=6)


        #todo hovmoeller plots !!!!


        #/// Reichler statistics ///
        Diag = Diagnostic(obs_alb,model_data)
        e2   = Diag.calc_reichler_index()
        Rplot.add(e2,model_data.label,color='red')

        #/// Gleckler plot ///
        e2a = GP.calc_index(obs_alb,model_data,model,'albedo')
        GP.add_data('albedo',model.name,e2a,pos=1)

        #/// report results
        report.subsubsection(model.name)
        report.figure(f_season,caption='Seasonal differences')
        report.figure(f_dif,caption='Mean and relative differences')


    del obs_monthly

    f_reich = Rplot.bar(title='relative model error: ALBEDO')
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()

#=======================================================================
# ALBEDO -- end
#=======================================================================

#=======================================================================
# TEMPERATURE -- begin
#=======================================================================


#=======================================================================
# TEMPERATURE -- end
#=======================================================================
def temperature_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):

    if report == None:
        raise ValueError, 'You need to specify report option!'

    report.section('Temperature')
    temperature_analysis_cru(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report)

def temperature_analysis_cru(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report=None):
    """
    units: K
    """
    print 'Doing Temperature analysis ...'

    vmin = 270; vmax = 320.

    #--- Gleckler plot
    model_names = []

    #--- T63 weights
    #~ t63_weights = get_T63_weights(shift_lon)

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- load CRU data

    if interval == 'season': #seasonal comparison
        #~ t2_file      = get_data_pool_directory() + 'variables/land/Ta_2m/CRUTEM3v.nc'
        #gpcp_file_std  = get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasstd.nc'

        s_start_time = '1979-01-01'; s_stop_time = '2006-12-31'
        obs_file_raw = get_data_pool_directory() + 'variables/land/Ta_2m/cru_ts_3_00.1901.2006.tmp_miss_t63.nc'

        tmp      = pyCDO(obs_file_raw,s_start_time,s_stop_time).seldate()
        tmp1     = pyCDO(tmp,s_start_time,s_stop_time).remap()
        tmp2     = pyCDO(tmp1,s_start_time,s_stop_time).seasmean()
        obs_file = pyCDO(tmp2,s_start_time,s_stop_time).yseasmean() #seasonal mean

        obs_file_std = pyCDO(tmp2,s_start_time,s_stop_time).yseasstd() #seasonal std
        obs_var = 'tmp'
    else:
        sys.exit('Unknown interval for rainfall_analyis()')

    T2 = Data(obs_file,obs_var,read=True,label='CRU',unit='K',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data,level=0)
    T2.data +=  273.15 #Kelvin
    T2_std = Data(obs_file_std,obs_var,read=True,label='CRU',unit='K',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data,level=0)
    T2.std = T2_std.data.copy(); del T2_std




    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    #--- get model field of precipitation
    for model in model_list:
        model_data = model.variables['temperature']
        GP.add_model(model.name)

        if model_data == None:
            continue

        model_names.append(model.name)

        if model_data.data.shape != T2.data.shape:
            print 'WARNING Inconsistent geometries for CRU temperature'
            print model_data.data.shape; print T2.data.shape

        dmin=-10.;dmax=10.
        f_dif = map_difference(model_data,T2,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0,5,10],cmap_difference='RdBu',show_zonal=True,zonal_timmean=False)

        #seasonal map
        f_season = map_season(model_data.sub(T2),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu',show_zonal=True,zonal_timmean=True)

        #--- calculate Reichler diagnostic for preciptation
        Diag = Diagnostic(T2,model_data); e2 = Diag.calc_reichler_index()
        Rplot.add(e2,model_data.label,color='red')

        #/// Gleckler plot ///
        e2a = GP.calc_index(T2,model_data,model,'T2')
        GP.add_data('T2',model.name,e2a,pos=1)

        #report
        report.subsection(model.name)
        report.figure(f_season,caption='Seasonal differences')
        report.figure(f_dif,caption='Mean and relative differences')

    f_reich = Rplot.bar()
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()





#=======================================================================
# RAINFALL -- begin
#=======================================================================
#def rainfall_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):
#
#    if report == None:
#        raise ValueError, 'You need to specify report option!'
#
#    report.section('Precipitation')
#
#    report.subsection('GPCP')
#    rainfall_analysis_template(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,obs_type='GPCP')
#
#    report.subsection('CRU')
#    rainfall_analysis_template(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,obs_type='CRU')
#
#    #todo add TMPA data
#
#    #~ report.subsection('GPCC')
#    #~ rainfall_analysis_template(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,obs_type='GPCC')



def xxxxxxxxxxxxxxxxxrainfall_analysis_template(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report=None,obs_type=None):
    '''
    units: mm/day
    '''
    print 'Doing rainfall analysis ...'

    if obs_type == None:
        raise ValueError, 'Can not do precipitation analysis:  missing observation type'

    vmin = 0.; vmax = 10.

    #--- Gleckler plot
    model_names = []

    #--- T63 weights
    #~ t63_weights = get_T63_weights(shift_lon)

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    if obs_type == 'GPCP':
        #--- load GPCP data
        gleckler_pos = 1
        if interval == 'season': #seasonal comparison
            obs_file      = get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasmean.nc'
            obs_file_std  = get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasstd.nc'
            obs_var = 'precip'
        else:
            sys.exit('Unknown interval for rainfall_analyis()')

        scale_data = 1.

    elif obs_type == 'CRU':
        gleckler_pos = 2
        s_start_time = '1979-01-01'; s_stop_time = '2006-12-31'

        obs_file_raw = get_data_pool_directory() + 'variables/land/precipitation/CRU/cru_ts_3_00.1901.2006.pre_miss.nc'

        tmp      = pyCDO(obs_file_raw,s_start_time,s_stop_time).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).remap()
        tmp2     = pyCDO(tmp1,s_start_time,s_stop_time).seasmean()
        obs_file = pyCDO(tmp2,s_start_time,s_stop_time).yseasmean() #seasonal mean

        obs_file_std = pyCDO(tmp2,s_start_time,s_stop_time).yseasstd() #standard deviation of seasonal mean

        scale_data = 12./365. #scale factor to scale from mm/month to mm/day #@todo: revise scale factor in precipitation analysis (not taking care yet for different month lengths!)

        obs_var = 'pre'

    elif obs_type == 'GPCC':

        raise ValueError, 'GPCC not supported yet, as no coordinates in OBS file!!' #@todo: coordinates in GPCC observational file!

        gleckler_pos = 3
        s_start_time = '1979-01-01'; s_stop_time = '2007-12-31'

        obs_file_raw = get_data_pool_directory() + 'variables/land/precipitation/GPCC/gpcc_full_vs4_1951-2007.nc'

        tmp      = pyCDO(obs_file_raw,s_start_time,s_stop_time).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).remap()
        tmp2     = pyCDO(tmp1,s_start_time,s_stop_time).seasmean()
        obs_file = pyCDO(tmp2,s_start_time,s_stop_time).yseasmean() #seasonal mean

        obs_file_std = pyCDO(tmp2,s_start_time,s_stop_time).yseasstd() #standard deviation of seasonal mean

        scale_data = 12./365. #scale factor to scale from mm/month to mm/day #@todo: revise scale factor in precipitation analysis (not taking care yet for different month lengths!)

        obs_var = 'var260'


    else:
        raise ValueError, 'UNKNOWN obs_type: ' + obs_type

    theobs = Data(obs_file,obs_var,read=True,label=obs_type,unit='mm/day',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data,scale_factor = scale_data)
    theobs_std = Data(obs_file_std,obs_var,read=True,label=obs_type,unit='mm/day',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    theobs.std = theobs_std.data.copy(); del theobs_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    #--- get model field of precipitation
    for model in model_list:
        model_data = model.variables['rain']
        GP.add_model(model.name)

        if model_data == None:
            continue

        model_names.append(model.name)

        if model_data.data.shape != theobs.data.shape:
            print 'WARNING Inconsistent geometries for GPCP'
            print model_data.data.shape; print theobs.data.shape

        dmin=-1.;dmax=1.
        f_dif = map_difference(model_data,theobs,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0,5,10],cmap_difference='RdBu',show_zonal=True,zonal_timmean=False)

        #seasonal map
        f_season = map_season(model_data.sub(theobs),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu',show_zonal=True,zonal_timmean=True)

        #--- calculate Reichler diagnostic for preciptation
        Diag = Diagnostic(theobs,model_data); e2 = Diag.calc_reichler_index()
        Rplot.add(e2,model_data.label,color='red')

        #/// Gleckler plot ///
        e2a = GP.calc_index(theobs,model_data,model,'rain')
        GP.add_data('rain',model.name,e2a,pos=gleckler_pos)

        #report
        report.subsubsection(model.name)
        report.figure(f_season,caption='Seasonal differences')
        report.figure(f_dif,caption='Mean and relative differences')

    f_reich = Rplot.bar()
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()


#=======================================================================
# RAINFALL -- end
#=======================================================================


#=======================================================================
# SIS -- begin
#=======================================================================


def sis_analysis(model_list,interval = 'season', GP=None,shift_lon=None,use_basemap=None,report=None):
    """
    main routine for SIS analysis

    calls currently 4 different analyses with different observational datasets
    """
    if shift_lon == None:
        raise ValueError, 'You need to specify shift_lon option!'
    if use_basemap == None:
        raise ValueError, 'You need to specify use_basemap option!'
    if report == None:
        raise ValueError, 'You need to specify report option!'

    vmin=0.;vmax=300;dmin=-18.;dmax = 18.

    print
    print '************************************************************'
    print '* BEGIN SIS analysis ...'
    print '************************************************************'

    report.section('Shortwave downwelling radiation (SIS)')
    fG = plt.figure(); axg = fG.add_subplot(211); axg1 = fG.add_subplot(212)
    GM = GlobalMeanPlot(ax=axg,ax1=axg1) #global mean plot

    #ISCCP
    report.subsection('ISCCP')
    generic_analysis(obs_dict, model_list, 'sis', 'ISCCP', GP = GP, GM = GM, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    #SRB
    report.subsection('SRB')
    generic_analysis(obs_dict, model_list, 'sis', 'SRB', GP = GP, GM = GM, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    #ceres
    report.subsection('CERES')
    generic_analysis(obs_dict, model_list, 'sis', 'CERES', GP = GP, GM = GM, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    #cm-saf
    report.subsection('CMSAF')
    report.write('Please note that the CMSAF analysis is limited to the Meteosat spatial domain!')
    generic_analysis(obs_dict, model_list, 'sis', 'CMSAF', GP = GP, GM = GM, report = report, use_basemap = use_basemap, shift_lon = shift_lon,interval=interval)

    report.figure(fG,caption='Global means for SIS ')

    print '************************************************************'
    print '* END SIS analysis ...'
    print '************************************************************'
    print


#-----------------------------------------------------------------------

def xxxxxxxxxxxsis_analysis_plots(model_list,interval = 'season',GP=None,GM=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.,obs_type=None,report=None):
    """
    model_list = list which contains objects of data type MODEL

    @param GM: global mean plot
    @type GM: C{GlobalMeanPlot}
    """

    print '    ... ' + obs_type

    #--- GlecklerPlot
    if GP == None:
        GP = GlecklerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)


    if obs_type == 'CERES':
        #todo EBAF data
        #raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres/T63_CERES__srbavg__surface_downwelling_shortwave_radiative_flux_in_air__1x1__2000_2004.nc'
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres_ebaf2.6/CERES_EBAF-Surface__Ed2.6r__sfc_sw_down_all_mon__1x1__200003-201002.nc'
        #y1 = '2001-01-01'; y2='2003-12-31'

        #obs_var = 'BfCER00'
        obs_var = 'sfc_sw_down_all_mon'
        gleckler_pos = 2

    else:
        print obs_type
        raise ValueError, 'Unknown observation type for SIS-analysis!'

    raise ValueError, 'This routine is outdated!'

    #/// do data preprocessing ///
    obs_sis, obs_monthly = preprocess_seasonal_data(raw_sis,interval=interval,themask = ls_mask,force=False,obs_var=obs_var,label=obs_type,shift_lon=shift_lon)

    if GM != None:
        GM.plot(obs_monthly,linestyle='--')

    #--- initialize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:
        GP.add_model(model.name) #register model name in GlecklerPlot

        if GM != None:
            if 'sis_org' in model.variables.keys():
                GM.plot(model.variables['sis_org'][2],label=model.name) #(time,meandata)

        #--- get model data
        model_data = model.variables['sis'] #model_data is a Data object!
        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != obs_sis.data.shape:
            print model_data.data.shape; print obs_sis.data.shape
            raise ValueError, 'Inconsistent geometries for SIS'


        #--- generate difference maps for each month/season

        #todo test for calculation ????

        #use welch test to calculate significant different areas
        #isdifferent , t1, t2 = welchs_approximate_ttest(model_data.n, model_data.data, model_data.std, obs_sis.n, obs_sis.data, obs_sis.std, 0.95)

        f_season1 = map_season(model_data,titlefontsize=10,cmap='jet',vmin=0.,vmax=350.,cticks=[0.,100.,200.,300.],nclasses=7,use_basemap=use_basemap)
        f_season2 = map_season(obs_sis,titlefontsize=10,cmap='jet',vmin=0.,vmax=350.,cticks=[0.,100.,200.,300.],nclasses=7,use_basemap=use_basemap)
        #f_season3 = map_season(model_data.sub(obs_sis),overlay=~isdifferent,titlefontsize=10,cmap='RdBu_r',vmin=-50.,vmax=50.,cticks=[-50.,-25.,0.,25.,50.],nclasses=8,use_basemap=use_basemap)
        f_season3 = map_season(model_data.sub(obs_sis),titlefontsize=10,cmap='RdBu_r',vmin=-50.,vmax=50.,cticks=[-50.,-25.,0.,25.,50.],nclasses=8,use_basemap=use_basemap)

        report.figure(f_season1,caption='SSI climatology of  ' + model.name)
        report.figure(f_season2,caption='SSI climatology of  ' + obs_type)
        report.figure(f_season3,caption='Difference between ' + model.name + ' and ' + obs_type + '; areas with significant differences ($p<0.05$) are shown, while areas with the same means are marked/shaded ' )
        del f_season1, f_season2, f_season3


        #--- generate hovmoeller plot ---
        if False:
            f_hov = plt.figure(figsize=(8,12))
            ax1=f_hov.add_subplot(4,1,1); ax2=f_hov.add_subplot(4,1,2)
            ax3=f_hov.add_subplot(4,1,3); ax4=f_hov.add_subplot(4,1,4)

            #hovmoeller for model
            start_time = pl.num2date(pl.datestr2num('1980-01-01')) #common period of data
            stop_time  = pl.num2date(pl.datestr2num('2012-12-31'))

            #generate a reference monthly timeseries (datetime)
            tref = rrule(MONTHLY, dtstart = start_time).between(start_time, stop_time, inc=True) #monthly timeseries

            #perform temporal subsetting and interpolation for hovmoeller plot
            tmp = model.variables['sis_org'][2]
            #i1,i2 = tmp._get_time_indices(start_time,stop_time)
            #tmp._temporal_subsetting(i1,i2)
            tmp = tmp.interp_time(pl.date2num(tref))
            print '      interpol done 1'

            hov_model = hovmoeller(num2date(tmp.time),None,rescaley=20,rescalex=20)
            hov_model.plot(climits=[0.,300.],input=tmp,xtickrotation=90,cmap='jet',ax=ax1,showcolorbar=True,showxticks=False)
            hov_model.hov = None
            hov_model.plot(climits=[-10.,10.],input=tmp.get_deseasonalized_anomaly(base='current'),xtickrotation=90,cmap='RdBu_r',ax=ax2,showcolorbar=True,showxticks=True)
            del hov_model, tmp

            #hovmoeller for observations
            tmp = obs_monthly.copy()
            #i1,i2 = tmp._get_time_indices(start_time,stop_time)
            #tmp._temporal_subsetting(i1,i2)
            tmp = tmp.interp_time(pl.date2num(tref))
            print 'interpol done 2'


            hov_obs = hovmoeller(num2date(tmp.time),None,rescaley=20,rescalex=20)
            hov_obs.plot(climits=[0.,300.],input=tmp,xtickrotation=90,cmap='jet',ax=ax3,showcolorbar=True,showxticks=False)
            hov_obs.hov = None
            hov_obs.plot(climits=[-10.,10.],input=tmp.get_deseasonalized_anomaly(base='current'),xtickrotation=90,cmap='RdBu_r',ax=ax4,showcolorbar=True)
            del hov_obs, tmp

            report.figure(f_hov,caption='Time-latitude diagram of SIS and SIS anomalies (top: ' + model.name + ', bottom: ' + obs_type + ')' )
            del f_hov

        #--- generate difference map
        f_dif  = map_difference(model_data,obs_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=6,show_zonal=True,zonal_timmean=False,cticks=[0.,50.,100.,150.,200.,250.,300.],cticks_diff=[-18.,-12.,-6.,0.,6.,12.,18.],rmin=-0.25,rmax=0.25)

        #/// Reichler statistics ///
        Diag = Diagnostic(obs_sis,model_data)
        e2   = Diag.calc_reichler_index()
        Rplot.add(e2,model_data.label,color='red')
        #print e2

        #/// Gleckler plot ///
        e2a = GP.calc_index(obs_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=gleckler_pos)

        #/// report results
        report.subsubsection(model.name)
        report.figure(f_dif,caption='Mean and relative differences ' + obs_type + ' ' + model.name)
        del f_dif

    del obs_monthly



    f_reich = Rplot.bar(title='relative model error: SIS')
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()

    sys.stdout.write('\n *** Processing finished. \n')






