# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""

from utils import *
from external_analysis import *
from datetime import *
from dateutil.rrule import *
from cdo import *
from matplotlib.font_manager import FontProperties
import pickle

import matplotlib.pylab as pl

matplotlib.rcParams['backend'] = 'Agg'



#///////////////////////////////////////////////////////////////////////



def preprocess_seasonal_data(raw_file, interval=None, themask=None,
                             force=False, obs_var=None, label='',
                             shift_lon=None,
                             start_date=None, stop_date=None,
                             target_grid='t63grid',
                             interpolation_method='remapcon'):
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

    @param start_date: date where the observations should start
    @type start_date: datetime.datetime

    @param stop_date: date where the observations should stop
    @type stop_date: datetime.datetime

    @param target_grid: target grid specification for remapping. This can be either a string specifying the grid following the cdo conventions or a file which contains the grid already
    @type target_grid: str

    """

    sys.stdout.write(' *** Preprocessing ' + raw_file + '\n')

    if obs_var is None:
        raise ValueError, 'Name of variable to be processed needs to be specified!'
    if shift_lon is None:
        raise ValueError, 'Lon shift parameter needs to be specified!'

    # PREPROCESSING of observational data  ---
    cdo = Cdo()

    # 1) generate monthly mean file projected to T63
    obs_mon_file = get_temporary_directory() + os.path.basename(raw_file)

    # construct string for seldate; it is assumed that it was already checked before that start_date,stop_date are valid datetime objects!

    if (start_date is not None) and (stop_date is not None):
        print 'Temporal subsetting for ' + raw_file + ' will be performed! ', start_date, stop_date
        seldate_str = ' -seldate,' + str(start_date)[0:10] + ',' + str(stop_date)[0:10]
        obs_mon_file = obs_mon_file[:-3] + '_' + str(start_date)[0:10] + '_' + str(stop_date)[0:10] + '_monmean.nc'
    else:
        seldate_str = ''
        obs_mon_file = obs_mon_file[:-3] + '_monmean.nc'

    cdo.monmean(options='-f nc', output=obs_mon_file,
                input='-' + interpolation_method + ',' + target_grid + seldate_str + ' ' + raw_file, force=force)

    # 2) generate monthly mean or seasonal mean climatology as well as standard deviation
    if interval == 'monthly':
        obs_ymonmean_file = obs_mon_file[:-3] + '_ymonmean.nc'
        obs_ymonstd_file = obs_mon_file[:-3] + '_ymonstd.nc'
        obs_ymonsum_file = obs_mon_file[:-3] + '_ymonsum.nc'
        obs_ymonN_file = obs_mon_file[:-3] + '_ymonN.nc'
        cdo.ymonmean(options='-f nc -b 32', output=obs_ymonmean_file, input=obs_mon_file, force=force)
        cdo.ymonsum(options='-f nc -b 32', output=obs_ymonsum_file, input=obs_mon_file, force=force)
        cdo.ymonstd(options='-f nc -b 32', output=obs_ymonstd_file, input=obs_mon_file, force=force)
        cdo.div(options='-f nc -b 32', output=obs_ymonN_file, input=obs_ymonsum_file + ' ' + obs_ymonmean_file,
                force=force) #number of samples
        time_cycle = 12
    elif interval == 'season':
        obs_ymonmean_file = obs_mon_file[:-3] + '_yseasmean.nc'
        obs_ymonstd_file = obs_mon_file[:-3] + '_yseasstd.nc'
        obs_ymonsum_file = obs_mon_file[:-3] + '_yseassum.nc'
        obs_ymonN_file = obs_mon_file[:-3] + '_yseasN.nc'
        cdo.yseasmean(options='-f nc -b 32', output=obs_ymonmean_file, input=obs_mon_file, force=force)
        cdo.yseassum(options='-f nc -b 32', output=obs_ymonsum_file, input=obs_mon_file, force=force)
        cdo.yseasstd(options='-f nc -b 32', output=obs_ymonstd_file, input=obs_mon_file, force=force)
        cdo.div(options='-f nc -b 32', output=obs_ymonN_file, input=obs_ymonsum_file + ' ' + obs_ymonmean_file,
                force=force)  # number of samples
        time_cycle = 4
    else:
        print interval
        raise ValueError('Unknown temporal interval. Can not perform preprocessing! ')


    #--- READ DATA ---
    obs = Data(obs_ymonmean_file, obs_var, read=True, label=label,
               lat_name='lat', lon_name='lon', shift_lon=shift_lon,
               time_cycle=time_cycle)
    obs_std = Data(obs_ymonstd_file, obs_var, read=True,
                   label=label + ' std', lat_name='lat', lon_name='lon',
                   shift_lon=shift_lon, time_cycle=time_cycle)
    obs.std = obs_std.data.copy()
    del obs_std
    obs_N = Data(obs_ymonN_file, obs_var, read=True,
                 label=label + ' N', unit='-', lat_name='lat',
                 lon_name='lon',
                 shift_lon=shift_lon, time_cycle=time_cycle)
    obs.n = obs_N.data.copy()
    del obs_N

    # sort climatology to be sure that it starts always with January
    obs.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
    obs.timsort()

    #/// read monthly data (needed for global means and hovmoeller plots) ///
    obs_monthly = Data(obs_mon_file, obs_var, read=True, label=label,
                       lat_name='lat', lon_name='lon',
                       shift_lon=shift_lon) #,mask=ls_mask.data.data)

    #try to ensure really monthly increasing time series
    if hasattr(obs_monthly, 'time_cycle'):
        if obs_monthly.time_cycle != 12:
            obs_monthly._pad_timeseries(fill_value=obs_monthly.fill_value)
    if obs_monthly.time_cycle != 12:
        raise ValueError, 'Timecycle could still not be set!!!' + obs_mon_file

    #/// center dates of months
    obs_monthly.adjust_time(day=15)
    obs.adjust_time(day=15)

    if themask is not None:
        obs._apply_mask(themask)
        obs_monthly._apply_mask(themask)

    return obs, obs_monthly

    #--- preprocessing END ---


#///////////////////////////////////////////////////////////////////////
# ANALYSIS SECTION
#///////////////////////////////////////////////////////////////////////

#-----------------------------------------------------------------------------------------------------------------------

def seaice_concentration_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False,
                                  report=None, plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='seaice_concentration', regions=regions)


def seaice_extent_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                           plot_options=None, regions=None):
    header = 'The sea ice extent is calculated based on sea-ice concentration. It is defined by all grid cells which have a sea ice concentration larger than 15 percent.'
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='seaice_extent', regions=regions, sectionheader=header)

#-----------------------------------------------------------------------------------------------------------------------

def evaporation_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                         plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='evap', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def rainfall_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                      plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='rain', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def wind_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                  plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='wind', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def twpa_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                  plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='twpa', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def wvpa_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                  plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='wvpa', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def hair_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                  plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='hair', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def late_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                  plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='late', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def budg_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                  plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='budg', regions=regions)

#-------------------------------------------------------------------------------------------------------------

def gpp_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                 plot_options=None, regions=None):
#    main_analysis(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report = report,plot_options=plot_options,actvar='gpp',regions=regions)
    beer_analysis(model_list, interval="season", GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options)

#-----------------------------------------------------------------------------------------------------------------------

def cfc_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                 plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='cfc', regions=regions)

#-----------------------------------------------------------------------------------------------------------------------

def nsJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='nsJch', regions=regions)


def acJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='acJch', regions=regions)


def cuJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='cuJch', regions=regions)


def cbJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='cbJch', regions=regions)


def stJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='stJch', regions=regions)


def ciJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='ciJch', regions=regions)


def asJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='asJch', regions=regions)


def csJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='csJch', regions=regions)


def scJch_analysis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                   plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='scJch', regions=regions)

#=======================================================================
# GENERIC - start
#=======================================================================

def generic_analysis(plot_options, model_list, obs_type, obs_name,
                     GP=None, GM=None, shift_lon=False,
                     use_basemap=False, report=None, interval=None,
                     regions=None):
    """
    function for performing common analysis actions
    it is not a parameter specific function
    use it as a template for specific analysis

    @param      plot_options: class of PlotOptions which specifies the
    options how plots shall look like!

    @param obs_type: type of observation (variable to be analyzed);
    needs to be consistent with the variables specified in main.py and
    the config file (e.g. 'sis','rain')

    @param obs_name: name of observational record as specified in the
    INI file e.g. HOAPS, CMSAF ...
    """

    #---- GENERAL CHECKS
    if interval not in ['monthly', 'season']:
        raise ValueError, 'invalid interval in generic_analysis() ' + interval

    if obs_type not in plot_options.options.keys():
        raise ValueError, 'No plot options available for the following data: ' + obs_type

    if report is None:
        raise ReportError("Report option was not enabled")

    #---- PLOT OPTIONS

    local_plot_options = plot_options.options[
        obs_type]  # gives a dictionary with all the options for the current variable
    if 'OPTIONS' not in local_plot_options.keys():
        raise ValueError, 'No OPTIONS specified for analysis of variable ' + obs_type
    if obs_name not in local_plot_options.keys():
        print 'Observational data not existing: ', obs_name, ' ... skipping analysis'
        return

    #... now everything should be fine and the plot options can be assigned locally

    #//////////////////////////////////////////////////////////////////
    # plot options which are SPECIFIC to observational data sets
    for_report = local_plot_options[obs_name]['add_to_report']  # add a certain observational dataset to report
    if not for_report:
        print '   The following data will not be included in the report as option for reporting is not set: ' + obs_name + ' ' + obs_type
        return

    obs_raw = local_plot_options[obs_name]['obs_file']
    obs_var = local_plot_options[obs_name]['obs_var']
    gleckler_pos = local_plot_options[obs_name]['gleckler_position']

    if 'start' in local_plot_options['OPTIONS'].keys():
        obs_start = local_plot_options['OPTIONS']['start']
        if isinstance(obs_start, datetime):
            pass
        else:
            try:
                obs_start = pl.num2date(pl.datestr2num(obs_start))
            except:
                print 'WARNING: INVALID start time: ', obs_start, ' SETTING to NONE!'
                obs_start = None
    else:
        obs_start = None

    if 'stop' in local_plot_options['OPTIONS'].keys():
        obs_stop = local_plot_options['OPTIONS']['stop']
        if isinstance(obs_stop, datetime):
            pass
        else:
            try:
                obs_stop = pl.num2date(pl.datestr2num(obs_stop))
            except:
                print 'WARNING: INVALID stop time: ', obs_stop, ' SETTING to NONE!'
                obs_stop = None
    else:
        obs_stop = None

    if 'scale_data' in local_plot_options[obs_name].keys():
        obs_scale_data = local_plot_options[obs_name]['scale_data']
    else:
        obs_scale_data = 1.
    if 'add_offset' in local_plot_options[obs_name].keys():
        obs_add_offset = local_plot_options[obs_name]['add_offset']
    else:
        obs_add_offset = 0.
    if 'valid_mask' in local_plot_options[obs_name].keys():
        valid_mask = local_plot_options[obs_name]['valid_mask']
    else:
        valid_mask = 'global'
    valid_mask = valid_mask.lower()

    #//////////////////////////////////////////////////////////////////
    #--- plot options which are the same for all observational datasets
    cticks = local_plot_options['OPTIONS']['cticks']
    f_mapdifference = local_plot_options['OPTIONS']['map_difference']
    f_mapseasons = local_plot_options['OPTIONS']['map_seasons']
    f_mapseason_difference = local_plot_options['OPTIONS']['map_season_difference']
    f_preprocess = local_plot_options['OPTIONS']['preprocess']
    f_reichler = local_plot_options['OPTIONS']['reichler_plot']
    f_gleckler = local_plot_options['OPTIONS']['gleckler_plot']
    f_hovmoeller = local_plot_options['OPTIONS']['hovmoeller_plot']
    f_regional_analysis = local_plot_options['OPTIONS']['regional_analysis']
    f_globalmeanplot = local_plot_options['OPTIONS']['global_mean']
    f_pattern_correlation = local_plot_options['OPTIONS']['pattern_correlation']
    interpolation_method = local_plot_options['OPTIONS']['interpolation']
    targetgrid = local_plot_options['OPTIONS']['targetgrid']
    projection = local_plot_options['OPTIONS']['projection']
    theunit = local_plot_options['OPTIONS']['units']

    if regions is None:
        f_regional_analysis = False

    if 'nclasses' in local_plot_options['OPTIONS'].keys():
        nclasses = local_plot_options['OPTIONS']['nclasses']
    else:
        nclasses = 6

    if 'vmin' in local_plot_options['OPTIONS'].keys():
        vmin = local_plot_options['OPTIONS']['vmin']
    else:
        vmin = None

    if 'vmax' in local_plot_options['OPTIONS'].keys():
        vmax = local_plot_options['OPTIONS']['vmax']
    else:
        vmax = None

    if 'dmin' in local_plot_options['OPTIONS'].keys():
        dmin = local_plot_options['OPTIONS']['dmin']
    else:
        dmin = None
    if 'dmax' in local_plot_options['OPTIONS'].keys():
        dmax = local_plot_options['OPTIONS']['dmax']
    else:
        dmax = None
    if 'cticks_diff' in local_plot_options['OPTIONS'].keys():
        cticks_diff = local_plot_options['OPTIONS']['cticks_diff']
    else:
        cticks_diff = None

    if 'cticks_rdiff' in local_plot_options['OPTIONS'].keys():
        cticks_rdiff = local_plot_options['OPTIONS']['cticks_rdiff']
    else:
        cticks_rdiff = None


    m_data_org = obs_type + '_org'  # name of original data field

    #/// land sea mask (options: land,ocean, global for parameter area)
    #mask_antarctica masks everything below 60°S.
    #here we only mask Antarctica, if only LAND points shall be used
    if valid_mask == 'land':
        mask_antarctica = True
    elif valid_mask == 'ocean':
        mask_antarctica = False
    else:
        mask_antarctica = False

    if targetgrid == 't63grid':
        ls_mask = get_T63_landseamask(shift_lon, area=valid_mask, mask_antarctica=mask_antarctica)
    else:
        ls_mask = get_generic_landseamask(shift_lon, area=valid_mask, target_grid=targetgrid,
                                          mask_antarctica=mask_antarctica)

    #####################################################################
    # OBSERVATION DATA PREPROCESSING
    #####################################################################
    #if f_preprocess == True: #always do preprocessing
    #obs_orig corresponds to the climatological mean field
    #obs_monthly is the monthly mean timeseries of the whole record (or seasonal mean)
    obs_orig, obs_monthly = preprocess_seasonal_data(obs_raw, interval=interval, themask=ls_mask,
                                                     force=False, obs_var=obs_var, label=obs_name,
                                                     shift_lon=shift_lon, start_date=obs_start, stop_date=obs_stop,
                                                     interpolation_method=interpolation_method,
                                                     target_grid=targetgrid)

    # rescale data following CF conventions
    obs_orig.mulc(obs_scale_data, copy=False)
    obs_monthly.mulc(obs_scale_data, copy=False)
    obs_orig.addc(obs_add_offset, copy=False)
    obs_monthly.addc(obs_add_offset, copy=False)

    # set unit
    obs_orig.unit = theunit
    obs_monthly.unit = theunit

    #### IDENTIFY AREAS WHERE THERE IS AT LEAST SOME VALID DATA ####
    valid_obs = ((~obs_orig.data.mask).sum(axis=0)) > 0 #find all data where there is at least SOME data
    obs_orig._apply_mask(valid_obs)
    obs_monthly._apply_mask(valid_obs)


    #obs_orig.save(filename='testfile.nc',delete=True)

    def mask_seaice_extent(x, msk):
        tmpmsk = x.data > 15.
        x.data.mask[:, :, :] = False  # set all to false to be able to put entire ocean to zero
        x.data[:, :, :] = 0.
        x.data[tmpmsk] = 1.
        x._apply_mask(msk)
        return x

    if obs_type == 'seaice_extent':
        stat_type = 'sum'  # show area weighted sum as statistic
        obs_orig = mask_seaice_extent(obs_orig, ls_mask)
        obs_monthly = mask_seaice_extent(obs_monthly, ls_mask)
        obs_orig._apply_mask(valid_obs)
        obs_monthly._apply_mask(valid_obs)
    else:
        stat_type = 'mean'


    #####################################################################
    # PLOTS
    #####################################################################

    # initialize Reichler plot
    if f_reichler == True:
        Rplot = ReichlerPlot()  # needed here, as it might include multiple model results

    if f_pattern_correlation:  # init plot for PatternCorrelation
        fig_pc = plt.figure()
        ax_pc = fig_pc.add_subplot(111)

    if f_globalmeanplot:
        if GM is None:
            fG = plt.figure();
            axg = fG.add_subplot(211);
            axg1 = fG.add_subplot(212)
            GM = GlobalMeanPlot(ax=axg, ax1=axg1, climatology=True)  # global mean plot
        else:
            if isinstance(GM, GlobalMeanPlot):
                pass
            else:
                raise ValueError('Global mean variable GM has invalid object type')
    else:
        GM = None

    if GM is not None:
        GM.plot(obs_monthly, linestyle='--', show_std=False, group='observations', stat_type=stat_type)

    if f_mapseasons == True:  #seasonal mean plot
        f_season = map_season(obs_orig, use_basemap=use_basemap, cmap_data='jet',
                              show_zonal=False, zonal_timmean=True, nclasses=nclasses,
                              vmin=vmin, vmax=vmax, cticks=cticks, proj=projection, stat_type=stat_type)
        if len(obs_orig.data) == 4:
            report.figure(f_season, caption='Seasonal mean ' + obs_name)
        else:
            report.figure(f_season, caption='Monthly mean ' + obs_name)



    ####################################################################################################################
    # MAIN LIST OVER MODELS
    ####################################################################################################################
    for model in model_list:
        sys.stdout.write('\n *** %s Analysis of model: ' % (obs_type.upper()) + model.name + "\n")

        if model.variables[obs_type] is None:
            sys.stdout.write('\n *** WARNING: No processing for %s possible for model (likely missing data!): ' % (
                obs_type) + model.name + "\n")
            continue
        else:
            model_data = model.variables[obs_type].copy()

        model_data.unit = theunit

        #todo: make this more flexible !!!
        if obs_type == 'seaice_extent':
            model_data = mask_seaice_extent(model_data, ls_mask)

        # land/sea mask
        if ls_mask is not None:
            actmask = ls_mask.data & valid_obs
        else:
            actmask = valid_obs
        model_data._apply_mask(actmask)

        GP.add_model(model._unique_name)  # register model name in GlecklerPlot

        if for_report == True:
            #/// report results
            sys.stdout.write('\n *** Making report figures. \n')
            report.subsubsection(model._unique_name)

        if GM is not None:
            if model.name == 'mean-model':
                pass
            else:
                if m_data_org in model.variables.keys():
                    if model.variables[m_data_org][2] is None:  # invalid data for mean_model --> skip
                        pass
                    else:
                        #print model.variables[m_data_org][2].label
                        tmp = model.variables[m_data_org][2]
                        tmp._apply_mask(actmask)
                        #GM.plot(model.variables[m_data_org][2],label=model._unique_name,show_std=False,group='models') #(time,meandata) replace rain_org with data_org
                        GM.plot(tmp, label=model._unique_name, show_std=False,
                                group='models')  # (time,meandata) replace rain_org with data_org
                        del tmp

        if model_data is None:
            sys.stdout.write('Data not existing for model %s' % model.name)
            continue

        if model_data.data.shape != obs_orig.data.shape:
            print 'Inconsistent geometries'  # add here parameter name
            print 'Model: ', model_data.data.shape
            print 'Observation: ', obs_orig.data.shape
            raise ValueError("Invalid geometries")

        if f_mapdifference == True:
            sys.stdout.write('\n *** Map difference plotting. \n')
            # generate difference map
            savegraphics_prefix = os.environ['PYCMBS_OUTPUTDIR'] + obs_orig._get_label() + '___' + model_data._get_label() + '.' + os.environ['PYCMBS_OUTPUTFORMAT']
            f_dif = map_difference(model_data, obs_orig,
                                   nclasses=nclasses,
                                   use_basemap=use_basemap,
                                   show_zonal=False, zonal_timmean=True,
                                   dmin=dmin, dmax=dmax, vmin=vmin,
                                   vmax=vmax, cticks=cticks,
                                   proj=projection, stat_type=stat_type,
                                   show_stat=True, drawparallels=False,
                                   colorbar_orientation='horizontal',
                                   savegraphics_prefix=savegraphics_prefix,
                                   cticks_diff=cticks_diff,
                                   cticks_rdiff=cticks_rdiff)
            report.figure(f_dif, caption='Temporal mean fields (top) and absolute and relative differences (bottom)')
            pl.close(f_dif.number)
            del f_dif

        f_compare_timeseries = False  # todo
        if f_compare_timeseries == True:
            """
            compare timeseries of two datasets
            1) apply same mask
            2) calculate regional statistics based on a predefined mask
            3) report back figure with timeseries
            """

            #--- in general it would be good to define a mask based on predefined regions.
            # in general, such a RegionalAnalysis is already implemented at 80% in pyCMBS, but not finished yet
            #I therefore propose to just use simple latitudinal bands here for the EvaClimod purposes at the moment

            themask = np.ones(model_data.timmean().shape) #generate a mask with proper geometry

            #1) lat-band 1
            tmp = (model_data.lat > 30.) #extratropic 1
            themask[tmp] = 2.
            del tmp

            tmp = (model_data.lat < -30.) #extratropic 2
            themask[tmp] = 3.;
            del tmp

            fig_comp_time = time_series_regional_analysis(model_data, obs_orig, themask)
            report.figure(fig_comp_time, caption='Here your user defined caption ... Have fun')
            pl.close(fig_comp_time.number);
            del fig_comp_time

        if f_mapseasons == True:
            sys.stdout.write('\n *** Seasonal maps plotting\n')

            # seasonal map
            f_season = map_season(model_data, use_basemap=use_basemap,
                                  cmap_data='jet',
                                  show_zonal=False, zonal_timmean=True,
                                  nclasses=nclasses,
                                  vmin=vmin, vmax=vmax, cticks=cticks,
                                  proj=projection, stat_type=stat_type,
                                  show_stat=True,
                                  drawparallels=False, titlefontsize=10)
            if len(model_data.data) == 4:
                report.figure(f_season, caption='Seasonal mean climatology for ' + model.name)
            else:
                report.figure(f_season, caption='Monthly mean climatology for ' + model.name)
            pl.close(f_season.number);
            del f_season

        if f_mapseason_difference:
            # generate seasonal plot of difference
            f_season_dif = map_season(model_data.sub(obs_orig),
                                      use_basemap=use_basemap,
                                      cmap_data='RdBu_r',
                                      show_zonal=False,
                                      zonal_timmean=True,
                                      nclasses=nclasses,
                                      vmin=dmin, vmax=dmax,
                                      proj=projection,
                                      stat_type=stat_type,
                                      show_stat=True,
                                      drawparallels=False,
                                      titlefontsize=10)

            if len(model_data.data) == 4:
                report.figure(f_season_dif,
                              caption='Seasonal mean climatology of difference between '
                                      + model.name.upper() + ' and ' + obs_orig.label.upper())
            else:
                report.figure(f_season_dif,
                              caption='Monthly mean climatology of difference between '
                                      + model.name.upper() + ' and ' + obs_orig.label.upper())
            pl.close(f_season_dif.number);
            del f_season_dif


        if f_pattern_correlation:
            # perform pattern correlation
            PC = PatternCorrelation(obs_orig, model_data, ax=ax_pc)
            xxx = PC.plot(label=model.name.upper())


        if f_hovmoeller == True:
            print('    Doing Hovmoeller plot ...')
            f_hov = plt.figure(figsize=(8, 12))
            ax1 = f_hov.add_subplot(4, 1, 1)
            ax2 = f_hov.add_subplot(4, 1, 2)
            ax3 = f_hov.add_subplot(4, 1, 3)
            ax4 = f_hov.add_subplot(4, 1, 4)

            s_start_time = '1979-01-01' #todo: use periods specified in options !!!!
            s_stop_time = '2012-12-31'
            start_time = pl.num2date(pl.datestr2num(s_start_time))
            stop_time = pl.num2date(pl.datestr2num(s_stop_time))

            #generate a reference monthly timeseries (datetime)
            tref = np.asarray(
                rrule(MONTHLY, dtstart=start_time).between(start_time, stop_time, inc=True))  # monthly timeseries

            ####################################################
            # HOVMOELLER PLOT FOR MODEL
            ####################################################

            #--- perform temporal subsetting and interpolation for hovmoeller plot
            tmp = model.variables[obs_type + '_org'][2]
            if tmp is not None:
                tmp = tmp.interp_time(pl.date2num(tref))
                if ls_mask is not None:
                    tmp._apply_mask(ls_mask)
                hov_model = hovmoeller(tmp.date, None, rescaley=20, rescalex=20)
                hov_model.plot(climits=[vmin, vmax], input=tmp, xtickrotation=90, cmap='jet', ax=ax1, showcolorbar=True,
                               showxticks=False)

                hov_model.hov = None
                hov_model.plot(climits=[dmin, dmax], input=tmp.get_deseasonalized_anomaly(base='current'),
                               xtickrotation=90, cmap='RdBu_r', ax=ax2, showcolorbar=True, showxticks=True)
                del hov_model, tmp

            ####################################################
            # HOVMOELLER PLOT FOR OBSERVATIONS
            ####################################################
            tmp = obs_monthly.copy()
            if tmp is not None:
                tmp = tmp.interp_time(pl.date2num(tref))
                if ls_mask is not None:
                    tmp._apply_mask(ls_mask)

                hov_obs = hovmoeller(tmp.date, None, rescaley=20, rescalex=20)
                hov_obs.plot(climits=[vmin, vmax], input=tmp, xtickrotation=90, cmap='jet', ax=ax3, showcolorbar=True,
                             showxticks=False)

                hov_obs.hov = None
                hov_obs.plot(climits=[dmin, dmax], input=tmp.get_deseasonalized_anomaly(base='current'),
                             xtickrotation=90, cmap='RdBu_r', ax=ax4, showcolorbar=True)
                del hov_obs, tmp

            report.figure(f_hov,
                          caption='Time-latitude diagram of SIS and SIS anomalies (top: ' + model.name
                                  + ', bottom: ' + obs_name.upper() + ')')
            pl.close(f_hov.number)
            del f_hov

        # Perform regional analysis
        if f_regional_analysis:

            #some dummy region
            #regions_for_stat = Data(None,None)
            #regions_for_stat.data = (np.random.random(obs_orig.shape[1:3])*5. + 1.).astype('int')

            region_file = local_plot_options['OPTIONS']['region_file']
            region_file_varname = local_plot_options['OPTIONS']['region_file_varname']
            # user needs to ensure that geometry matches!
            regions_for_stat = Data(region_file, region_file_varname, read=True)

            REGSTAT = RegionalAnalysis(obs_orig, model_data, regions_for_stat)
            REGSTAT.calculate()  # after that we have the regional statistics already calculated for the current model

            # save results
            sprefix = obs_type.upper() + '_' + obs_name.upper() + '_' + model._unique_name
            REGSTAT.save(prefix = sprefix, dir = report.outdir, format='pkl')
            REGSTAT.save(prefix = sprefix, dir = report.outdir, format='txt')

            del REGSTAT

        if f_pattern_correlation:
            ax_pc.set_title('Pattern correlation: ' + obs_type.upper())
            report.figure(ax_pc.figure, caption='Pattern correlation for ' + obs_type.upper())
            report.newpage()

        if f_reichler == True:
            sys.stdout.write('\n *** Computing diagnostics (Reichler index). \n')
            Diag = Diagnostic(obs_orig, model_data)
            e2 = Diag.calc_reichler_index()
            Rplot.add(e2, model_data.label, color='red')

        if f_gleckler == True:
            sys.stdout.write('\n *** Glecker plot. \n')
            e2a = GP.calc_index(obs_orig, model_data, model, obs_type)
            GP.add_data(obs_type, model._unique_name, e2a, pos=gleckler_pos)

    del obs_monthly

    if f_reichler == True:
        sys.stdout.write('\n *** Reichler plot.\n')
        f_reich = Rplot.bar(title='relative model error: %s' % obs_type.upper())
        report.figure(f_reich, caption='Relative model performance for ' + obs_type.upper())
        report.newpage()
        del Rplot

    sys.stdout.write('\n *** Processing finished. \n')

#=======================================================================
# GENERIC - end
#=======================================================================


#=======================================================================
# VEGETATION COVER FRACTION -- begin
#=======================================================================

def phenology_faPAR_analysis(model_list, GP=None, shift_lon=None, use_basemap=False, report=None):
    """
    Analysis of vegetation phenology

    @param model_list: list of model C{Data} objects
    @param GP: Gleckler plot instance
    @param shift_lon: shift data in longitude direction
    @param use_basemap: use basemap for the analysis
    @param report: instance of Report

    Reference:
    Dahlke, C., Loew, A. & Reick, C., 2012. Robust identification of global greening phase patterns from remote
           sensing vegetation products. Journal of Climate, p.120726140719003.
           Available at: http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-11-00319.1 [Accessed December 6, 2012].

    """

    #/// OPTIONS
    f_execute = True #execute scripts

    #script names
    template1 = './external/phenology_benchmarking/Greening_Phase_Analysis_I.m'
    template2 = './external/phenology_benchmarking/Greening_Phase_Analysis_II.m'
    template3 = './external/phenology_benchmarking/Greening_Phase_Analysis_III.m'
    template4 = './external/phenology_benchmarking/Greening_Phase_Analysis_IV.m'

    if GP is None:
        GP = GlecklerPlot()

    for model in model_list:
        print '*** PHENOLOGY analysis for model: ' + model.name

        fapar_data = model.variables['phenology_faPAR']  # data object for phenology variable
        snow_data = model.variables['snow']             # data object for snow variable

        #/// output directory for analysis results of current model
        outdir = os.getcwd() + '/output/GPA-01-' + model.name.replace(' ', '')

        #/// directory where model simulations are performed and templates are copied to
        rundir = './output/tmp_' + model.name.replace(' ', '') + '/'

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
        tags = [{'tag': '<STARTYEAR>', 'value': '1995'}, {'tag': '<STOPYEAR>', 'value': '2005'},
                {'tag': '<OUTDIR>', 'value': outdir},
                {'tag': '<FFTFIGNAME>', 'value': 'FFT-Mask-' + model.name.replace(' ', '')},
                {'tag': '<INPUTDATAFILE>', 'value': data_file}, {'tag': '<DATAVARNAME>', 'value': varname},
                {'tag': '<FFTMASKFILE>', 'value': outdir + '/fft_mask.mat'}]
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>', template1, tags, options=',quit',
                             output_directory=rundir) #temporary scripts are stored in directories that have same name as model
        E.run(execute=f_execute, remove_extension=True)

        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP2
        #//////////////////////////////////////
        tags.append({'tag': '<INPUTSNOWFILE>', 'value': snowfilename})
        tags.append({'tag': '<SNOWVARNAME>', 'value': snowvar})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>', template2, tags, options=',quit',
                             output_directory=rundir)
        E.run(execute=f_execute, remove_extension=True)

        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP3
        #//////////////////////////////////////
        tags.append({'tag': '<INPUTDIRECTORY>', 'value': outdir + '/'})
        tags.append(
            {'tag': '<SENSORMASKFILE>', 'value': os.getcwd() + '/external/phenology_benchmarking/sensor_mask.mat'})
        tags.append({'tag': '<RESULTDIR_AVHRR>',
                     'value': '/net/nas2/export/eo/workspace/m300028/GPA/sensors/AVH_T63'}) #@todo this needs to compe from POOL SEP !!!!
        tags.append(
            {'tag': '<RESULTDIR_SEAWIFS>', 'value': '/net/nas2/export/eo/workspace/m300028/GPA/sensors/SEA_T63'})
        tags.append(
            {'tag': '<RESULTDIR_CYCLOPES>', 'value': '/net/nas2/export/eo/workspace/m300028/GPA/sensors/CYC_T63'})
        tags.append({'tag': '<RESULTDIR_MODIS>', 'value': '/net/nas2/export/eo/workspace/m300028/GPA/sensors/MCD_T63'})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>', template3, tags, options=',quit',
                             output_directory=rundir)
        E.run(execute=f_execute, remove_extension=True)

        #//////////////////////////////////////
        # GREENING PHASE ANALYSIS - STEP4
        #//////////////////////////////////////
        tags.append({'tag': '<FIG1CAPTION>', 'value': 'Fig1_FFT_Mask_Evergreen_Overview'})
        tags.append({'tag': '<FIG2CAPTION>', 'value': 'Fig2_GPAII_Results'})
        tags.append({'tag': '<FIG3CAPTION>', 'value': 'Fig3_GPAIII_Shifts_between_Model_Sensors'})
        tags.append({'tag': '<FIG4CAPTION>', 'value': 'Fig4_Time_series_from_various_biomes'})
        tags.append({'tag': '<RUNDIR>', 'value': rundir})
        E = ExternalAnalysis('matlab -nosplash -nodesktop -r <INPUTFILE>', template4, tags, options=',quit',
                             output_directory=rundir)
        E.run(execute=f_execute, remove_extension=True)

        #/// Gleckler plot ///
        #~ e2a = GP.calc_index(gpcp,model_data,model,'pheno_faPAR')
        e2a = 0. #@todo
        GP.add_data('pheno_faPAR', model.name, e2a, pos=1)


def grass_fraction_analysis(model_list):
    #use same analysis script as for trees, but with different argument
    tree_fraction_analysis(model_list, pft='grass')


def tree_fraction_analysis(model_list, pft='tree'):
    def fraction2netcdf(pft):
        raise ValueError, 'try to avoid using Nio!'
        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/' + pft + '_05.dat'
        #save 0.5 degree data to netCDF file
        t = pl.loadtxt(filename)
        outname = filename[:-4] + '.nc'
        if os.path.exists(outname):
            os.remove(outname)
        F = Nio.open_file(outname, 'w')
        F.create_dimension('lat', t.shape[0])
        F.create_dimension('lon', t.shape[1])
        var = F.create_variable(pft + '_fraction', 'd', ('lat', 'lon'))
        lat = F.create_variable('lat', 'd', ('lat',))
        lon = F.create_variable('lon', 'd', ('lon',))
        lon.standard_name = "grid_longitude"
        lon.units = "degrees"
        lon.axis = "X"
        lat.standard_name = "grid_latitude"
        lat.units = "degrees"
        lat.axis = "Y"
        var._FillValue = -99.
        var.assign_value(t)
        lat.assign_value(
            np.linspace(90. - 0.25, -90. + 0.25, t.shape[0])) #todo: center coordinates or corner coordinates?
        lon.assign_value(np.linspace(-180. + 0.25, 180. - 0.25, t.shape[1]))

        #todo: lon 0 ... 360 ****

        F.close()

        return outname


    #//// load tree fraction observations ////
    #1) convert to netCDF
    fraction_file = fraction2netcdf(pft)
    #2) remap to T63
    y1 = '2001-01-01';
    y2 = '2001-12-31'
    cdo = pyCDO(fraction_file, y1, y2)
    t63_file = cdo.remap(method='remapcon')

    #3) load data
    ls_mask = get_T63_landseamask()
    hansen = Data(t63_file, pft + '_fraction', read=True, label='MODIS VCF ' + pft + ' fraction', unit='-',
                  lat_name='lat', lon_name='lon', shift_lon=shift_lon, mask=ls_mask.data.data)

    #~ todo now remap to T63 grid; which remapping is the most useful ???


    #//// PERFORM ANALYSIS ///
    vmin = 0.;
    vmax = 1.
    for model in model_list:

        model_data = model.variables[pft]

        if model_data.data.shape != hansen.data.shape:
            print 'WARNING Inconsistent geometries for GPCP'
            print model_data.data.shape;
            print hansen.data.shape

        dmin = -1;
        dmax = 1.
        dif = map_difference(model_data, hansen, vmin=vmin, vmax=vmax,
                             dmin=dmin, dmax=dmax, use_basemap=use_basemap, cticks=[0., 0.25, 0.5, 0.75, 1.])

        #/// ZONAL STATISTICS
        #--- append zonal plot to difference map
        ax1 = dif.get_axes()[0];
        ax2 = dif.get_axes()[1];
        ax3 = dif.get_axes()[2]

        #create net axis for zonal plot
        #http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
        divider1 = make_axes_locatable(ax1);
        zax1 = divider1.new_horizontal("50%", pad=0.05, axes_class=maxes.Axes, pack_start=True)
        dif.add_axes(zax1)

        divider2 = make_axes_locatable(ax2);
        zax2 = divider2.new_horizontal("50%", pad=0.05, axes_class=maxes.Axes, pack_start=True)
        dif.add_axes(zax2)

        divider3 = make_axes_locatable(ax3);
        zax3 = divider3.new_horizontal("50%", pad=0.05, axes_class=maxes.Axes, pack_start=True)
        dif.add_axes(zax3)

        #--- calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1);
        zon1.plot(model_data, xlim=[vmin, vmax])
        zon2 = ZonalPlot(ax=zax2);
        zon2.plot(hansen, xlim=[vmin, vmax])
        zon3 = ZonalPlot(ax=zax3);
        zon3.plot(model_data.sub(hansen))
        zon3.ax.plot([0., 0.], zon3.ax.get_ylim(), color='k')

        #--- calculate RMSE statistics etc.
        ES = CorrelationAnalysis(model_data, hansen)
        ES.do_analysis()


#=======================================================================
# VEGETATION COVER FRACTION -- end
#=======================================================================





#=======================================================================
# SURFACE UPWARD FLUX -- begin
#=======================================================================
def surface_upward_flux_analysis(model_list, interval='season', GP=None,
                                 shift_lon=False, use_basemap=False,
                                 report=None, plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon,
                  use_basemap=use_basemap, report=report,
                  plot_options=plot_options,
                  actvar='surface_upward_flux', regions=regions)

#=======================================================================
# SURFACE UPWARD FLUX -- end
#=======================================================================



#=======================================================================
# ALBEDO -- begin
#=======================================================================



def albedo_analysis_vis(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                        plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='albedo_vis', regions=regions)


def albedo_analysis_nir(model_list, interval='season', GP=None, shift_lon=False, use_basemap=False, report=None,
                        plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon, use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='albedo_nir', regions=regions)


def albedo_analysis(model_list, GP=None, shift_lon=None, use_basemap=False, report=None, interval='season',
                    plot_options=None, regions=None):
    if shift_lon is None:
        raise ValueError('You need to specify shift_lon option!')
    if use_basemap is None:
        raise ValueError('You need to specify use_basemap option!')
    if report is None:
        raise ValueError('You need to specify report option!')

    print
    print '************************************************************'
    print '* BEGIN ALBEDO analysis ...'
    print '************************************************************'

    report.section('Surface albedo')

    # Figures
    fG = plt.figure()
    axg = fG.add_subplot(211)
    axg1 = fG.add_subplot(212)

    # GlobalMean plot
    GM = GlobalMeanPlot(ax=axg, ax1=axg1, climatology=True)

    # MODIS white sky albedo
    report.subsection('MODIS WSA')
    report.write('MODIS albedo is based on the MODIS white-sky \
                    albedo product. Snow covered areas remain in the \
                    data product, but all pixels flagged as invalid \
                    was discarded.')
    generic_analysis(plot_options, model_list, 'albedo', 'MODIS',
                     GP=GP, GM=GM, report=report,
                     use_basemap=use_basemap, shift_lon=shift_lon,
                     interval=interval, regions=regions)

    # AVHRR GAC SAL black sky albedo
    report.subsection('AVHRR CLARASAL')
    report.write('AVHRR SAL is a black-sky albedo product generated in \
                    the frame of the CM-SAF')
    generic_analysis(plot_options, model_list, 'albedo', 'CLARASAL',
                     GP=GP, GM=GM, report=report,
                     use_basemap=use_basemap, shift_lon=shift_lon,
                     interval=interval, regions=regions)

    #- CERES surface albedo from all sky fluxes
    #~ report.subsection('CERES albedo')
    #~ report.write('The CERES surface albedo is calculated as the ratio of the upward and downward surface all sky shortwave radiation fluxes based on CERES EBAF v2.6.' )
    #~ generic_analysis(plot_options, model_list, 'albedo', 'CERES',
    #~ GP = GP, GM = GM, report = report,
    #~ use_basemap=use_basemap, shift_lon = shift_lon,
    #~ interval=interval, regions=regions)

    #~ albedo_analysis_plots(model_list, GP=GP, shift_lon=shift_lon,
    #~ use_basemap=use_basemap, report=report,
    #~ interval=interval, obs_type='CERES', GM=GM,
    #~ regions=regions)


    #climatological means
    fGa = GM.plot_mean_result(dt=5.,
                              colors={'observations': 'blue', 'models': 'red'})
    fGb = GM.plot_mean_result(dt=0.1,
                              colors={'observations': 'blue', 'models': 'red'},
                              plot_clim=True)

    report.figure(fG, caption='Global means for land surface albedo', bbox_inches=None)
    report.figure(fGa, caption='Global means for land surface albedo (summary)', bbox_inches=None)
    report.figure(fGb, caption='Global means for land surface albedo (climatological summary)', bbox_inches=None)

    pl.close(fGa.number)
    pl.close(fGb.number)
    pl.close(fG.number)
    del fG, fGa, fGb

    print('************************************************************')
    print('* END ALBEDO analysis ...')
    print('************************************************************')
    print('')

    del GM

#=======================================================================
# ALBEDO -- end
#=======================================================================

#=======================================================================
# TEMPERATURE -- begin
#=======================================================================

def temperature_analysis(model_list, interval='season', GP=None,
                         shift_lon=False, use_basemap=False, report=None,
                         plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon,
                  use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='temperature',
                  regions=regions)

#=======================================================================
# TEMPERATURE -- end
#=======================================================================


#=======================================================================
# MAIN analysis
#=======================================================================

def main_analysis(model_list, interval='season', GP=None, shift_lon=False,
                  use_basemap=False, report=None, plot_options=None,
                  actvar=None, regions=None, sectionheader=''):
    """
    model_list : xxxxxxxxxxxx
        xxxxxxxxxxxxxxxxxxxxxx
    interval : str
        specifies the main analysis interval ['monthly','season']
    GP : GlecklerPlot
        instance of GlecklerPlot class; default=None, then new instance
        is generated

    actvar: variable to analyze
    """

    if interval not in ['monthly', 'season']:
        raise ValueError('Invalid interval option! %s' % interval)
    if shift_lon is None:
        raise ValueError('You need to specify shift_lon option!')
    if use_basemap is None:
        raise ValueError('You need to specify use_basemap option!')
    if report is None:
        raise ValueError('You need to specify report option!')
    if plot_options is None:
        raise ValueError('No plot options are specified. No further processing possible!')
    if actvar is None:
        raise ValueError('No name for actual variable specified! Please correct!')
    if GP is None:
        raise ValueError('Gleckler plot not specified!')

    thevar = actvar
    if thevar not in plot_options.options.keys():
        raise ValueError, 'The variable is not existing in the plot_options: ' + thevar
    thelabel = plot_options.options[thevar]['OPTIONS']['label']
    thelabel = thelabel.upper()

    print('')
    print('************************************************************')
    print('* BEGIN ' + thelabel + ' analysis ...')
    print('************************************************************')

    # section header for current variable
    report.section(thelabel)
    if sectionheader != '':
        report.write(sectionheader)

    # figure
    fG = plt.figure()
    axg = fG.add_subplot(211)
    axg1 = fG.add_subplot(212)
    GM = GlobalMeanPlot(ax=axg, ax1=axg1, climatology=True)

    if thevar in plot_options.options.keys():
        #do analysis for all observational datasets specified in INI file
        for k in plot_options.options[thevar].keys():
            if k == 'OPTIONS':
                continue
            else:
                print 'IN MAIN ANALYSIS: ', thevar, k
                if plot_options.options[thevar][k]['add_to_report']:
                    report.subsection(k)
                    generic_analysis(plot_options, model_list, thevar, k, GP=GP, GM=GM, report=report,
                                     use_basemap=use_basemap, shift_lon=shift_lon, interval=interval, regions=regions)
    else:
        raise ValueError, 'Can not do analysis for ' + thelabel + ' for some reason! Check config and plot option files!'

    report.figure(fG, caption='Global means for ' + thelabel, bbox_inches=None)
    fGa = GM.plot_mean_result(dt=5., colors={'observations': 'blue', 'models': 'red'})
    fGb = GM.plot_mean_result(dt=0., colors={'observations': 'blue', 'models': 'red'}, plot_clim=True)

    report.figure(fGa,
                  caption='Global means for ' + thelabel + ' (summary); areas indicate $\pm 1 \sigma$ as derived from'
                                                           ' the ensemble of models or observations', bbox_inches=None)
    report.figure(fGb, caption='Global means for ' + thelabel + ' (summary climatology)', bbox_inches=None)

    pl.close(fG.number)
    pl.close(fGa.number)
    pl.close(fGb.number)

    del GM, fG, fGa, fGb

    print('')
    print('************************************************************')
    print('* END ' + thelabel + ' analysis ...')
    print('************************************************************')


#=======================================================================
# SIS -- begin
#=======================================================================
def sis_analysis(model_list, interval='season', GP=None,
                 shift_lon=None, use_basemap=None, report=None,
                 plot_options=None, regions=None):
    main_analysis(model_list, interval=interval, GP=GP, shift_lon=shift_lon,
                  use_basemap=use_basemap, report=report,
                  plot_options=plot_options, actvar='sis', regions=regions)

#=======================================================================
# SIS -- end
#=======================================================================


#=======================================================================
# Beer -- begin
#=======================================================================

def beer_analysis(model_list, GP=None, shift_lon=None, use_basemap=False, report=None, interval='season',
                  plot_options=None):
    if shift_lon == None:
        raise ValueError, 'You need to specify shift_lon option!'
    if use_basemap == None:
        raise ValueError, 'You need to specify use_basemap option!'
    if report == None:
        raise ValueError, 'You need to specify report option!'

    print
    print '************************************************************'
    print '* BEGIN Beer analysis ...'
    print '************************************************************'

    report.section('GPP')

    local_plot_options = plot_options.options["gpp"]
    cticks = (0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0)

    fG = figure()
    ax = fG.add_subplot(111,
                        xlabel='lat',
                        ylabel='GPP [gC m-2 a-1]',
                        xlim=[90.0, -90.0],
                        xticks=[90, 60, 30, 0, -30, -60, -90])

    for k in plot_options.options['gpp'].keys():
        if k == 'BEER':
            f_beer = figure()
            ax_Beer = f_beer.add_subplot(111)
            ls_mask = get_T63_landseamask(shift_lon)

            obs_raw = local_plot_options['BEER']['obs_file']
            obs_var = local_plot_options['BEER']['obs_var']
            GPP_BEER = Data(obs_raw, obs_var, read=True, lat_name='lat', lon_name='lon', label='GPP Beer_etal_2010',
                            mask=ls_mask.data.data, unit='gC m-2 a-1')
            Beer_zonal = GPP_BEER.get_zonal_mean(return_object=True)
            ax.plot(Beer_zonal.lat, Beer_zonal.data.data, 'k', label='Beer etal T63')
            map_plot(GPP_BEER, title='Beer etal T63', vmin=0.0, vmax=3000.0, use_basemap=use_basemap, cticks=cticks,
                     show_zonal=True, ax=ax_Beer)

            report.subsection('Beer')
            report.figure(f_beer, caption='Beer etal 2010', bbox_inches=None)

        if k == 'OPTIONS':
            for model in model_list:
                if model.name != 'mean-model':
                    f_model = figure()
                    ax_model = f_model.add_subplot(111)
                    obs_mon_file = get_temporary_directory()
                    f = obs_mon_file + model.experiment + '_' + str(model.start_time)[0:4] + '-' + str(model.stop_time)[
                                                                                                   0:4] + '_ymonmean.nc'
                    GPP_Model = Data4D(f, 'var167', read=True, lat_name='lat', lon_name='lon',
                                       label=model.name + '(' + str(model.start_time)[0:4] + '-' + str(model.stop_time)[
                                                                                                   0:4] + ')',
                                       unit='gC m-2 a-1',
                                       scale_factor=3600. * 24. * 30. / 0.083).sum_data4D()

                    GPP_Model2 = GPP_Model.timmean(return_object=True)
                    GPP_Model2.data.mask[less(GPP_Model2.data, 0.1)] = True
                    map_plot(GPP_Model2, title=model.name, vmin=0.0, vmax=3000.0, use_basemap=use_basemap,
                             cticks=cticks, show_zonal=True, ax=ax_model)
                    JSBACH_zonal = GPP_Model2.get_zonal_mean(return_object=True)
                    ax.plot(JSBACH_zonal.lat, JSBACH_zonal.data.data,
                            label=model.name + '(' + str(model.start_time)[0:4] + '-' + str(model.stop_time)[0:4] + ')')
                    report.subsection('Model (' + model.name + ')')
                    report.figure(f_model,
                                  caption=model.name + '(' + str(model.start_time)[0:4] + '-' + str(model.stop_time)[
                                                                                                0:4] + ')',
                                  bbox_inches=None)

    fontP = FontProperties()
    fontP.set_size('xx-small')
    ax.legend(prop=fontP, frameon=True, labelspacing=0.2, loc='upper left')

    report.subsection('Zonal')
    report.figure(fG, caption='Zonal mean GPP for model and Beer from ' + '(' + str(model.start_time)[0:4] + '-' + str(
        model.stop_time)[0:4] + ')', bbox_inches=None)

    print '************************************************************'
    print '* END Beer analysis ...'
    print '************************************************************'
    print

    del f_model
    del fG


def time_series_regional_analysis(x, y, msk):
    """
    This is a dummy routine to perform a comparison of timeseries on a regional basis
    The routine does
    1) apply common mask; prerequesite is that both datasets have same geometry
    2) calculate timeseries
    3) plots results
    """

    if x.shape != y.shape:
        raise ValueError, 'Datasets have not the same geometry'

    X = x.copy();
    Y = y.copy() #copy data to not work on the original datasets

    #--- mask data as invalid ---
    invalid = x.data.mask | y.data.mask #mask as invalid if one of the datasets is invalid
    X.data.mask[invalid] = True
    Y.data.mask[invalid] = True

    X.weighting_type = 'all' #whole globe for weighting
    Y.weighting_type = 'all'

    L = LinePlot(title='mytitle')

    for v in np.unique(msk): #analysis for all unique values in mask

        #--- now do analysis for each region ---

        #1) calculate mask
        m = msk == v
        tmpx = X.copy();
        tmpx._apply_mask(m)
        tmpy = Y.copy();
        tmpy._apply_mask(m)

        #2) plot results
        L.plot(tmpx, label='X region ' + str(int(v)).zfill(4))  # calculates automatically the fldmean()
        L.plot(tmpy, label='X region ' + str(int(v)).zfill(4))

        del tmpx, tmpy

    return L.figure





