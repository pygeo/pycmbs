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

#///////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////


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
    - Dahlke, Loew, Reick: ???? @todo: complete reference when paper published
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

def albedo_analysis(model_list,GP=None,shift_lon=None,use_basemap=False,report=None):
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
    report.subsection('MODIS WSA')
    albedo_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report)

    print '************************************************************'
    print '* END ALBEDO analysis ...'
    print '************************************************************'
    print


def albedo_analysis_plots(model_list,GP=None,shift_lon=None,use_basemap=False,report=None):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 0.6

    print 'Doing ALBEDO analysis ...'

    #--- GlecklerPlot
    if GP == None:
        GP = GlecklerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    #~ t63_weights = get_T63_weights(shift_lon)

    #--- load MODIS data

    modis_file     = get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged_2001_2010_seas_mean.nc'
    modis_file_std = get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged_2001_2010_seas_std.nc'
    albedo=Data(modis_file,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)

    albedo_std=Data(modis_file_std,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    albedo.std = albedo_std.data.copy(); del albedo_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        print '    ALBEDO analysis of model: ', model.name

        GP.add_model(model.name) #register model for Gleckler Plot

        #--- get model data
        model_data = model.variables['albedo']

        #--- use only valid albedo data (invalid values might be due to polar night effects)
        model_data.data = np.ma.array(model_data.data,mask = ((model_data.data<0.) | (model_data.data > 1.)) )

        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != albedo.data.shape:
            print 'Inconsistent geometries for ALBEDO'
            print model_data.data.shape; print albedo.data.shape
            raise ValueError, "Invalid geometries"

        #--- generate difference map
        dmin = -0.09; dmax = 0.09
        f_dif  = map_difference(model_data ,albedo,nclasses=6,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,show_zonal=True,zonal_timmean=False,vmin_zonal=0.,vmax_zonal=0.7,cticks=[0.,0.1,0.2,0.3,0.4,0.5,0.6],cticks_diff=[-0.09,-0.06,-0.03,0.,0.03,0.06,0.09])

        #seasonal map
        f_season = map_season(model_data.sub(albedo),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu_r',show_zonal=True,zonal_timmean=True,cticks=[-0.09,-0.06,-0.03,0.,0.03,0.06,0.09],nclasses=6)

        #/// Reichler statistics ///
        Diag = Diagnostic(albedo,model_data)
        e2   = Diag.calc_reichler_index()
        Rplot.add(e2,model_data.label,color='red')

        #/// Gleckler plot ///
        e2a = GP.calc_index(albedo,model_data,model,'albedo')
        GP.add_data('albedo',model.name,e2a,pos=1)

        #/// report results
        report.subsubsection(model.name)
        report.figure(f_season,caption='Seasonal differences')
        report.figure(f_dif,caption='Mean and relative differences')

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
    '''
    units: K
    '''
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
    T2.data = T2.data + 273.15 #Kelvin
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
def rainfall_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):

    if report == None:
        raise ValueError, 'You need to specify report option!'

    report.section('Precipitation')

    report.subsection('GPCP')
    rainfall_analysis_template(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,obs_type='GPCP')

    report.subsection('CRU')
    rainfall_analysis_template(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,obs_type='CRU')

    #~ report.subsection('GPCC')
    #~ rainfall_analysis_template(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report,obs_type='GPCC')



def rainfall_analysis_template(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report=None,obs_type=None):
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
    '''
    main routine for SIS analysis

    calls currently 4 different analyses with different observational datasets
    '''
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
    fG = plt.figure(); axg = fG.add_subplot(111)
    GM = GlobalMeanPlot(ax=axg) #global mean plot

    #isccp
    report.subsection('ISCCP')
    sis_analysis_plots(model_list,interval=interval,GP=GP,GM=GM,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='ISCCP',report=report,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)
    #srb
    report.subsection('SRB')
    sis_analysis_plots(model_list,interval=interval,GP=GP,GM=GM,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='SRB',report=report,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)
    #ceres
    report.subsection('CERES')
    sis_analysis_plots(model_list,interval=interval,GP=GP,GM=GM,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='CERES',report=report,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)
    #cm-saf
    report.subsection('CMSAF')
    sis_analysis_plots(model_list,interval=interval,GP=GP,GM=GM,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='CMSAF',report=report,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)

    report.figure(fG,caption='Global means for SIS ')

    print '************************************************************'
    print '* END SIS analysis ...'
    print '************************************************************'
    print


#-----------------------------------------------------------------------

def sis_analysis_plots(model_list,interval = 'season',GP=None,GM=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.,obs_type=None,report=None):
    '''
    model_list = list which contains objects of data type MODEL


    @param GM: global mean plot
    @type GM: C{GlobalMeanPlot}


    '''

    print '    ... ' + obs_type

    #--- GlecklerPlot
    if GP == None:
        GP = GlecklerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    #~ t63_weights = get_T63_weights(shift_lon)

    if obs_type == 'ISCCP':
        #--- load ISCCP-SIS data
        f1 = 'T63_ISCCP__versD1__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc'
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/isccp/' + f1
        y1 = '1984-01-01'; y2='2005-12-31' #todo: specifiy externally

        obs_var = 'BfISC84'
        gleckler_pos = 3

    elif obs_type == 'SRB':
        f1 = 'T63_SRB__vers28__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc'
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/srb/' + f1
        y1 = '1984-01-01'; y2='2005-12-31'

        obs_var = 'BfSRB84'
        gleckler_pos = 4

    elif obs_type == 'CERES':
        #todo EBAF data
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres/T63_CERES__srbavg__surface_downwelling_shortwave_radiative_flux_in_air__1x1__2000_2004.nc'
        y1 = '2001-01-01'; y2='2003-12-31'

        obs_var = 'BfCER00'
        gleckler_pos = 2

    elif obs_type == 'CMSAF':
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/cmsaf_sis/SISmm_all_t63.nc'
        y1 = '1984-01-01'; y2='2005-12-31'

        obs_var = 'SIS'
        gleckler_pos = 1

    else:
        print obs_type
        raise ValueError, 'Unknown observation type for SIS-analysis!'


    #--- PREPROCESSING ---
    if interval == 'season':
        #aggregate to seasons
        cdo = pyCDO(raw_sis,y1,y2) #todo: start/stop years dynamically !!!
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            cdo = pyCDO(seasfile,y1,y2)
            obs_sis_file = cdo.yseasmean()
            obs_sis_std_file  = cdo.yseasstd()
        else:
            raise ValueError, 'Invalid interval option ', interval

    #--- READ DATA ---
    obs_sis     = Data(obs_sis_file,obs_var,read=True,label=obs_type,unit = '$W m^{-2}$',lat_name='lat',lon_name='lon',shift_lon=shift_lon) #,mask=ls_mask.data.data)
    obs_sis_std = Data(obs_sis_std_file,obs_var,read=True,label=obs_type + ' std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon) #,mask=ls_mask.data.data)
    obs_sis.std = obs_sis_std.data.copy(); del obs_sis_std

    #read monthly data if global means desired
    if GM != None:
        obs_monthly = Data(raw_sis,obs_var,read=True,label=obs_type,unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon) #,mask=ls_mask.data.data)
        GM.plot(obs_monthly,linestyle='--')
        del obs_monthly



    #--- initialize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:
        GP.add_model(model.name) #register model name in GlecklerPlot

        if GM != None:
            if 'sis_org' in model.variables.keys():
                GM.plot(model.variables['sis_org'],label=model.name) #(time,meandata)

        #--- get model data
        model_data = model.variables['sis']
        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != obs_sis.data.shape:
            print model_data.data.shape; print obs_sis.data.shape
            raise ValueError, 'Inconsistent geometries for SIS'

        #--- generate difference map
        f_dif  = map_difference(model_data,obs_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=6,show_zonal=True,zonal_timmean=False,cticks=[0.,50.,100.,150.,200.,250.,300.],cticks_diff=[-18.,-12.,-6.,0.,6.,12.,18.],rmin=-0.25,rmax=0.25)

        #/// Reichler statistics ///
        Diag = Diagnostic(obs_sis,model_data)
        e2   = Diag.calc_reichler_index()
        Rplot.add(e2,model_data.label,color='red')

        #/// Gleckler plot ///
        e2a = GP.calc_index(obs_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=gleckler_pos)

        #/// report results
        report.subsubsection(model.name)
        report.figure(f_dif,caption='Mean and relative differences ' + obs_type + ' ' + model.name)

    f_reich = Rplot.bar(title='relative model error: SIS')
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()

#-----------------------------------------------------------------------


#=======================================================================
# SIS -- end
#=======================================================================










