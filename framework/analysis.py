#!/usr/bin/env python
# -*- coding: utf-8 -*-

from utils import *
from external_analysis import *

'''
@todo: implement temperature analysis
'''


#=======================================================================
# VEGETATION COVER FRACTION -- begin
#=======================================================================


def phenology_faPAR_analysis(model_list,GP=None,shift_lon=None,use_basemap=False,report=None):
    '''
    Analysis of vegetation phenology

    Reference:
    - Dahlke, Loew, Reick: ???? @todo: complete reference when paper published
    '''

    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)


    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        GP.add_model(model.name) #register model for Glecker Plot

        #/// GREENING PHASE ANALYSIS STEP 1 ... ///
        ddir = './external/phenology_benchmarking/'
        template = ddir + 'Greening_Phase_Analysis_I_CYC.m'

        tags = [{'tag':'C:\Matlab\data\Cyclopes','value':'model-save-path'},{'tag':'FFT-Mask','value':'asdasdadsadsa'} ]

        E=ExternalAnalysis('matlab -r <INPUTFILE>',template,tags,output_directory='./tmp_' + model.name.replace(' ','') + '/' ) #temporary scripts are stored in directories that have same name as model
        E.run()

        #/// GREENING PHASE ANALYSIS STEP 2 ... ///




        #--- calculate Reichler diagnostic for preciptation
        #~ Diag = Diagnostic(gpcp,model_data); e2 = Diag.calc_reichler_index(t63_weights)
        #~ Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        #~ e2a = GP.calc_index(gpcp,model_data,model,'pheno_faPAR')
        e2a = np.random.rand() #todo
        GP.add_data('pheno_faPAR',model.name,e2a,pos=1)









def grass_fraction_analysis(model_list):
    #use same analysis script as for trees, but with different argument
    tree_fraction_analysis(model_list,pft='grass')



def tree_fraction_analysis(model_list,pft='tree'):
    '''
    '''


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


        #~ print t.shape
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

        dmin=-1;dmax=1.
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
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax])  #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(hansen,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(hansen),None)  #None == no area weighting performed
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

    report.section('Surface albedo')
    albedo_analysis_plots(model_list,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report)

def albedo_analysis_plots(model_list,GP=None,shift_lon=None,use_basemap=False,report=None):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 0.6

    print 'Doing albedo analysis ...'

    #--- GleckerPlot
    #~ GP = global_glecker
    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    #--- load MODIS data
    modis_file     = get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged_2001_2010_seas_mean.nc'
    modis_file_std = get_data_pool_directory() + 'variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged_2001_2010_seas_std.nc'
    albedo=Data(modis_file,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    albedo_std=Data(modis_file_std,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    albedo.std = albedo_std.data.copy(); del albedo_std


    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        GP.add_model(model.name) #register model for Glecker Plot

        #--- get model data
        model_data = model.variables['albedo']

        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != albedo.data.shape:
            print 'Inconsistent geometries for ALBEDO'
            print model_data.data.shape; print albedo.data.shape; stop

        #--- generate difference map
        dmin = -0.1; dmax = 0.1
        f_dif  = map_difference(model_data ,albedo,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10,show_zonal=True,zonal_timmean=False)

        #seasonal map
        f_season = map_season(model_data.sub(albedo),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu_r',show_zonal=True,zonal_timmean=True)


        #/// Reichler statistics ///
        Diag = Diagnostic(albedo,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(albedo,model_data,model,'albedo')
        GP.add_data('albedo',model.name,e2a,pos=1)

        #/// report results
        report.subsection(model.name)
        report.figure(f_season,caption='Seasonal differences')
        report.figure(f_dif,caption='Mean and relative differences')



    f_reich = Rplot.bar()
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()

#=======================================================================
# ALBEDO -- end
#=======================================================================


#=======================================================================
# RAINFALL -- begin
#=======================================================================
def rainfall_analysis(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report = None):

    if report == None:
        raise ValueError, 'You need to specify report option!'

    report.section('Precipitation')
    rainfall_analysis_gpcp(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,report=report)




def rainfall_analysis_gpcp(model_list,interval='season',GP=None,shift_lon=False,use_basemap=False,report=None):
    '''
    units: mm/day
    '''
    print 'Doing rainfall analysis ...'

    vmin = 0.; vmax = 10.

    #--- Glecker plot
    model_names = []

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- load GPCP data

    if interval == 'season': #seasonal comparison
        gpcp_file      = get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasmean.nc'
        gpcp_file_std  = get_data_pool_directory() + 'variables/land/precipitation/GPCP/GPCP__V2_2dm__PRECIP__2.5x2.5__197901-201012_T63_seasmean_yseasstd.nc'
    else:
        sys.exit('Unknown interval for rainfall_analyis()')

    gpcp = Data(gpcp_file,'precip',read=True,label='GPCP',unit='mm',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    gpcp_std = Data(gpcp_file_std,'precip',read=True,label='GPCP',unit='mm',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    gpcp.std = gpcp_std.data.copy(); del gpcp_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    #--- get model field of precipitation
    for model in model_list:
        model_data = model.variables['rain']
        GP.add_model(model.name)

        if model_data == None:
            continue

        model_names.append(model.name)

        if model_data.data.shape != gpcp.data.shape:
            print 'WARNING Inconsistent geometries for GPCP'
            print model_data.data.shape; print gpcp.data.shape

        dmin=-1.;dmax=1.
        f_dif = map_difference(model_data,gpcp,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0,5,10],cmap_difference='RdBu',show_zonal=True,zonal_timmean=False)

        #seasonal map
        f_season = map_season(model_data.sub(gpcp),vmin=dmin,vmax=dmax,use_basemap=use_basemap,cmap_data='RdBu',show_zonal=True,zonal_timmean=True)

        #--- calculate Reichler diagnostic for preciptation
        Diag = Diagnostic(gpcp,model_data); e2 = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(gpcp,model_data,model,'rain')
        GP.add_data('rain',model.name,e2a,pos=1)

        #report
        report.subsection(model.name)
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

    vmin=0.;vmax=300;dmin=-20.;dmax = 20.

    print '   SIS analysis ...'

    #isccp
    report.section('Shortwave downwelling radiation - ISCCP')
    sis_analysis_plots(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='ISCCP',report=report)
    #srb
    report.section('Shortwave downwelling radiation - SRB')
    sis_analysis_plots(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='SRB',report=report)
    #ceres
    report.section('Shortwave downwelling radiation - CERES')
    sis_analysis_plots(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='CERES',report=report)
    #cm-saf
    report.section('Shortwave downwelling radiation - CMSAF')
    sis_analysis_plots(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,obs_type='CMSAF',report=report)

#-----------------------------------------------------------------------

def sis_analysis_plots(model_list,interval = 'season',GP=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.,obs_type=None,report=None):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    print '      ... ' + obs_type

    #--- GleckerPlot
    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    if obs_type == 'ISCCP':
        #--- load ISCCP-SIS data
        f1 = 'T63_ISCCP__versD1__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc'
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/isccp/' + f1
        y1 = '1984-01-01'; y2='2005-12-31' #todo: specifiy externally

        obs_var = 'BfISC84'
        glecker_pos = 3

    elif obs_type == 'SRB':
        f1 = 'T63_SRB__vers28__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc'
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/srb/' + f1
        y1 = '1984-01-01'; y2='2005-12-31'

        obs_var = 'BfSRB84'
        glecker_pos = 4

    elif obs_type == 'CERES':
        #todo EBAF data
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres/T63_CERES__srbavg__surface_downwelling_shortwave_radiative_flux_in_air__1x1__2000_2004.nc'
        y1 = '2001-01-01'; y2='2003-12-31'

        obs_var = 'BfCER00'
        glecker_pos = 2

    elif obs_type == 'CMSAF':
        raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/cmsaf_sis/SISmm_all_t63.nc'
        y1 = '1984-01-01'; y2='2005-12-31'

        obs_var = 'SIS'
        glecker_pos = 1

    else:
        print obs_type
        raise ValueError, 'Unknown observation type for SIS-analysis!'


    #--- PREPROESSING ---
    if interval == 'season':
        #aggregate to seasons
        cdo = pyCDO(raw_sis,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            obs_sis_file = cdo.yseasmean()
            obs_sis_std_file  = cdo.yseasstd()
        else:
            raise ValueError, 'Invalid interval option ', interval

    #--- READ DATA ---
    obs_sis     = Data(obs_sis_file,obs_var,read=True,label=obs_type,unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    obs_sis_std = Data(obs_sis_std_file,obs_var,read=True,label=obs_type + ' std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    obs_sis.std = obs_sis_std.data.copy(); del obs_sis_std

    #--- initialize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:
        GP.add_model(model.name) #register model name in GleckerPlot

        #--- get model data
        model_data = model.variables['sis']
        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name; continue

        if model_data.data.shape != obs_sis.data.shape:
            print model_data.data.shape; print obs_sis.data.shape
            raise ValueError, 'Inconsistent geometries for SIS'

        #--- generate difference map
        f_dif  = map_difference(model_data,obs_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10,show_zonal=True,zonal_timmean=False)

        #/// Reichler statistics ///
        Diag = Diagnostic(obs_sis,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(obs_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=glecker_pos)

        #/// report results
        report.subsection(model.name)
        report.figure(f_dif,caption='Mean and relative differences ' + obs_type + ' ' + model.name)

    f_reich = Rplot.bar()
    report.figure(f_reich,caption='Relative model performance after Reichler and Kim, 2008')
    report.newpage()

#-----------------------------------------------------------------------


#=======================================================================
# SIS -- end
#=======================================================================










