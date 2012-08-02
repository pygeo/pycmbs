#!/usr/bin/env python
# -*- coding: utf-8 -*-

from utils import *

#=======================================================================
# VEGETATION COVER FRACTION -- begin
#=======================================================================


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


def albedo_analysis(model_list):
    albedo_analysis_modis(model_list)



def albedo_analysis_modis(model_list):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 0.6

    print 'Doing albedo analysis ...'

    #--- GleckerPlot
    GP = global_glecker
    if GP == None:
        GP = GleckerPlot()


    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- T63 weights
    t63_weights = get_T63_weights()

    #--- load MODIS data
    modis_file     = '/home/m300028/shared/data/SEP/variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged_2001_2010_seas_mean.nc'
    modis_file_std = '/home/m300028/shared/data/SEP/variables/land/surface_albedo/modis/with_snow/T63_MCD43C3-QC_merged_2001_2010_seas_std.nc'
    albedo=Data(modis_file,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    albedo_std=Data(modis_file_std,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    albedo.std = albedo_std.data.copy(); del albedo_std


    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        GP.add_model(model.name) #register model for Glecker Plot

        #--- get model data
        model_data = model.variables['albedo']

        if model_data.data.shape != albedo.data.shape:
            print 'Inconsistent geometries for ALBEDO'
            print model_data.data.shape; print albedo.data.shape; stop

        #--- generate difference map
        dmin = -0.1; dmax = 0.1
        dif  = map_difference(model_data ,albedo,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10)

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

        #calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(albedo,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(albedo),None,xlim=[-0.2,0.2]) #None == no area weighting performed
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #/// Reichler statistics ///
        Diag = Diagnostic(albedo,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = global_glecker.calc_index(albedo,model_data,model,'sis')
        global_glecker.add_data('albedo',model.name,e2a,pos=1)

    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()

#=======================================================================
# ALBEDO -- end
#=======================================================================


#=======================================================================
# RAINFALL -- begin
#=======================================================================


def rainfall_analysis(model_list,interval='season'):
    '''
    units: mm/day
    '''


    print 'Doing rainfall analysis ...'

    vmin = 0.; vmax = 10.

    #--- Glecker plot
    GP = global_glecker
    model_names = []


    #--- T63 weights
    t63_weights = get_T63_weights()

    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- load GPCP data

    #todo GPCP v2.2 data

    if interval == 'season': #seasonal comparison
        gpcp_file  = data_pool_directory + 'variables/land/precipitation/GPCP__V2_1dm__PRECIP__2.5x2.5__197901-200612_T63_seasmean_yseasmean.nc' #todo unit ??
        gpcp_file_std  = data_pool_directory + 'variables/land/precipitation/GPCP__V2_1dm__PRECIP__2.5x2.5__197901-200612_T63_seasmean_yseasstd.nc'
    else:
        sys.exit('Unknown interval for rainfall_analyis()')

    print 'GPCP data'
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
        dif = map_difference(model_data,gpcp,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0,5,10],cmap_difference='RdBu')

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
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(gpcp,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(gpcp),None) #None == no area weighting performed
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #--- calculate Reichler diagnostic for preciptation
        Diag = Diagnostic(gpcp,model_data); e2 = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

    Rplot.bar()


    for i in range(len(Rplot.e2_norm)):
            GP.add_data('rain',model_names[i],Rplot.e2_norm[i],pos=1)


#=======================================================================
# RAINFALL -- end
#=======================================================================


#=======================================================================
# SIS -- begin
#=======================================================================


def sis_analysis(model_list,interval = 'season', GP=None,shift_lon=None,use_basemap=None):
    '''
    main routine for SIS analysis

    calls currently 4 different analyses with different observational datasets
    '''
    if shift_lon == None:
        raise ValueError, 'You need to specify shift_lon option!'
    if use_basemap == None:
        raise ValueError, 'You need to specify use_basemap option!'

    vmin=0.;vmax=300;dmin=-20.;dmax = 20.

    print '   SIS analysis ...'

    isccp_sis_analysis(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)
    srb_sis_analysis(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)
    cmsaf_sis_analysis(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)
    ceres_sis_analysis(model_list,interval=interval,GP=GP,shift_lon=shift_lon,use_basemap=use_basemap,vmin=vmin,vmax=vmax,dmin=dmin,dmax = dmax)



def cmsaf_sis_analysis(model_list,interval = 'season',GP=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    print '      ... CM-SAF'

    #--- GleckerPlot
    #~ GP = global_glecker
    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    #--- load CMSAF-SIS data
    raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/cmsaf_sis/SISmm_all_t63.nc'
    y1 = '1984-01-01'; y2='2005-12-31'

    if interval == 'season':
        #aggregate to seasons
        cdo = pyCDO(raw_sis,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            cmsaf_sis_file = cdo.yseasmean()
            cmsaf_sis_std_file  = cdo.yseasstd()
        else:
            raise ValueError, 'Invalid interval option ', interval

    cmsaf_sis     = Data(cmsaf_sis_file,'SIS',read=True,label='cmsaf-sis',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    cmsaf_sis_std = Data(cmsaf_sis_std_file,'SIS',read=True,label='cmsaf-sis std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    cmsaf_sis.std = cmsaf_sis_std.data.copy(); del cmsaf_sis_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        #--- glecker plot
        GP.add_model(model.name)

        #--- get model data
        model_data = model.variables['sis']
        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name
            stop
            continue

        if model_data.data.shape != cmsaf_sis.data.shape:
            print 'Inconsistent geometries for SIS'
            print model_data.data.shape; print cmsaf_sis.data.shape; stop

        #--- generate difference map
        dif  = map_difference(model_data ,cmsaf_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10)

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

        #calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(cmsaf_sis,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(cmsaf_sis),None,xlim=[-0.2,0.2]) #None == no area weighting performed
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #/// Reichler statistics ///
        Diag = Diagnostic(cmsaf_sis,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(cmsaf_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=1)

    Rplot.bar()

#-----------------------------------------------------------------------

def isccp_sis_analysis(model_list,interval = 'season',GP=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.):
    '''
    model_list = list which contains objects of data type MODEL
    '''


    print '      ... ISCCP'

    #--- GleckerPlot
    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    #--- load ISCCP-SIS data
    f1 = 'T63_ISCCP__versD1__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc'
    raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/isccp/' + f1
    y1 = '1984-01-01'; y2='2005-12-31' #todo: specifiy externally

    #--- PREPROESSING ---
    if interval == 'season':
        #aggregate to seasons
        cdo = pyCDO(raw_sis,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            isccp_sis_file = cdo.yseasmean()
            isccp_sis_std_file  = cdo.yseasstd()
        else:
            raise ValueError, 'Invalid interval option ', interval

    #--- READ DATA ---
    isccp_sis     = Data(isccp_sis_file,'BfISC84',read=True,label='isccp-sis',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    isccp_sis_std = Data(isccp_sis_std_file,'BfISC84',read=True,label='isccp-sis std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    isccp_sis.std = isccp_sis_std.data.copy(); del isccp_sis_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        GP.add_model(model.name) #register model name in GleckerPlot

        #--- get model data
        model_data = model.variables['sis']
        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name
            continue

        if model_data.data.shape != isccp_sis.data.shape:
            print 'Inconsistent geometries for SIS'
            print model_data.data.shape; print isccp_sis.data.shape; stop

        #--- generate difference map
        dif  = map_difference(model_data ,isccp_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10)

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

        #calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(isccp_sis,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(isccp_sis),None,xlim=[-0.2,0.2]) #None == no area weighting performed
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #/// Reichler statistics ///
        Diag = Diagnostic(isccp_sis,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(isccp_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=3)

    Rplot.bar()

#-----------------------------------------------------------------------

def srb_sis_analysis(model_list,interval = 'season',GP=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.):
    '''
    model_list = list which contains objects of data type MODEL
    '''


    print '      ... SRB'

    #--- GleckerPlot
    #~ GP = global_glecker
    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    #--- load CMSAF-SIS data
    f1 = 'T63_SRB__vers28__surface_downwelling_shortwave_radiative_flux_in_air__1x1__all.nc'
    raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/srb/' + f1
    y1 = '1984-01-01'; y2='2005-12-31'

    if interval == 'season':
        #aggregate to seasons
        cdo = pyCDO(raw_sis,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            srb_sis_file = cdo.yseasmean()
            srb_sis_std_file  = cdo.yseasstd()
        else:
            raise ValueError, 'Invalid interval option ', interval

    srb_sis     = Data(srb_sis_file,'BfSRB84',read=True,label='srb-sis',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    srb_sis_std = Data(srb_sis_std_file,'BfSRB84',read=True,label='srb-sis std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    srb_sis.std = srb_sis_std.data.copy(); del srb_sis_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:

        #--- glecker plot
        GP.add_model(model.name)

        #--- get model data
        model_data = model.variables['sis']
        if model_data == None: #data file was not existing
            print 'Data not existing for model: ', model.name
            stop
            continue

        if model_data.data.shape != srb_sis.data.shape:
            print 'Inconsistent geometries for SIS'
            print model_data.data.shape; print srb_sis.data.shape; stop

        #--- generate difference map

        dif  = map_difference(model_data ,srb_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10)

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

        #calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(srb_sis,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(srb_sis),None,xlim=[-0.2,0.2]) #None == no area weighting performed
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #/// Reichler statistics ///
        Diag = Diagnostic(srb_sis,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(srb_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=4)

    Rplot.bar()


#-----------------------------------------------------------------------

def ceres_sis_analysis(model_list,interval = 'season',GP=None,shift_lon=None,use_basemap=False,vmin=0.,vmax=300,dmin=-20.,dmax = 20.):
    '''
    model_list = list which contains objects of data type MODEL

    @todo: implement new CERES EBAF data

    '''

    #--- GleckerPlot
    if GP == None:
        GP = GleckerPlot()


    #~ model_names = []

    print '      ... CERES'

    #--- get land sea mask
    ls_mask = get_T63_landseamask(shift_lon)

    #--- T63 weights
    t63_weights = get_T63_weights(shift_lon)

    #--- load CMSAF-SIS data

    raw_sis        = get_data_pool_directory() + 'variables/land/surface_radiation_flux_in_air/ceres/T63_CERES__srbavg__surface_downwelling_shortwave_radiative_flux_in_air__1x1__2000_2004.nc'
    y1 = '2001-01-01'; y2='2003-12-31'

    if interval == 'season':
        #aggregate to seasons
        cdo = pyCDO(raw_sis,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            ceres_sis_file = cdo.yseasmean()
            ceres_sis_std_file  = cdo.yseasstd()
        else:
            raise ValueError, 'Invalid interval option ', interval

    ceres_sis     = Data(ceres_sis_file,'BfCER00',read=True,label='ceres-sis',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    ceres_sis_std = Data(ceres_sis_std_file,'BfCER00',read=True,label='ceres-sis std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    ceres_sis.std = ceres_sis_std.data.copy(); del ceres_sis_std


    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:
        GP.add_model(model.name)

        #--- get model data
        model_data = model.variables['sis']

        if model_data == None:
            continue

        if model_data.data.shape != ceres_sis.data.shape:
            print 'Inconsistent geometries for SIS'
            print model_data.data.shape; print ceres_sis.data.shape; stop

        #--- generate difference map
        dif  = map_difference(model_data ,ceres_sis,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,nclasses=10)

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

        #calculate zonal statistics and plot
        zon1 = ZonalPlot(ax=zax1); zon1.plot(model_data,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon2 = ZonalPlot(ax=zax2); zon2.plot(ceres_sis,None,xlim=[vmin,vmax]) #None == no area weighting performed
        zon3 = ZonalPlot(ax=zax3); zon3.plot(model_data.sub(ceres_sis),None,xlim=[-0.2,0.2]) #None == no area weighting performed
        zon3.ax.plot([0.,0.],zon3.ax.get_ylim(),color='k')

        #/// Reichler statistics ///
        Diag = Diagnostic(ceres_sis,model_data)
        e2   = Diag.calc_reichler_index(t63_weights)
        Rplot.add(e2,model_data.label,color='red')

        #/// Glecker plot ///
        e2a = GP.calc_index(ceres_sis,model_data,model,'sis')
        GP.add_data('sis',model.name,e2a,pos=2)


    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()

    #~ for i in range(len(Rplot.e2_norm)):
            #~ GP.add_data('sis',model_names[i],Rplot.e2_norm[i],pos=2)

#=======================================================================
# SIS -- end
#=======================================================================










