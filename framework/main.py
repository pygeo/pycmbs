#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
USAGE:

    main.py [-h, --help]

DESCRIPTION

    programm to plot difference in rainfall between two data sets
'''

#todo
#
# implement temperature analysis

#- implement different options for temporal aggregation/plotting of map differences
#- make routines for access of JSBACH output more flexible; what about required data pre-processing?
# why is SIS model too low ???? --> wrong data! CMIP5 seems o.k.
#
# implement interface for reading CMOR data
#
# - regional analysis based on an input mask
# - correlation and RMSE analysis and Taylor plotting
#
# other colorbar for vegetation fraction analysis
# implement grass cover fraction analysis
#
# documentation!
#
# todo: interpolation method for fractional coverage !!!
# area weighting for correlation analysis

# - significance tests of differences
# - pre-processing scripts from external configuration file
# - regional subsetting options ??

# - hovmoeller plots
#
# - get temporal subset of data and then generate e.g. seasonal difference plots
#
# - seasonals vs yearly vs individual year analysis
#
__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"


from pyCMBS import *

import matplotlib.pylab as pl



#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes





#--- global variables
#data_pool_directory = '/home/m300028/shared/dev/svn/alex/sahel_albedo_jsbach/'
if 'SEP' in os.environ.keys():
    data_pool_directory = os.environ['SEP'] #get directory of pool/SEP
else:
    data_pool_directory = '/pool/SEP/'

print 'SEP directory: ' + data_pool_directory

model_directory = '/home/m300028/shared/dev/svn/alex/sahel_albedo_jsbach/'


f_fast=True
shift_lon = use_basemap = not f_fast








class Model(Data):
    '''
    This class is the main class, specifying a climate model or a particular run
    Sub-classes for particular models or experiments are herited from this class
    '''
    def __init__(self,data_dir,dic_variables,name='',**kwargs):
        '''
        constructor for Model class

        INPUT
        -----
        filename: name of the file to read data from (single file currently)
        could be improved later by more abstract class definitions

        dic_variables: dictionary specifiying variable names for a model
        e.g. 'rainfall','var4'
        '''

        #--- set a list with different datasets for different models
        self.dic_vars = dic_variables

        #--- set some metadata
        self.name = name
        self.data_dir = data_dir

        if 'start_time' in kwargs.keys():
            self.start_time = kwargs['start_time']
        else:
            self.start_time = None
        if 'stop_time' in kwargs.keys():
            self.stop_time = kwargs['stop_time']
        else:
            self.stop_time = None


    def get_data(self):
        '''
        central routine to extract data for all variables
        using functions specified in derived class
        '''

        self.variables={}
        for k in self.dic_vars.keys():
            routine = self.dic_vars[k] #get name of routine to perform data extraction
            cmd = 'dat = self.' + routine
            #~ print cmd
            #~ if hasattr(self,routine):
            if hasattr(self,routine[0:routine.index('(')]): #check if routine name is there
                exec(cmd)
                self.variables.update({ k : dat })
            else:
                print 'WARNING: unknown function to read data (skip!) ', routine
                self.variables.update({ k : None })
                #~ sys.exit()




class CMIP5Data(Model):
    def __init__(self,data_dir,model,experiment,dic_variables,name='',**kwargs):
        Model.__init__(self,None,dic_variables,name=model,**kwargs)

        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir


    def get_surface_shortwave_radiation(self):
        filename1 = self.data_dir +  self.model + '/' + 'rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        if s_start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if s_stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        tmp  = pyCDO(filename1,s_start_time,s_stop_time).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).seasmean()
        filename = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        if not os.path.exists(filename):
            return None

        sis = Data(filename,'rsds',read=True,label=self.model,unit='W/m**2',lat_name='lat',lon_name='lon',shift_lon=False)
        print 'Data read!'

        return sis




#####################################################




class JSBACH(Model):
    '''
    '''

    def __init__(self,filename,dic_variables,experiment,name='',**kwargs):

        Model.__init__(self,filename,dic_variables,name=name,**kwargs)
        self.experiment = experiment
        self.get_data()



    def get_albedo_data(self):
        '''
        get albedo data for JSBACH

        returns Data object
        '''

        v = 'var176'

        filename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1979-2006_albedo_yseasmean.nc' #todo: proper files
        ls_mask = get_T63_landseamask()

        albedo = Data(filename,v,read=True,
        label='MPI-ESM albedo ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data)

        return albedo



    def get_tree_fraction(self):
        '''
        todo implement this for data from a real run !!!
        '''

        ls_mask = get_T63_landseamask()

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_forest_shrub.nc'
        v = 'var12'
        tree = Data(filename,v,read=True,
        label='MPI-ESM tree fraction ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data,start_time = pl.num2date(pl.datestr2num('2001-01-01')),stop_time=pl.num2date(pl.datestr2num('2001-12-31')))

        return tree

    def get_grass_fraction(self):
        '''
        todo implement this for data from a real run !!!
        '''

        ls_mask = get_T63_landseamask()

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_grass_crop_pasture_2001.nc'
        v = 'var12'
        grass = Data(filename,v,read=True,
        label='MPI-ESM tree fraction ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data,start_time = pl.num2date(pl.datestr2num('2001-01-01')),stop_time=pl.num2date(pl.datestr2num('2001-12-31')) , squeeze=True  )


        return grass









    def get_surface_shortwave_radiation(self,interval = 'season'):
        '''
        get surface shortwave incoming radiation data for JSBACH

        returns Data object
        '''

        v = 'var176'

        y1 = '1979-01-01'; y2 = '2006-12-31'
        rawfilename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1979-2006_srads.nc'

        if not os.path.exists(rawfilename):
            return None


        #--- read data
        cdo = pyCDO(rawfilename,y1,y2)
        if interval == 'season':
            seasfile = cdo.seasmean(); del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile,y1,y2)
            filename = cdo.yseasmean()
        else:
            raise ValueError, 'Invalid interval option ', interval

        #--- read land-sea mask
        ls_mask = get_T63_landseamask()

        #--- read SIS data
        sis = Data(filename,v,read=True,
        label='MPI-ESM SIS ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,
        mask=ls_mask.data.data)

        return sis


    def get_rainfall_data(self,interval='season'):
        '''
        get rainfall data for JSBACH

        returns Data object
        '''

        #todo: implement preprocessing here

        v = 'var4' #todo is this really precip ???
        if interval == 'season':
            filename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1982-2006_sel_yseasmean.nc'
        else:
            raise ValueError, 'Invalid value for interval: ' + interval

        ls_mask = get_T63_landseamask()

        if not os.path.exists(filename):
            stop
            return None

        try: #todo this is silly
            rain = Data(filename,v,read=True,scale_factor = 86400.,
            label='MPI-ESM ' + self.experiment, unit = 'mm/day',lat_name='lat',lon_name='lon',
            shift_lon=shift_lon,
            mask=ls_mask.data.data)
        except:
            v='var142'
            rain = Data(filename,v,read=True,scale_factor = 86400.,
            label='MPI-ESM ' + self.experiment, unit = 'mm/day',lat_name='lat',lon_name='lon',
            shift_lon=shift_lon,
            mask=ls_mask.data.data)



        return rain



def get_script_names():
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

def get_T63_landseamask():
    '''
    get JSBACH T63 land sea mask
    the LS mask is read from the JSBACH init file

    todo put this to the JSBACH model class
    '''
    ls_file = data_pool_directory + 'variables/land/land_sea_mask/jsbach_T63_GR15_4tiles_1992.nc'
    ls_mask = Data(ls_file,'slm',read=True,label='T63 land-sea mask',lat_name='lat',lon_name='lon',shift_lon=shift_lon)
    msk=ls_mask.data>0.; ls_mask.data[~msk] = 0.; ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool') #convert to bool

    return ls_mask

def get_T63_weights():
    '''
    get JSBACH T63 cell weights

    todo put this to the JSBACH model class
    '''
    w_file = data_pool_directory + 'variables/land/land_sea_mask/t63_weights.nc'
    weight = Data(w_file,'cell_weights',read=True,label='T63 cell weights',lat_name='lat',lon_name='lon',shift_lon=shift_lon)

    return weight.data


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
        dif = map_difference(model_data,gpcp,vmin=vmin,vmax=vmax,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0,5,10])

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


######### END TREE COVER ANALYSIS ##############









def albedo_analysis(model_list):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 0.6

    print 'Doing albedo analysis ...'

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

    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()


def sis_analysis(model_list,interval = 'season', GP=None):
    cmsaf_sis_analysis(model_list,interval=interval,GP=GP)
    ceres_sis_analysis(model_list,interval=interval,GP=GP)



def cmsaf_sis_analysis(model_list,interval = 'season',GP=None):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 300

    #~ model_names = []

    print 'Doing SIS analysis ...'

    GP = global_glecker

    #--- GleckerPlot
    if GP == None:
        GP = GleckerPlot()

    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- T63 weights
    t63_weights = get_T63_weights()

    #--- load CMSAF-SIS data
    raw_sis        = data_pool_directory + 'variables/land/surface_radiation_flux_in_air/cmsaf_sis/SISmm_all_t63.nc'
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

        GP.add_model(model.name)

        #~ model_names.append(model.name)

        #~ if model.name == None:
            #~ raise ValueError, 'Model needs to have a name attribute!'
        #~ GP.add_model(model.name)

        #--- get model data
        model_data = model.variables['sis']
        if model_data == None: #data file was not existing
            continue

        if model_data.data.shape != cmsaf_sis.data.shape:
            print 'Inconsistent geometries for SIS'
            print model_data.data.shape; print cmsaf_sis.data.shape; stop

        #--- generate difference map
        dmin = -20.; dmax = 20.
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
        #~ print e2

        #/// Glecker plot ///
        e2a = global_glecker.calc_index(cmsaf_sis,model_data,model,'sis')
        global_glecker.add_data('sis',model.name,e2a,pos=1)

    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()
    #print Rplot.e2_norm #the normalized Recihler index contains temporal aggregated and relative results
    #~ for i in range(len(Rplot.e2_norm)):
        #~ GP.add_data('sis',model_names[i],Rplot.e2_norm[i],pos=1)


def ceres_sis_analysis(model_list,interval = 'season',GP=None):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    GP = global_glecker

    #--- GleckerPlot
    if GP == None:
        GP = GleckerPlot()

    vmin = 0.; vmax = 300

    #~ model_names = []

    print 'Doing CERES SIS analysis ...'

    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- T63 weights
    t63_weights = get_T63_weights()

    #--- load CMSAF-SIS data

    raw_sis        = data_pool_directory + 'variables/land/surface_radiation_flux_in_air/ceres/T63_CERES__srbavg__surface_downwelling_shortwave_radiative_flux_in_air__1x1__2000_2004.nc'
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

        #~ model_names.append(model.name)


        #--- get model data
        model_data = model.variables['sis']

        if model_data == None:
            continue

        if model_data.data.shape != ceres_sis.data.shape:
            print 'Inconsistent geometries for SIS'
            print model_data.data.shape; print ceres_sis.data.shape; stop

        #--- generate difference map
        dmin = -20.; dmax = 20.
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
        e2a = global_glecker.calc_index(ceres_sis,model_data,model,'sis')
        global_glecker.add_data('sis',model.name,e2a,pos=2)


    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()

    #~ for i in range(len(Rplot.e2_norm)):
            #~ GP.add_data('sis',model_names[i],Rplot.e2_norm[i],pos=2)











#####################################################################
#####################################################################
#####################################################################









#one class that implements models and contains routines to extract data



from pyCMBS import *



#~ def main():

data_dir = '/home/m300028/shared/data/CMIP5/EvaCliMod/rsds/' #CMIP5 data directory


pl.close('all')

global s_start_time
global s_stop_time

s_start_time = '1983-01-01' #todo where is this used ?
s_stop_time  = '2005-12-31'




#--- specify variables to analyze
#~ variables = ['rain','albedo','sis']
#variables = ['tree','albedo']
variables = ['sis'] #sis

#--- specify mapping of variable to analysis script name
scripts = get_script_names()


#################################################

#-----
start_time = plt.num2date(plt.datestr2num(s_start_time))
stop_time  = plt.num2date(plt.datestr2num(s_stop_time ))

#--- get model results (needs to be already pre-processed)
jsbach_variables={}
hlp={}
hlp.update({'rain' : 'get_rainfall_data()'})
hlp.update({'albedo' : 'get_albedo_data()'})
hlp.update({'sis' : 'get_surface_shortwave_radiation()'})
hlp.update({'tree' : 'get_tree_fraction()'})
hlp.update({'grass' : 'get_grass_fraction()'})


cmip_model_list = ['CSIRO-Mk3-6-0','MPI-ESM-LR','MPI-ESM-MR','HadGEM2-A','AGCM3-2H','AGCM3-2S','bcc-csm1-1','CGCM3','CNRM-CM5','GFDL-HIRAM-C180','GFDL-HIRAM-C360','GISS-E2-R','inmcm4','IPSL-CM5A-LR','MIROC5','NorESM1-M']


cmip_model_list = ['CSIRO-Mk3-6-0','MPI-ESM-LR','MIROC5']
cmip_models = []




for k in hlp.keys(): #only use the variables that should be analyzed!
    print k, hlp[k]
    if k in variables:
        jsbach_variables.update({k:hlp[k]})

experiment = 'amip'
cmip_cnt = 1
for model in cmip_model_list:
    cmip = CMIP5Data(data_dir,model,experiment,jsbach_variables,unit='W/m**2',lat_name='lat',lon_name='lon',label=model)
    cmip.get_data()
    cmd = 'cmip' + str(cmip_cnt).zfill(4) + ' = ' + 'cmip.copy(); del cmip'
    exec(cmd) #store copy of cmip5 model in separate variable
    cmip_models.append('cmip' + str(cmip_cnt).zfill(4))


    cmip_cnt += 1

print str(cmip_models)





# TODO CMIP5:
# - seldate appropriately
# - check timestamp!
# seasmean calculations automatically!!

#todo check datetime; something is wrong! as data starts in Dcember 1978!
# do maps per season ???

# TODO: zonal means ??
# area weigting of zonal means and also area weighting of RMS errors etc.

#
# TODO seasmean yseasmean !!!!




jsbach72 = JSBACH(model_directory,jsbach_variables,'tra0072',start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
jsbach72.get_data()

#~ jsbach73 = JSBACH(model_directory,jsbach_variables,'tra0073',start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
#~ jsbach73.get_data()

#~ jsbach74 = JSBACH(model_directory,jsbach_variables,'tra0074',start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
#~ jsbach74.get_data()


skeys = scripts.keys()

#generate a global variable for glecker plot!
global global_glecker
global_glecker = GleckerPlot()

for variable in variables:

    global_glecker.add_variable(variable) #register current variable in Glecker Plot

    #--- call analysis scripts for each variable
    for k in range(len(skeys)):
        if variable == skeys[k]:

            print 'Doing analysis for variable ... ', variable
            print scripts[variable]
            #~ eval(scripts[variable]+'([jsbach72])') #here one can put a multitude of model output for comparison in the end
            #~ eval(scripts[variable]+'([cmip,jsbach72])') #here one can put a multitude of model output for comparison in the end
            model_list = str(cmip_models).replace("'","") #cmip5 model list
            #~ model_list = model_list.replace(']',', ') + 'jsbach72' + ']'
            eval(scripts[variable]+'(' + model_list + ')') #here one can put a multitude of model output for comparison in the end

#~ if __name__ == '__main__':
    #~ main()

global_glecker.plot(vmin=-0.8,vmax=0.8,nclasses=15)

