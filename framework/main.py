#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
USAGE:

    main.py [-h, --help]

DESCRIPTION
    
    programm to plot difference in rainfall between two data sets
'''

#todo
#- implement different options for temporal aggregation/plotting of map differences
#- make routines for access of JSBACH output more flexible; what about required data pre-processing?
# - units precipitation analysis

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
            print cmd
            #~ if hasattr(self,routine):
            if hasattr(self,routine[0:routine.index('(')]): #check if routine name is there
                exec(cmd)
                self.variables.update({ k : dat })
            else:
                print 'WARNING: unknown function to read data (skip!) ', routine
                sys.exit()
           
            
            


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



    def get_surface_shortwave_radiation(self,interval = 'season'):
        '''
        get surface shortwave incoming radiation data for JSBACH
        
        returns Data object
        '''

        v = 'var176'
        
        y1 = '1979-01-01'; y2 = '2006-12-31'
        rawfilename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1979-2006_srads.nc' 
        
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
            filename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1979-2006_precip_yseasmean.nc'
        else:
            raise ValueError, 'Invalid value for interval: ' + interval
            
        ls_mask = get_T63_landseamask()
        
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
    d.update({'sis':'cmsaf_sis_analysis'})
    d.update({'sis':'ceres_sis_analysis'})

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
        
        model_data = model.variables['rainfall']

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
    





def cmsaf_sis_analysis(model_list,interval = 'season'):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 300

    print 'Doing SIS analysis ...'
    
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
        #--- get model data
        model_data = model.variables['sis']

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
        
    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()
    


def ceres_sis_analysis(model_list,interval = 'season'):
    '''
    model_list = list which contains objects of data type MODEL
    '''

    vmin = 0.; vmax = 300

    print 'Doing CERES SIS analysis ...'
    
    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- T63 weights
    t63_weights = get_T63_weights()

    #--- load CMSAF-SIS data
    raw_sis        = data_pool_directory + 'variables/land/surface_radiation_flux_in_air/ceres/monthly2/T63_CERES__srbavg__surface_downwelling_shortwave_radiative_flux_in_air__1x1__2000mm-2003mm.nc'
    y1 = '2000-01-01'; y2='2004-12-31'
    
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
    
    ceres_sis     = Data(ceres_sis_file,'BfCER4e',read=True,label='ceres-sis',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    ceres_sis_std = Data(ceres_sis_std_file,'BfCER4e',read=True,label='ceres-sis std',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,mask=ls_mask.data.data)
    ceres_sis.std = ceres_sis_std.data.copy(); del ceres_sis_std

    #--- initailize Reichler plot
    Rplot = ReichlerPlot() #needed here, as it might include multiple model results

    for model in model_list:
        #--- get model data
        model_data = model.variables['sis']

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
        
    #~ Rplot.simple_plot()
    #~ Rplot.circle_plot()
    Rplot.bar()











#####################################################################
#####################################################################
#####################################################################









#one class that implements models and contains routines to extract data



from pyCMBS import *



def main():
    
    pl.close('all')


    s_start_time = '1983-01-01' #todo where is this used ?
    s_stop_time  = '2005-12-31'


    start_time = plt.num2date(plt.datestr2num(s_start_time))
    stop_time  = plt.num2date(plt.datestr2num(s_stop_time ))


    #--- specify mapping of variable to analysis script name
    scripts = get_script_names()

    #--- specify variables to analyze
    #~ variables = ['rain','albedo']
    variables = ['sis']



    #--- get model results (needs to be already pre-processed)
    jsbach_variables={}
    jsbach_variables.update({'rainfall' : 'get_rainfall_data()'})
    jsbach_variables.update({'albedo' : 'get_albedo_data()'})
    jsbach_variables.update({'sis' : 'get_surface_shortwave_radiation()'})

    jsbach72 = JSBACH(model_directory,jsbach_variables,'tra0072',start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
    jsbach72.get_data()
    
    jsbach73 = JSBACH(model_directory,jsbach_variables,'tra0073',start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
    jsbach73.get_data()
#~ 
    jsbach74 = JSBACH(model_directory,jsbach_variables,'tra0074',start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
    jsbach74.get_data()

    skeys = scripts.keys()
    print skeys
    print variables
    stop
    for variable in variables:
        #--- call analysis scripts for each variable
        for k in range(len(skeys)):
            if variable == skeys[k]:
                
                print 'Doing analysis for variable ... ', variable
                print scripts[variable]
                #~ eval(scripts[variable]+'([jsbach72])') #here one can put a multitude of model output for comparison in the end
                eval(scripts[variable]+'([jsbach72,jsbach73,jsbach74])') #here one can put a multitude of model output for comparison in the end


if __name__ == '__main__':
    main()




