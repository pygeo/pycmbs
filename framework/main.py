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


__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"


from pyCMBS import *

#--- global variables 
data_pool_directory = '/home/m300028/shared/dev/svn/alex/sahel_albedo_jsbach/'
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
    
    def __init__(self,filename,dic_variables,name='',**kwargs):
        
        Model.__init__(self,filename,dic_variables,name='',**kwargs)
        self.get_data()
        
        

    def get_albedo_data(self):
        '''
        get albedo data for JSBACH
        
        returns Data object
        '''

        v = 'var176'
        filename = self.data_dir + 'data/model/tra0072_echam6_BOT_mm_1983-2006_albedo_JAS.nc' #todo: proper files
        ls_mask = get_T63_landseamask()
        
        albedo = Data(filename,v,read=True,
        label='MPI-ESM albedo', unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,start_time=self.start_time,stop_time=self.stop_time,
        mask=ls_mask.data.data)

        return albedo


        
        
    def get_rainfall_data(self):
        '''
        get rainfall data for JSBACH
        
        returns Data object
        '''
        
        v = 'var4'
        filename = self.data_dir + 'data/model/tra0072_echam6_BOT_mm_1983-2006_4_JAS.nc'
        ls_mask = get_T63_landseamask()
        
        rain = Data(filename,v,read=True,scale_factor = 86400.,
        label='MPI-ESM', unit = 'mm',lat_name='lat',lon_name='lon',
        shift_lon=shift_lon,start_time=self.start_time,stop_time=self.stop_time,
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

    return d

def get_T63_landseamask():
    '''
    get JSBACH T63 land sea mask as a DATA object
    '''
    ls_file = data_pool_directory + 'data/model/jsbach_T63_GR15_4tiles_1992.nc'
    ls_mask = Data(ls_file,'slm',read=True,label='T63 land-sea mask',lat_name='lat',lon_name='lon',shift_lon=shift_lon)
    msk=ls_mask.data>0.; ls_mask.data[~msk] = 0.; ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool') #convert to bool

    return ls_mask

def rainfall_analysis(model):
    print 'Doing rainfall analysis ...'
    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- load GPCP data
    gpcp_file  = data_pool_directory + 'data/gpcp/GPCP__V2_1dm__PRECIP__2.5x2.5__198301-200612_T63_JAS.nc'
    gpcp = Data(gpcp_file,'precip',read=True,label='GPCP',unit='mm',lat_name='lat',lon_name='lon',shift_lon=shift_lon,start_time=model.start_time,stop_time=model.stop_time,mask=ls_mask.data.data)

    model_data = model.variables['rainfall']

    if model_data.data.shape != gpcp.data.shape:
        print 'Inconsistent geometries for GPCP'
        print model_data.data.shape
        print gpcp.data.shape

    dmin=-1.;dmax=1.
    dif = map_difference(model_data ,gpcp,vmin=0.,vmax=10.,dmin=dmin,dmax=dmax,use_basemap=use_basemap,cticks=[0,5,10])


def albedo_analysis(model):

    print 'Doing albedo analysis ...'
    #--- get land sea mask
    ls_mask = get_T63_landseamask()

    #--- load MODIS data
    modis_file = '/home/m300028/shared/data/CMIP5/modis/new/T63_MCD43C3-QCSnow_merged_JAS_ymean.nc'
    albedo=Data(modis_file,'surface_albedo_WSA',read=True,label='albedo',unit = '-',lat_name='lat',lon_name='lon',shift_lon=shift_lon,start_time=model.start_time,stop_time=model.stop_time,mask=ls_mask.data.data)

    model_data = model.variables['albedo']

    if model_data.data.shape != albedo.data.shape:
        print 'Inconsistent geometries for GPCP'
        print model_data.data.shape
        print albedo.data.shape

    dmin = -0.01; dmax = 0.01
    dif  = map_difference(model_data ,albedo,vmin=0.,vmax=0.6,dmin=dmin,dmax=dmax,use_basemap=use_basemap)


#one class that implements models and contains routines to extract data



from pyCMBS import *



def main():


    s_start_time = '1983-01-01'
    s_stop_time  = '2005-12-31'


    start_time = plt.num2date(plt.datestr2num(s_start_time))
    stop_time  = plt.num2date(plt.datestr2num(s_stop_time ))


    #--- specify mapping of variable to analysis script name
    scripts = get_script_names()

    #--- specify variables to analyze
    variables = ['rain','albedo']



    #--- get model results (needs to be already pre-processed)
    jsbach_variables={}
    jsbach_variables.update({'rainfall' : 'get_rainfall_data()'})
    jsbach_variables.update({'albedo' : 'get_albedo_data()'})

    jsbach = JSBACH(data_pool_directory,jsbach_variables,start_time=start_time,stop_time=stop_time,name='jsbach') #model output is in kg/m**2 s --> mm
    jsbach.get_data()

    for variable in variables:
        #--- call analysis scripts for each variable
        if variable in scripts.keys():
            print 'Doing analysis for variable ... ', variable
            print scripts[variable]
            eval(scripts[variable]+'(jsbach)')
        else:
            print 'No analysis script found for variable ... ', variable
        
        


if __name__ == '__main__':
    main()




