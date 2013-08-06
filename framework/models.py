#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012-2013 Alexander Loew, alexander.loew@zmaw.de
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



from pyCMBS import *
from utils import *
from cdo import *
import tempfile as tempfile
import copy

'''
@todo: implement reading of air temperature fields
'''




class Model(Data):
    """
    This class is the main class, specifying a climate model or a particular run
    Sub-classes for particular models or experiments are herited from this class
    """
    def __init__(self,data_dir,dic_variables,name='',intervals=None,**kwargs):
        """
        constructor for Model class

        @param intervals: a dictionary from configuration, that specifies the temporal interval to be used within each analyis
        @type intervals: dict



        INPUT
        -----
        filename: name of the file to read data from (single file currently)
        could be improved later by more abstract class definitions

        dic_variables: dictionary specifiying variable names for a model
        e.g. 'rainfall','var4'
        """

        #--- check
        if intervals == None:
            raise ValueError, 'Invalid intervals for Model data: needs specification!'
        #--- set a list with different datasets for different models
        self.dic_vars = dic_variables
        self.intervals = intervals

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
        """
        central routine to extract data for all variables
        using functions specified in derived class
        """

        self.variables={}
        for k in self.dic_vars.keys():
            self._actplot_options = self.plot_options.options[k]['OPTIONS'] #set variable specific options (needed for interpolation when reading the data)

            routine  = self.dic_vars [k] #get name of routine to perform data extraction
            interval = self.intervals[k]
            cmd = 'dat = self.' + routine

            if hasattr(self,routine[0:routine.index('(')]): #check if routine name is there
                exec(cmd)

                #--- if a tuple is returned, then it is the data + a tuple for the original global mean field
                if 'tuple' in str(type(dat)):
                    self.variables.update({ k : dat[0] }) #update field with data
                    self.variables.update({ k + '_org' : dat[1]}) #(time, meanfield, originalfield)
                else:
                    self.variables.update({ k : dat }) #update field with data

            else:
                print 'WARNING: unknown function to read data (skip!), variable: ', k
                self.variables.update({ k : None })








class CMIP5Data(Model):
    """
    Class for CMIP5 model simulations. This class is derived from C{Model}.
    """
    def __init__(self,data_dir,model,experiment,dic_variables,name='',shift_lon=False,**kwargs):
        """

        @param data_dir: directory that specifies the root directory where the data is located
        @param model: TBD tood
        @param experiment: specifies the ID of the experiment (str)
        @param dic_variables:
        @param name: TBD todo
        @param shift_lon: specifies if longitudes of data need to be shifted
        @param kwargs: other keyword arguments
        @return:
        """
        Model.__init__(self,None,dic_variables,name=model,shift_lon=shift_lon,**kwargs)

        self.model      = model; self.experiment = experiment
        self.data_dir   = data_dir; self.shift_lon  = shift_lon
        self.type       = 'CMIP5'

        self._unique_name = self._get_unique_name()



    def _get_unique_name(self):
        """
        get unique name from model and experiment
        @return: string with unique combination of models and experiment
        """
        return self.model.replace(' ','') + '-' + self.experiment.replace(' ','')

#-----------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------

    def get_model_data_generic(self,interval='season', **kwargs):
        """
        unique parameters are:
            filename - file basename
            variable - name of the variable as the short_name in the netcdf file


            kwargs is a dictionary with keys for each model. Then a dictionary with properties follows

        """

        if not self.type in kwargs.keys():
            print 'WARNING: it is not possible to get data using generic function, as method missing: ', self.type, kwargs.keys()
            return None

        locdict = kwargs[self.type]

        # read settings and details from the keyword arguments
        # no defaults; everything should be explicitely specified in either the config file or the dictionaries
        varname  = locdict.pop('variable')
        units    = locdict.pop('unit', 'Crazy Unit')
        #interval = kwargs.pop('interval') #, 'season') #does not make sense to specifiy a default value as this option is specified by configuration file!

        lat_name     = locdict.pop('lat_name', 'lat')
        lon_name     = locdict.pop('lon_name', 'lon')
        model_suffix = locdict.pop('model_suffix')
        model_prefix = locdict.pop('model_prefix')
        file_format  = locdict.pop('file_format')
        scf = locdict.pop('scale_factor')
        valid_mask    = locdict.pop('valid_mask')
        custom_path  = locdict.pop('custom_path', None)
        thelevel  = locdict.pop('level', None)

        target_grid   = self._actplot_options['targetgrid']
        interpolation = self._actplot_options['interpolation']


        if custom_path == None:
            filename1 = ("%s%s/merged/%s_%s_%s_%s_%s.%s" %
                        (self.data_dir, varname, varname, model_prefix, self.model, self.experiment, model_suffix, file_format))
        else:
            if self.type == 'CMIP5':
                filename1 = ("%s/%s_%s_%s_%s_%s.%s" %
                        (custom_path, varname, model_prefix, self.model, self.experiment, model_suffix, file_format))
            elif self.type == 'CMIP3':
                filename1 = ("%s/%s_%s_%s_%s.%s" %
                        (custom_path, self.experiment, self.model, varname,  model_suffix, file_format))
            else:
                print self.type
                raise ValueError, 'Can not generate filename: invalid model type!'



        force_calc = False

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        #/// PREPROCESSING ///
        cdo = Cdo()
        s_start_time = str(self.start_time)[0:10]
        s_stop_time  = str(self.stop_time )[0:10]

        #1) select timeperiod and generate monthly mean file
        if target_grid == 't63grid':
            gridtok = 'T63'
        else:
            gridtok = 'SPECIAL_GRID'

        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_' + gridtok + '_monmean.nc' #target filename
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)

        sys.stdout.write('\n *** Model file monthly: %s\n' % file_monthly)

        if not os.path.exists(filename1):
            print 'WARNING: File not existing: ' + filename1
            return None

        cdo.monmean(options='-f nc',output=file_monthly,input = '-' + interpolation + ',' + target_grid + ' -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1, force=force_calc)


        sys.stdout.write('\n *** Reading model data... \n')
        sys.stdout.write('     Interval: ' + interval + '\n')


        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            mdata_clim_file     = file_monthly[:-3] + '_ymonmean.nc'
            mdata_sum_file      = file_monthly[:-3] + '_ymonsum.nc'
            mdata_N_file        = file_monthly[:-3] + '_ymonN.nc'
            mdata_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32',output = mdata_clim_file,input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32',output = mdata_sum_file,input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32',output = mdata_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc',output = mdata_N_file,input=mdata_sum_file + ' ' + mdata_clim_file, force=force_calc) #number of samples
        elif interval == 'season':
            mdata_clim_file     = file_monthly[:-3] + '_yseasmean.nc'
            mdata_sum_file      = file_monthly[:-3] + '_yseassum.nc'
            mdata_N_file        = file_monthly[:-3] + '_yseasN.nc'
            mdata_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32',output = mdata_clim_file,input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32',output = mdata_sum_file,input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32',output = mdata_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32',output = mdata_N_file,input=mdata_sum_file + ' ' + mdata_clim_file, force=force_calc) #number of samples
        else:
            raise ValueError, 'Unknown temporal interval. Can not perform preprocessing! '

        if not os.path.exists(mdata_clim_file):
            return None

        #3) read data
        if interval == 'monthly':
            thetime_cylce = 12
        elif interval == 'season':
            thetime_cylce = 4
        else:
            print interval
            raise ValueError, 'Unsupported interval!'
        mdata = Data(mdata_clim_file,varname,read=True,label=self.model,unit=units,lat_name=lat_name,lon_name=lon_name,shift_lon=False, scale_factor=scf,level=thelevel,time_cycle=thetime_cylce)
        mdata_std = Data(mdata_clim_std_file,varname,read=True,label=self.model+ ' std',unit='-',lat_name=lat_name,lon_name=lon_name,shift_lon=False,level=thelevel,time_cycle=thetime_cylce)
        mdata.std = mdata_std.data.copy(); del mdata_std
        mdata_N = Data(mdata_N_file,varname,read=True,label=self.model+ ' std',unit='-',lat_name=lat_name,lon_name=lon_name,shift_lon=False, scale_factor=scf,level=thelevel)
        mdata.n = mdata_N.data.copy(); del mdata_N

        #ensure that climatology always starts with J  anuary, therefore set date and then sort
        mdata.adjust_time(year=1700,day=15) #set arbitrary time for climatology
        mdata.timsort()

        #4) read monthly data
        mdata_all = Data(file_monthly,varname,read=True,label=self.model,unit=units,lat_name=lat_name,lon_name=lon_name,shift_lon=False,time_cycle=12,scale_factor=scf,level=thelevel)
        mdata_all.adjust_time(day=15)

        if target_grid == 't63grid':
            mdata._apply_mask(get_T63_landseamask(False, area = valid_mask))
            mdata_all._apply_mask(get_T63_landseamask(False, area = valid_mask))
        else:
            tmpmsk = get_generic_landseamask(False,area=valid_mask,target_grid=target_grid)
            mdata._apply_mask(tmpmsk)
            mdata_all._apply_mask(tmpmsk)
            del tmpmsk

        mdata_mean = mdata_all.fldmean()

        #/// return data as a tuple list
        retval = (mdata_all.time,mdata_mean,mdata_all)

        del mdata_all


        return mdata,retval


#-----------------------------------------------------------------------

    def get_snow_fraction(self):
        """
        Specifies for CMIP5 class how to read SNOWFRACTION

        @return: C{Data} object for snow
        """
        data_file = '/net/nas2/export/eo/workspace/m300028/GPA/input/historical_r1i1p1-LR_snow_fract.nc' #todo change this !!!

        #todo: which temporal resolution is needed?? preprocessing with CDO's needed ??? --> monthly

        return Data(data_file,'snow_fract')


#-----------------------------------------------------------------------

    def get_faPAR(self):
        """
        Specifies how to read faPAR information for CMIP5 data
        @return: C{Data} object for faPAR
        """

        ddir = '/net/nas2/export/eo/workspace/m300028/GPA/'   #<<< todo: change this output directory !!!
        data_file = ddir + 'input/historical_r1i1p1-LR_fapar.nc' #todo set inputfilename interactiveley !!!! DUMMY so far for testnig


        #todo: which temporal resolution is needed?? preprocessing with CDO's needed ??? --> monthly

        return Data(data_file,'fapar')

#-----------------------------------------------------------------------




    def get_rainfall_data(self,interval=None):

        '''
        return data object of
        a) seasonal means for precipitation
        b) global mean timeseries for PR at original temporal resolution
        '''

        if interval != 'season':
            raise ValueError, 'Other data than seasonal not supported at the moment for CMIP5 data and rainfall!'


    #original data
        #filename1 = self.data_dir + 'pr/' +  self.model + '/' + 'pr_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
        filename1 = self.data_dir + 'pr/merged/'  + 'pr_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        force_calc = False

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        tmp  = pyCDO(filename1,s_start_time,s_stop_time,force=force_calc).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).seasmean()
        filename = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        if not os.path.exists(filename):
            return None

        pr = Data(filename,'pr',read=True,label=self.model,unit='mm/day',lat_name='lat',lon_name='lon',shift_lon=False,scale_factor=86400.)

        prall  = Data(filename1,'pr',read=True,label=self.model,unit='mm/day',lat_name='lat',lon_name='lon',shift_lon=False,scale_factor=86400.,time_cycle=12)
        prmean = prall.fldmean()

        retval = (prall.time,prmean,prall); del prall

        pr.data = np.ma.array(pr.data,mask=pr.data < 0.)


        return pr,retval

#-----------------------------------------------------------------------

    def get_temperature_2m(self,interval=None):

        """
        return data object of
        a) seasonal means for air temperature
        b) global mean timeseries for TAS at original temporal resolution
        """

        if interval != 'season':
            raise ValueError, 'Other data than seasonal not supported at the moment for CMIP5 data and temperature!'

        #original data
        filename1 = self.data_dir + 'tas/' +  self.model + '/' + 'tas_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        force_calc = False

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        tmp  = pyCDO(filename1,s_start_time,s_stop_time,force=force_calc).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).seasmean()
        filename = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        if not os.path.exists(filename):
            print 'WARNING: Temperature file not found: ', filename
            return None

        tas = Data(filename,'tas',read=True,label=self.model,unit='K',lat_name='lat',lon_name='lon',shift_lon=False)

        tasall = Data(filename1,'tas',read=True,label=self.model,unit='K',lat_name='lat',lon_name='lon',shift_lon=False,time_cycle=12)
        #todo: automatic check that data is monthly! use get_months for that!

        tasmean = tasall.fldmean()

        retval = (tasall.time,tasmean,tasall); del tasall

        tas.data = np.ma.array(tas.data,mask=tas.data < 0.)

        return tas,retval


#-----------------------------------------------------------------------



    def get_surface_shortwave_radiation_down(self,interval = 'season'):

        """
        return data object of
        a) seasonal means for SIS
        b) global mean timeseries for SIS at original temporal resolution
        """

        #original data
        #filename1 = self.data_dir + 'rsds/' +  self.model + '/' + 'rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
        filename1 = self.data_dir + 'rsds/' +  self.experiment + '/ready/' + self.model + '/rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        force_calc = False

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        #/// PREPROCESSING ///
        cdo = Cdo()
        s_start_time = str(self.start_time)[0:10]
        s_stop_time  = str(self.stop_time)[0:10]

        #1) select timeperiod and generatget_she monthly mean file
        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_T63_monmean.nc'
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)

        sys.stdout.write('\n *** Model file monthly: %s\n' % file_monthly)
        cdo.monmean(options='-f nc',output=file_monthly,input = '-remapcon,t63grid -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1, force=force_calc)

        sys.stdout.write('\n *** Reading model data... \n')
        sys.stdout.write('     Interval: ' + interval + '\n')


        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            sis_clim_file     = file_monthly[:-3] + '_ymonmean.nc'
            sis_sum_file      = file_monthly[:-3] + '_ymonsum.nc'
            sis_N_file        = file_monthly[:-3] + '_ymonN.nc'
            sis_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32',output = sis_clim_file,input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32',output = sis_sum_file,input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32',output = sis_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc',output = sis_N_file,input=sis_sum_file + ' ' + sis_clim_file, force=force_calc) #number of samples
        elif interval == 'season':
            sis_clim_file     = file_monthly[:-3] + '_yseasmean.nc'
            sis_sum_file      = file_monthly[:-3] + '_yseassum.nc'
            sis_N_file        = file_monthly[:-3] + '_yseasN.nc'
            sis_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32',output = sis_clim_file,input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32',output = sis_sum_file,input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32',output = sis_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32',output = sis_N_file,input=sis_sum_file + ' ' + sis_clim_file, force=force_calc) #number of samples
        else:
            print interval
            raise ValueError, 'Unknown temporal interval. Can not perform preprocessing! '

        if not os.path.exists(sis_clim_file):
            return None

        #3) read data
        sis = Data(sis_clim_file,'rsds',read=True,label=self.model,unit='$W m^{-2}$',lat_name='lat',lon_name='lon',shift_lon=False)
        sis_std = Data(sis_clim_std_file,'rsds',read=True,label=self.model+ ' std',unit='-',lat_name='lat',lon_name='lon',shift_lon=False)
        sis.std = sis_std.data.copy(); del sis_std
        sis_N = Data(sis_N_file,'rsds',read=True,label=self.model+ ' std',unit='-',lat_name='lat',lon_name='lon',shift_lon=False)
        sis.n = sis_N.data.copy(); del sis_N

        #ensure that climatology always starts with January, therefore set date and then sort
        sis.adjust_time(year=1700,day=15) #set arbitrary time for climatology
        sis.timsort()

        #4) read monthly data
        sisall = Data(file_monthly,'rsds',read=True,label=self.model,unit='W m^{-2}',lat_name='lat',lon_name='lon',shift_lon=False,time_cycle=12) #todo check timecycle
        sisall.adjust_time(day=15)

        sis._apply_mask(get_T63_landseamask(False))
        sisall._apply_mask(get_T63_landseamask(False))

        sismean = sisall.fldmean()

        #/// return data as a tuple list
        retval = (sisall.time,sismean,sisall); del sisall

        #/// mask areas without radiation (set to invalid): all data < 1 W/m**2
        sis.data = np.ma.array(sis.data,mask=sis.data < 1.)

        return sis,retval

#-----------------------------------------------------------------------

    def get_surface_shortwave_radiation_up(self,interval='season'):

        #original data
        filename1 = self.data_dir + 'rsus/' +  self.experiment + '/ready/' + self.model + '/rsus_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        print 'READING SURFACE UPWARD file: ', filename1

        force_calc = False

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        #/// PREPROCESSING ///
        cdo = Cdo()
        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        #1) select timeperiod and generate monthly mean file
        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_T63_monmean.nc'
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)
        cdo.monmean(options='-f nc',output=file_monthly,input = '-remapcon,t63grid -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1,force=force_calc)

        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            sup_clim_file     = file_monthly[:-3] + '_ymonmean.nc'
            sup_sum_file      = file_monthly[:-3] + '_ymonsum.nc'
            sup_N_file        = file_monthly[:-3] + '_ymonN.nc'
            sup_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32',output = sup_clim_file,input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32',output = sup_sum_file,input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32',output = sup_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc',output = sup_N_file,input=sup_sum_file + ' ' + sup_clim_file, force=force_calc) #number of samples
        elif interval == 'season':
            sup_clim_file     = file_monthly[:-3] + '_yseasmean.nc'
            sup_sum_file      = file_monthly[:-3] + '_yseassum.nc'
            sup_N_file        = file_monthly[:-3] + '_yseasN.nc'
            sup_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32',output = sup_clim_file,input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32',output = sup_sum_file,input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32',output = sup_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32',output = sup_N_file,input=sup_sum_file + ' ' + sup_clim_file, force=force_calc) #number of samples
        else:
            print interval
            raise ValueError, 'Unknown temporal interval. Can not perform preprocessing! '



        if not os.path.exists(sup_clim_file):
            print 'File not existing (sup_clim_file): ' + sup_clim_file
            return None


        #3) read data
        sup = Data(sup_clim_file,'rsus',read=True,label=self.model,unit='$W m^{-2}$',lat_name='lat',lon_name='lon',shift_lon=False)
        sup_std = Data(sup_clim_std_file,'rsus',read=True,label=self.model+ ' std',unit='-',lat_name='lat',lon_name='lon',shift_lon=False)
        sup.std = sup_std.data.copy(); del sup_std
        sup_N = Data(sup_N_file,'rsus',read=True,label=self.model+ ' std',unit='-',lat_name='lat',lon_name='lon',shift_lon=False)
        sup.n = sup_N.data.copy(); del sup_N

        #ensure that climatology always starts with January, therefore set date and then sort
        sup.adjust_time(year=1700,day=15) #set arbitrary time for climatology
        sup.timsort()

        #4) read monthly data
        supall = Data(file_monthly,'rsus',read=True,label=self.model,unit='$W m^{-2}$',lat_name='lat',lon_name='lon',shift_lon=False,time_cycle=12) #todo check timecycle
        supmean = supall.fldmean()

        #/// return data as a tuple list
        retval = (supall.time,supmean,supall); del supall

        #/// mask areas without radiation (set to invalid): all data < 1 W/m**2
        #sup.data = np.ma.array(sis.data,mask=sis.data < 1.)

        return sup,retval


#-------------------------------------------------------------------------------------------------------------


    def get_albedo_data(self,interval='season'):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        """

        force_calc = False

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'


        #--- get fluxes
        Fu = self.get_surface_shortwave_radiation_up  (interval=interval)
        if Fu == None:
            print 'File not existing for UPWARD flux!: ', self.name
            return None
        else:
            Fu_i = Fu[0]
        lab = Fu_i.label
        Fd = self.get_surface_shortwave_radiation_down(interval=interval)
        if Fd == None:
            print 'File not existing for DOWNWARD flux!: ', self.name
            return None
        else:
            Fd_i = Fd[0]

        #albedo for chosen interval as caluclated as ratio of means of fluxes in that interval (e.g. season, months)
        Fu_i.div(Fd_i,copy=False); del Fd_i #Fu contains now the albedo
        Fu_i._apply_mask(ls_mask.data)

        #albedo for monthly data (needed for global mean plots )
        Fu_m = Fu[1][2]; del Fu
        Fd_m = Fd[1][2]; del Fd

        #print Fu_m.data.shape
        #print Fd_m.data.shape
        #stop
        Fu_m.div(Fd_m,copy=False); del Fd_m
        Fu_m._apply_mask(ls_mask.data)
        Fu_m._set_valid_range(0.,1.)
        Fu_m.label = lab + ' albedo'
        Fu_i.label = lab + ' albedo'
        Fu_m.unit = '-'
        Fu_i.unit = '-'

        #/// center dates of months
        Fu_m.adjust_time(day=15)
        Fu_i.adjust_time(day=15)

        #/// return data as a tuple list
        retval = (Fu_m.time,Fu_m.fldmean(),Fu_m)



        #downward geht nur bis Jahr 2000: warum ????
        #und upward bis 2008 ???? warum ????


        return Fu_i, retval







#####################################################




class JSBACH_BOT(Model):

    def __init__(self,filename,dic_variables,experiment,name='',shift_lon=False,**kwargs):

        Model.__init__(self,filename,dic_variables,name=name,**kwargs)
        self.experiment = experiment
        self.shift_lon = shift_lon
        self.type = 'JSBACH_BOT'

        self._unique_name = self._get_unique_name()



    def _get_unique_name(self):
        """
        get unique name from model and experiment
        @return: string with unique combination of models and experiment
        """
        return self.name.replace(' ','') + '-' + self.experiment.replace(' ','')


    def get_albedo_data(self,interval='season'):
        """
        get albedo data for JSBACH

        returns Data object
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry'

        v = 'var176'

        filename = self.data_dir + 'data/model1/' + self.experiment + '_echam6_BOT_mm_1979-2006_albedo_yseasmean.nc' #todo: proper files
        ls_mask = get_T63_landseamask(self.shift_lon)

        albedo = Data(filename,v,read=True,
        label='MPI-ESM albedo ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data)

        return albedo



    def get_tree_fraction(self,interval='season'):
        """
        todo implement this for data from a real run !!!
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry'

        ls_mask = get_T63_landseamask(self.shift_lon)

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_forest_shrub.nc'
        v = 'var12'
        tree = Data(filename,v,read=True,
        label='MPI-ESM tree fraction ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data,start_time = pl.num2date(pl.datestr2num('2001-01-01')),stop_time=pl.num2date(pl.datestr2num('2001-12-31')))

        return tree

    def get_grass_fraction(self,interval='season'):
        """
        todo implement this for data from a real run !!!
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry'

        ls_mask = get_T63_landseamask(self.shift_lon)

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_grass_crop_pasture_2001.nc'
        v = 'var12'
        grass = Data(filename,v,read=True,
        label='MPI-ESM tree fraction ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        #shift_lon=shift_lon,
        mask=ls_mask.data.data,start_time = pl.num2date(pl.datestr2num('2001-01-01')),stop_time=pl.num2date(pl.datestr2num('2001-12-31')) , squeeze=True  )


        return grass









    def get_surface_shortwave_radiation_down(self,interval = 'season'):
        """
        get surface shortwave incoming radiation data for JSBACH

        returns Data object
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry'

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
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SIS data
        sis = Data(filename,v,read=True,
        label='MPI-ESM SIS ' + self.experiment, unit = '-',lat_name='lat',lon_name='lon',
        #shift_lon=shift_lon,
        mask=ls_mask.data.data)

        return sis


    def get_rainfall_data(self,interval='season'):
        """
        get rainfall data for JSBACH
        returns Data object
        """

        if interval == 'season':
            pass
        else:
            raise ValueError, 'Invalid value for interval: ' + interval

        #/// PREPROCESSING: seasonal means ///
        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        filename1 = self.data_dir +  self.experiment + '_echam6_BOT_mm_1980_sel.nc'
        tmp  = pyCDO(filename1,s_start_time,s_stop_time).seldate()
        tmp1 = pyCDO(tmp,s_start_time,s_stop_time).seasmean()
        filename = pyCDO(tmp1,s_start_time,s_stop_time).yseasmean()

        #/// READ DATA ///

        #1) land / sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #2) precipitation data
        try:
            v = 'var4'
            rain = Data(filename,v,read=True,scale_factor = 86400.,
            label='MPI-ESM ' + self.experiment, unit = 'mm/day',lat_name='lat',lon_name='lon',
            shift_lon=self.shift_lon,
            mask=ls_mask.data.data)
        except:
            v='var142'
            rain = Data(filename,v,read=True,scale_factor = 86400.,
            label='MPI-ESM ' + self.experiment, unit = 'mm/day',lat_name='lat',lon_name='lon',
            shift_lon=self.shift_lon,
            mask=ls_mask.data.data)

        return rain



#-----------------------------------------------------------------------

class JSBACH_RAW2(Model):
    """
    Class for RAW JSBACH model output
    """

    def __init__(self,filename,dic_variables,experiment,name='',shift_lon=False,model_dict=None,**kwargs):

        Model.__init__(self,filename,dic_variables,name=name,**kwargs)

        self.experiment = experiment
        self.shift_lon = shift_lon
        #self.get_data()
        self.type = 'JSBACH_RAW2'

        self._unique_name = self._get_unique_name()

        #--- do preprocessing of streams (only needed once!) ---
        self.files={}
        self._preproc_streams()
        self.model_dict=copy.deepcopy(model_dict)

        self.model = 'JSBACH'



    def _preproc_streams(self):
        """
        It is assumed that the standard JSBACH postprocessing scripts have been applied.
        Thus monthly mean data is available for each stream and code tables still need to be applied.

        This routine does the following:
        1) merge all times from individual (monthly mean) output files
        2) assign codetables to work with proper variable names
        3) aggregate data from tiles to gridbox values
        """

        print 'Preprocessing JSBACH raw data streams (may take a while) ...'

        cdo = Cdo()

        # --- jsbach stream
        print '   JSBACH stream ...'
        outfile = get_temporary_directory()  + self.experiment + '_jsbach_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach.codes'
            tmp = tempfile.mktemp(suffix='.nc',prefix=self.experiment + '_jsbach_',dir=get_temporary_directory()) #temporary file
            cdo.mergetime(options='-f nc',output=tmp,input=self.data_dir + 'outdata/jsbach/' + self.experiment + '_jsbach_main_mm_*.grb')
            cdo.monmean(options='-f nc' ,output=outfile,input='-setpartab,' + codetable + ' ' + tmp) #monmean needed here, as otherwise interface does not work
            os.remove(tmp)
        self.files.update({'jsbach':outfile})

        #--- veg stream
        print '   VEG stream ...'
        outfile = get_temporary_directory()  + self.experiment + '_jsbach_veg_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach_veg.codes'
            tmp = tempfile.mktemp(suffix='.nc',prefix=self.experiment + '_jsbach_veg_',dir=get_temporary_directory()) #temporary file
            cdo.mergetime(options='-f nc',output=tmp,input=self.data_dir + 'outdata/jsbach/' + self.experiment + '_jsbach_veg_mm_*.grb')
            cdo.monmean(options='-f nc' ,output=outfile,input='-setpartab,' + codetable + ' ' + tmp) #monmean needed here, as otherwise interface does not work
            os.remove(tmp)
        self.files.update({'veg':outfile})

        #--- veg land
        print '   LAND stream ...'
        outfile = get_temporary_directory()  + self.experiment + '_jsbach_land_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach_land.codes'
            tmp = tempfile.mktemp(suffix='.nc',prefix=self.experiment + '_jsbach_land_',dir=get_temporary_directory()) #temporary file
            cdo.mergetime(options='-f nc',output=tmp,input=self.data_dir + 'outdata/jsbach/' + self.experiment + '_jsbach_land_mm_*.grb')
            cdo.monmean(options='-f nc' ,output=outfile,input='-setpartab,' + codetable + ' ' + tmp) #monmean needed here, as otherwise interface does not work
            os.remove(tmp)
        self.files.update({'land':outfile})

        #--- surf stream
        print '   SURF stream ...'
        outfile = get_temporary_directory()  + self.experiment + '_jsbach_surf_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach_surf.codes'
            tmp = tempfile.mktemp(suffix='.nc',prefix=self.experiment + '_jsbach_surf_',dir=get_temporary_directory()) #temporary file
            cdo.mergetime(options='-f nc',output=tmp,input=self.data_dir + 'outdata/jsbach/' + self.experiment + '_jsbach_surf_mm_*.grb')
            cdo.monmean(options='-f nc' ,output=outfile,input='-setpartab,' + codetable + ' ' + tmp) #monmean needed here, as otherwise interface does not work
            os.remove(tmp)
        self.files.update({'surf':outfile})

        #--- ECHAM BOT stream
        print '   BOT stream ...'
        outfile = get_temporary_directory()  + self.experiment + '_echam6_echam_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_echam6_echam.codes'
            tmp = tempfile.mktemp(suffix='.nc',prefix=self.experiment + '_echam6_echam_',dir=get_temporary_directory()) #temporary file
            cdo.mergetime(options='-f nc',output=tmp,input=self.data_dir + 'outdata/echam6/' + self.experiment + '_echam6_BOT_mm_*.sz')
            cdo.monmean(options='-f nc' ,output=outfile,input='-setpartab,' + codetable + ' ' + tmp) #monmean needed here, as otherwise interface does not work
            os.remove(tmp)
        self.files.update({'echam':outfile})




    def _get_unique_name(self):
        """
        get unique name from model and experiment
        @return: string with unique combination of models and experiment
        """
        return self.name.replace(' ','') + '-' + self.experiment.replace(' ','')


    def get_albedo_data(self, interval='season'):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        """

        if self.start_time is None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time is None:
            raise ValueError, 'Stop time needs to be specified'

        print ''
        print 'in get_albedo() before call: ', self.model_dict['sis']
        print ''

        sw_down = self.get_surface_shortwave_radiation_down(interval=interval)
        sw_up   = self.get_surface_shortwave_radiation_up  (interval=interval)

        #climatological mean
        alb     = sw_up[0].div(sw_down[0])
        alb.label = self.experiment + ' albedo'
        alb.unit = '-'

        #original data
        alb_org  = sw_up[1][2].div(sw_down[1][2])
        alb_org.label = self.experiment + ' albedo'
        alb_org.unit = '-'

        retval = (alb_org.time,alb_org.fldmean(),alb_org)

        return alb, retval

    def get_surface_shortwave_radiation_up(self,interval='season'):
        tmpdict = copy.deepcopy(self.model_dict['surface_upward_flux'])
        return self.get_jsbach_data_generic(interval=interval,**tmpdict)

    def get_surface_shortwave_radiation_down(self,interval='season'):
        tmpdict = copy.deepcopy(self.model_dict['sis'])
        return self.get_jsbach_data_generic(interval=interval,**tmpdict)

    def get_rainfall_data(self,interval='season'):
        tmpdict = copy.deepcopy(self.model_dict['rain'])
        return self.get_jsbach_data_generic(interval=interval,**tmpdict)

    def get_temperature_2m(self,interval='season'):
        tmpdict = copy.deepcopy(self.model_dict['temperature'])
        return self.get_jsbach_data_generic(interval=interval,**tmpdict)

    def get_jsbach_data_generic(self,interval='season', **kwargs):
        """
        unique parameters are:
            filename - file basename
            variable - name of the variable as the short_name in the netcdf file

            kwargs is a dictionary with keys for each model. Then a dictionary with properties follows

        """

        if not self.type in kwargs.keys():
            print 'WARNING: it is not possible to get data using generic function, as method missing: ', self.type, kwargs.keys()
            return None

        print self.type
        print kwargs

        locdict = kwargs[self.type]

        # read settings and details from the keyword arguments
        # no defaults; everything should be explicitely specified in either the config file or the dictionaries

        varname  = locdict.pop('variable')
        units    = locdict.pop('unit', 'Crazy Unit')
        #interval = kwargs.pop('interval') #, 'season') #does not make sense to specifiy a default value as this option is specified by configuration file!

        lat_name     = locdict.pop('lat_name', 'lat')
        lon_name     = locdict.pop('lon_name', 'lon')
        #model_suffix = locdict.pop('model_suffix')
        #model_prefix = locdict.pop('model_prefix')
        file_format  = locdict.pop('file_format')
        scf = locdict.pop('scale_factor')
        valid_mask    = locdict.pop('valid_mask')
        custom_path  = locdict.pop('custom_path', None)
        thelevel  = locdict.pop('level', None)

        target_grid   = self._actplot_options['targetgrid']
        interpolation = self._actplot_options['interpolation']

        if self.type != 'JSBACH_RAW2':
            print self.type
            raise ValueError, 'Invalid data format here!'


        if varname in ['swdown_acc','swdown_reflect_acc']:
            filename1= self.files['jsbach']
        elif varname in ['precip_acc']:
            filename1 = self.files['land']
        elif varname in ['temp2']:
            filename1 = self.files['echam']
        else:
            print varname
            raise ValueError, 'Unknown variable type for JSBACH_RAW2 processing!'


        force_calc = False

        if self.start_time == None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time == None:
            raise ValueError, 'Stop time needs to be specified'

        #/// PREPROCESSING ///
        cdo = Cdo()
        s_start_time = str(self.start_time)[0:10]
        s_stop_time  = str(self.stop_time )[0:10]

        #1) select timeperiod and generate monthly mean file
        if target_grid == 't63grid':
            gridtok = 'T63'
        else:
            gridtok = 'SPECIAL_GRID'

        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_' + gridtok + '_monmean.nc' #target filename
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)

        sys.stdout.write('\n *** Model file monthly: %s\n' % file_monthly)

        if not os.path.exists(filename1):
            print 'WARNING: File not existing: ' + filename1
            return None

        cdo.monmean(options='-f nc',output=file_monthly,input = '-' + interpolation + ',' + target_grid + ' -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1, force=force_calc)

        sys.stdout.write('\n *** Reading model data... \n')
        sys.stdout.write('     Interval: ' + interval + '\n')

        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            mdata_clim_file     = file_monthly[:-3] + '_ymonmean.nc'
            mdata_sum_file      = file_monthly[:-3] + '_ymonsum.nc'
            mdata_N_file        = file_monthly[:-3] + '_ymonN.nc'
            mdata_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32',output = mdata_clim_file,input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32',output = mdata_sum_file,input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32',output = mdata_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc',output = mdata_N_file,input=mdata_sum_file + ' ' + mdata_clim_file, force=force_calc) #number of samples
        elif interval == 'season':
            mdata_clim_file     = file_monthly[:-3] + '_yseasmean.nc'
            mdata_sum_file      = file_monthly[:-3] + '_yseassum.nc'
            mdata_N_file        = file_monthly[:-3] + '_yseasN.nc'
            mdata_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32',output = mdata_clim_file,input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32',output = mdata_sum_file,input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32',output = mdata_clim_std_file,input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32',output = mdata_N_file,input=mdata_sum_file + ' ' + mdata_clim_file, force=force_calc) #number of samples
        else:
            raise ValueError, 'Unknown temporal interval. Can not perform preprocessing! '

        if not os.path.exists(mdata_clim_file):
            return None

        #3) read data
        if interval == 'monthly':
            thetime_cylce = 12
        elif interval == 'season':
            thetime_cylce = 4
        else:
            print interval
            raise ValueError, 'Unsupported interval!'
        mdata = Data(mdata_clim_file,varname,read=True,label=self.model,unit=units,lat_name=lat_name,lon_name=lon_name,shift_lon=False, scale_factor=scf,level=thelevel,time_cycle=thetime_cylce)
        mdata_std = Data(mdata_clim_std_file,varname,read=True,label=self.model+ ' std',unit='-',lat_name=lat_name,lon_name=lon_name,shift_lon=False,level=thelevel,time_cycle=thetime_cylce)
        mdata.std = mdata_std.data.copy(); del mdata_std
        mdata_N = Data(mdata_N_file,varname,read=True,label=self.model+ ' std',unit='-',lat_name=lat_name,lon_name=lon_name,shift_lon=False, scale_factor=scf,level=thelevel)
        mdata.n = mdata_N.data.copy(); del mdata_N

        #ensure that climatology always starts with J  anuary, therefore set date and then sort
        mdata.adjust_time(year=1700,day=15) #set arbitrary time for climatology
        mdata.timsort()

        #4) read monthly data
        mdata_all = Data(file_monthly,varname,read=True,label=self.model,unit=units,lat_name=lat_name,lon_name=lon_name,shift_lon=False,time_cycle=12,scale_factor=scf,level=thelevel)
        mdata_all.adjust_time(day=15)

        if target_grid == 't63grid':
            mdata._apply_mask(get_T63_landseamask(False, area = valid_mask))
            mdata_all._apply_mask(get_T63_landseamask(False, area = valid_mask))
        else:
            tmpmsk = get_generic_landseamask(False,area=valid_mask,target_grid=target_grid)
            mdata._apply_mask(tmpmsk)
            mdata_all._apply_mask(tmpmsk)
            del tmpmsk

        mdata_mean = mdata_all.fldmean()

        #/// return data as a tuple list
        retval = (mdata_all.time,mdata_mean,mdata_all)

        del mdata_all


        return mdata,retval





#-----------------------------------------------------------------------


class JSBACH_RAW(Model):
    """
    Class for RAW JSBACH model output
    """

    def __init__(self,filename,dic_variables,experiment,name='',shift_lon=False,**kwargs):

        Model.__init__(self,filename,dic_variables,name=name,**kwargs)

        self.experiment = experiment
        self.shift_lon = shift_lon
        #self.get_data()
        self.type = 'JSBACH_RAW'

        self._unique_name = self._get_unique_name()

    def _get_unique_name(self):
        """
        get unique name from model and experiment
        @return: string with unique combination of models and experiment
        """
        return self.name.replace(' ','') + '-' + self.experiment.replace(' ','')


    def get_temperature_2m(self,interval = 'season'):
        """
        get surface temperature (2m) from JSBACH model results

        @param interval: specifies the aggregation interval. Possible options: ['season']
        @type interval: str

        @return: returns a C{Data} object
        @rtype: C{Data}
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH RAW files, sorry'

        v = 'surface_temperature'

        y1 = '1979-01-01'; y2 = '2010-12-31'
        rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_land_' + y1[0:4] + '_' + y2[0:4] + '.nc'


        if not os.path.exists(rawfilename):
            print 'File not existing (rawfile): ', rawfilename
            return None

        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SIS data
        t2m = Data(filename,v,read=True,
            label=self.experiment + ' ' + v, unit = 'K',lat_name='lat',lon_name='lon',
            shift_lon=self.shift_lon,
            mask=ls_mask.data.data)


        retval = (t2m.time,t2m,t2m)

        return t2m, retval



#-----------------------------------------------------------------------

    def get_albedo_data(self, interval='season'):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH RAW files, sorry'

        if self.start_time is None:
            raise ValueError, 'Start time needs to be specified'
        if self.stop_time is None:
            raise ValueError, 'Stop time needs to be specified'

        sw_down = self.get_surface_shortwave_radiation_down()
        sw_up   = self.get_surface_shortwave_radiation_up()
        alb     = sw_up.div(sw_down)
        alb.label = self.experiment + ' albedo'
        alb.unit = '-'

        return alb

#-----------------------------------------------------------------------


    def get_surface_shortwave_radiation_down(self,interval = 'season'):
        """
        get surface shortwave incoming radiation data for JSBACH

        @param interval: specifies the aggregation interval. Possible options: ['season']
        @type interval: str

        @return: returns a C{Data} object
        @rtype: C{Data}
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH RAW files, sorry'

        v = 'swdown_acc'

        #y1 = '1979-01-01'; y2 = '2010-12-31'
        y1 = '1980-01-01'; y2 = '2010-12-31'
        rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '.nc'
        #rawfilename = self.data_dir +  self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '_yseasmean.nc'

        if not os.path.exists(rawfilename):
            print 'File not existing: ', rawfilename
            return None

        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SIS data
        sw_down = Data(filename,v,read=True,
        label=self.experiment + ' ' + v, unit = 'W/m**2',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data)

        return sw_down


#-----------------------------------------------------------------------

    def get_surface_shortwave_radiation_up(self,interval = 'season'):
        """
        get surface shortwave upward radiation data for JSBACH

        returns Data object

        todo CDO preprocessing of seasonal means
        todo temporal aggregation of data --> or leave it to the user!
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH RAW files, sorry'

        v = 'swdown_reflect_acc'

        #y1 = '1979-01-01'; y2 = '2010-12-31'
        y1 = '1980-01-01'; y2 = '2010-12-31'
        rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '.nc'
        #rawfilename = self.data_dir +  self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '_yseasmean.nc'

        if not os.path.exists(rawfilename):
            print 'File not existing: ', rawfilename
            return None

        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SW up data
        sw_up = Data(filename,v,read=True,
        label=self.experiment + ' ' + v, unit = 'W/m**2',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data)

        return sw_up

#-----------------------------------------------------------------------

    def get_rainfall_data(self,interval = 'season'):
        """
        get surface rainfall data for JSBACH

        returns Data object

        todo CDO preprocessing of seasonal means
        todo temporal aggregation of data --> or leave it to the user!
        """

        if interval != 'season':
            raise ValueError, 'Other temporal sampling than SEASON not supported yet for JSBACH RAW files, sorry'

        v = 'precip_acc'

        y1 = '1979-01-01'; y2 = '2010-12-31'
        #rawfilename = self.data_dir + 'yseasmean_' + self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '.nc'
        rawfilename = self.data_dir +  self.experiment + '_jsbach_' + y1[0:4] + '_' + y2[0:4] + '_yseasmean.nc'

        if not os.path.exists(rawfilename):
            return None

        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SW up data
        rain = Data(filename,v,read=True,
        label=self.experiment + ' ' + v, unit = 'mm/day',lat_name='lat',lon_name='lon',
        shift_lon=self.shift_lon,
        mask=ls_mask.data.data,scale_factor = 86400.)

        return rain


#-----------------------------------------------------------------------



class CMIP3Data(CMIP5Data):
    """
    Class for CMIP3 model simulations. This class is derived from C{Model}.
    """
    def __init__(self,data_dir,model,experiment,dic_variables,name='',shift_lon=False,**kwargs):
        """

        @param data_dir: directory that specifies the root directory where the data is located
        @param model: TBD tood
        @param experiment: specifies the ID of the experiment (str)
        @param dic_variables:
        @param name: TBD todo
        @param shift_lon: specifies if longitudes of data need to be shifted
        @param kwargs: other keyword arguments
        @return:
        """
        Model.__init__(self,None,dic_variables,name=model,shift_lon=shift_lon,**kwargs)

        self.model      = model; self.experiment = experiment
        self.data_dir   = data_dir; self.shift_lon  = shift_lon
        self.type       = 'CMIP3'

        self._unique_name = self._get_unique_name()


