# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from cdo import Cdo
from pycmbs.data import Data
import tempfile as tempfile
import copy
import glob
import os
import sys
import numpy as np

from pycmbs.benchmarking import preprocessor
from pycmbs.benchmarking.utils import get_T63_landseamask, get_temporary_directory
from pycmbs.benchmarking.models.model_basic import *


class CMIP5Data(Model):
    """
    Class for CMIP5 model simulations. This class is derived from C{Model}.
    """
    def __init__(self, data_dir, model, experiment, dic_variables, name='', shift_lon=False, **kwargs):
        """
        Parameters
        ----------
        data_dir : str
            directory that specifies the root directory where the data is located
        model : TBD todo
        experiment : str
            specifies the ID of the experiment
        dic_variables : TODO
        name : str
            name of model
        shift_lon : bool
            specifies if longitudes of data need to be shifted
        kwargs : dict
            other keyword arguments
        """
        super(CMIP5Data, self).__init__(data_dir, dic_variables, name=model, shift_lon=shift_lon, **kwargs)

        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon
        self.type = 'CMIP5'
        self._unique_name = self._get_unique_name()

    def _get_unique_name(self):
        """
        get unique name from model and experiment

        Returns
        -------
        string with unique combination of models and experiment
        """
        return self.model.replace(' ', '') + '-' + self.experiment.replace(' ', '')

#-----------------------------------------------------------------------

    def get_model_data_generic(self, interval='season', **kwargs):
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
        varname = locdict.pop('variable')
        units = locdict.pop('unit', 'Crazy Unit')
        #interval = kwargs.pop('interval') #, 'season') #does not make sense to specifiy a default value as this option is specified by configuration file!

        lat_name = locdict.pop('lat_name', 'lat')
        lon_name = locdict.pop('lon_name', 'lon')
        model_suffix = locdict.pop('model_suffix')
        model_prefix = locdict.pop('model_prefix')
        file_format = locdict.pop('file_format')
        scf = locdict.pop('scale_factor')
        valid_mask = locdict.pop('valid_mask')
        custom_path = locdict.pop('custom_path', None)
        thelevel = locdict.pop('level', None)

        target_grid = self._actplot_options['targetgrid']
        interpolation = self._actplot_options['interpolation']

        if custom_path is None:
            filename1 = ("%s%s/merged/%s_%s_%s_%s_%s.%s" %
                        (self.data_dir, varname, varname, model_prefix, self.model, self.experiment, model_suffix, file_format))
        else:
            if self.type == 'CMIP5':
                filename1 = ("%s/%s_%s_%s_%s_%s.%s" %
                             (custom_path, varname, model_prefix, self.model, self.experiment, model_suffix, file_format))
            elif self.type == 'CMIP5RAW':
                filename1 = ("%s/%s_%s_%s_%s_%s.%s" %
                             (custom_path, varname, model_prefix, self.model, self.experiment, model_suffix, file_format))
            elif self.type == 'CMIP3':
                filename1 = ("%s/%s_%s_%s_%s.%s" %
                             (custom_path, self.experiment, self.model, varname, model_suffix, file_format))
            else:
                print self.type
                raise ValueError('Can not generate filename: invalid model type! %s' % self.type)

        force_calc = False

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        #/// PREPROCESSING ///
        cdo = Cdo()
        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        #1) select timeperiod and generate monthly mean file
        if target_grid == 't63grid':
            gridtok = 'T63'
        else:
            gridtok = 'SPECIAL_GRID'

        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_' + gridtok + '_monmean.nc'  # target filename
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)

        sys.stdout.write('\n *** Model file monthly: %s\n' % file_monthly)

        if not os.path.exists(filename1):
            print 'WARNING: File not existing: ' + filename1
            return None

        cdo.monmean(options='-f nc', output=file_monthly, input='-' + interpolation + ',' + target_grid + ' -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1, force=force_calc)

        sys.stdout.write('\n *** Reading model data... \n')
        sys.stdout.write('     Interval: ' + interval + '\n')

        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            mdata_clim_file = file_monthly[:-3] + '_ymonmean.nc'
            mdata_sum_file = file_monthly[:-3] + '_ymonsum.nc'
            mdata_N_file = file_monthly[:-3] + '_ymonN.nc'
            mdata_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32', output=mdata_clim_file, input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32', output=mdata_sum_file, input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32', output=mdata_clim_std_file, input=file_monthly, force=force_calc)
            cdo.div(options='-f nc', output=mdata_N_file, input=mdata_sum_file + ' ' + mdata_clim_file, force=force_calc)  # number of samples
        elif interval == 'season':
            mdata_clim_file = file_monthly[:-3] + '_yseasmean.nc'
            mdata_sum_file = file_monthly[:-3] + '_yseassum.nc'
            mdata_N_file = file_monthly[:-3] + '_yseasN.nc'
            mdata_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32', output=mdata_clim_file, input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32', output=mdata_sum_file, input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32', output=mdata_clim_std_file, input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32', output=mdata_N_file, input=mdata_sum_file + ' ' + mdata_clim_file, force=force_calc)  # number of samples
        else:
            raise ValueError('Unknown temporal interval. Can not perform preprocessing!')

        if not os.path.exists(mdata_clim_file):
            return None

        #3) read data
        if interval == 'monthly':
            thetime_cylce = 12
        elif interval == 'season':
            thetime_cylce = 4
        else:
            print interval
            raise ValueError('Unsupported interval!')
        mdata = Data(mdata_clim_file, varname, read=True, label=self.model, unit=units, lat_name=lat_name, lon_name=lon_name, shift_lon=False, scale_factor=scf, level=thelevel, time_cycle=thetime_cylce)
        mdata_std = Data(mdata_clim_std_file, varname, read=True, label=self.model + ' std', unit='-', lat_name=lat_name, lon_name=lon_name, shift_lon=False, level=thelevel, time_cycle=thetime_cylce)
        mdata.std = mdata_std.data.copy()
        del mdata_std
        mdata_N = Data(mdata_N_file, varname, read=True, label=self.model + ' std', unit='-', lat_name=lat_name, lon_name=lon_name, shift_lon=False, scale_factor=scf, level=thelevel)
        mdata.n = mdata_N.data.copy()
        del mdata_N

        #ensure that climatology always starts with January, therefore set date and then sort
        mdata.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
        mdata.timsort()

        #4) read monthly data
        mdata_all = Data(file_monthly, varname, read=True, label=self.model, unit=units, lat_name=lat_name, lon_name=lon_name, shift_lon=False, time_cycle=12, scale_factor=scf, level=thelevel)
        mdata_all.adjust_time(day=15)

        #mask_antarctica masks everything below 60 degrees S.
        #here we only mask Antarctica, if only LAND points shall be used
        if valid_mask == 'land':
            mask_antarctica = True
        elif valid_mask == 'ocean':
            mask_antarctica = False
        else:
            mask_antarctica = False

        if target_grid == 't63grid':
            mdata._apply_mask(get_T63_landseamask(False, area=valid_mask, mask_antarctica=mask_antarctica))
            mdata_all._apply_mask(get_T63_landseamask(False, area=valid_mask, mask_antarctica=mask_antarctica))
        else:
            tmpmsk = get_generic_landseamask(False, area=valid_mask, target_grid=target_grid, mask_antarctica=mask_antarctica)
            mdata._apply_mask(tmpmsk)
            mdata_all._apply_mask(tmpmsk)
            del tmpmsk

        mdata_mean = mdata_all.fldmean()

        # return data as a tuple list
        retval = (mdata_all.time, mdata_mean, mdata_all)

        del mdata_all
        return mdata, retval

#-----------------------------------------------------------------------

    def get_snow_fraction(self):
        """
        Specifies for CMIP5 class how to read SNOWFRACTION

        @return: C{Data} object for snow
        """
        data_file = '/net/nas2/export/eo/workspace/m300028/GPA/input/historical_r1i1p1-LR_snow_fract.nc'  # todo change this !!!

        #todo: which temporal resolution is needed?? preprocessing with CDO's needed ??? --> monthly

        return Data(data_file, 'snow_fract')

#-----------------------------------------------------------------------

    def get_faPAR(self):
        """
        Specifies how to read faPAR information for CMIP5 data
        @return: C{Data} object for faPAR
        """

        ddir = '/net/nas2/export/eo/workspace/m300028/GPA/'  # TODO <<< todo: change this output directory !!!
        data_file = ddir + 'input/historical_r1i1p1-LR_fapar.nc'  # TODO todo set inputfilename interactiveley !!!! DUMMY so far for testnig

        #todo: which temporal resolution is needed?? preprocessing with CDO's needed ??? --> monthly
        return Data(data_file, 'fapar')

#-----------------------------------------------------------------------

    def get_temperature_2m(self, interval=None):
        """
        return data object of
        a) seasonal means for air temperature
        b) global mean timeseries for TAS at original temporal resolution
        """

        if interval != 'season':
            raise ValueError('Other data than seasonal not supported at the moment for CMIP5 data and temperature!')

        #original data
        filename1 = self.data_dir + 'tas/' + self.model + '/' + 'tas_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'

        force_calc = False

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        tmp = pyCDO(filename1, s_start_time, s_stop_time, force=force_calc).seldate()
        tmp1 = pyCDO(tmp, s_start_time, s_stop_time).seasmean()
        filename = pyCDO(tmp1, s_start_time, s_stop_time).yseasmean()

        if not os.path.exists(filename):
            print 'WARNING: Temperature file not found: ', filename
            return None

        tas = Data(filename, 'tas', read=True, label=self.model, unit='K', lat_name='lat', lon_name='lon', shift_lon=False)

        tasall = Data(filename1, 'tas', read=True, label=self.model, unit='K', lat_name='lat', lon_name='lon', shift_lon=False)
        if tasall.time_cycle != 12:
            raise ValueError('Timecycle of 12 expected here!')

        tasmean = tasall.fldmean()
        retval = (tasall.time, tasmean, tasall)
        del tasall

        tas.data = np.ma.array(tas.data, mask=tas.data < 0.)

        return tas, retval

#-----------------------------------------------------------------------

    def get_surface_shortwave_radiation_down(self, interval='season', force_calc=False, **kwargs):
        """
        return data object of
        a) seasonal means for SIS
        b) global mean timeseries for SIS at original temporal resolution
        """

        the_variable = 'rsds'

        locdict = kwargs[self.type]
        valid_mask = locdict.pop('valid_mask')

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        if self.type == 'CMIP5':
            filename1 = self.data_dir + 'rsds/' + self.experiment + '/ready/' + self.model + '/rsds_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
        elif self.type == 'CMIP5RAW':  # raw CMIP5 data based on ensembles
            filename1 = self._get_ensemble_filename('rsds')
        else:
            raise ValueError('Unknown model type! not supported here!')

        if not os.path.exists(filename1):
            print ('WARNING file not existing: %s' % filename1)
            return None

        #/// PREPROCESSING ///
        cdo = Cdo()

        #1) select timeperiod and generatget_she monthly mean file
        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_T63_monmean.nc'
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)

        print file_monthly

        sys.stdout.write('\n *** Model file monthly: %s\n' % file_monthly)
        cdo.monmean(options='-f nc', output=file_monthly, input='-remapcon,t63grid -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1, force=force_calc)

        sys.stdout.write('\n *** Reading model data... \n')
        sys.stdout.write('     Interval: ' + interval + '\n')

        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            sis_clim_file = file_monthly[:-3] + '_ymonmean.nc'
            sis_sum_file = file_monthly[:-3] + '_ymonsum.nc'
            sis_N_file = file_monthly[:-3] + '_ymonN.nc'
            sis_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32', output=sis_clim_file, input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32', output=sis_sum_file, input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32', output=sis_clim_std_file, input=file_monthly, force=force_calc)
            cdo.div(options='-f nc', output=sis_N_file, input=sis_sum_file + ' ' + sis_clim_file, force=force_calc)  # number of samples
        elif interval == 'season':
            sis_clim_file = file_monthly[:-3] + '_yseasmean.nc'
            sis_sum_file = file_monthly[:-3] + '_yseassum.nc'
            sis_N_file = file_monthly[:-3] + '_yseasN.nc'
            sis_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32', output=sis_clim_file, input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32', output=sis_sum_file, input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32', output=sis_clim_std_file, input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32', output=sis_N_file, input=sis_sum_file + ' ' + sis_clim_file, force=force_calc)  # number of samples
        else:
            print interval
            raise ValueError('Unknown temporal interval. Can not perform preprocessing!')

        if not os.path.exists(sis_clim_file):
            return None

        #3) read data
        sis = Data(sis_clim_file, 'rsds', read=True, label=self.model, unit='$W m^{-2}$', lat_name='lat', lon_name='lon', shift_lon=False)
        sis_std = Data(sis_clim_std_file, 'rsds', read=True, label=self.model + ' std', unit='-', lat_name='lat', lon_name='lon', shift_lon=False)
        sis.std = sis_std.data.copy()
        del sis_std
        sis_N = Data(sis_N_file, 'rsds', read=True, label=self.model + ' std', unit='-', lat_name='lat', lon_name='lon', shift_lon=False)
        sis.n = sis_N.data.copy()
        del sis_N

        #ensure that climatology always starts with January, therefore set date and then sort
        sis.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
        sis.timsort()

        #4) read monthly data
        sisall = Data(file_monthly, 'rsds', read=True, label=self.model, unit='W m^{-2}', lat_name='lat', lon_name='lon', shift_lon=False)
        if not sisall._is_monthly():
            raise ValueError('Timecycle of 12 expected here!')
        sisall.adjust_time(day=15)

        # land/sea masking ...
        if valid_mask == 'land':
            mask_antarctica = True
        elif valid_mask == 'ocean':
            mask_antarctica = False
        else:
            mask_antarctica = False

        sis._apply_mask(get_T63_landseamask(False, mask_antarctica=mask_antarctica, area=valid_mask))
        sisall._apply_mask(get_T63_landseamask(False, mask_antarctica=mask_antarctica, area=valid_mask))
        sismean = sisall.fldmean()

        # return data as a tuple list
        retval = (sisall.time, sismean, sisall)
        del sisall

        # mask areas without radiation (set to invalid): all data < 1 W/m**2
        sis.data = np.ma.array(sis.data, mask=sis.data < 1.)

        return sis, retval

#-----------------------------------------------------------------------

    def get_surface_shortwave_radiation_up(self, interval='season', force_calc=False, **kwargs):

        if self.type == 'CMIP5':
            filename1 = self.data_dir + 'rsus/' + self.experiment + '/ready/' + self.model + '/rsus_Amon_' + self.model + '_' + self.experiment + '_ensmean.nc'
        elif self.type == 'CMIP5RAW':  # raw CMIP5 data based on ensembles
            filename1 = self._get_ensemble_filename('rsus')
        else:
            raise ValueError('Unknown type! not supported here!')

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        if not os.path.exists(filename1):
            print ('WARNING file not existing: %s' % filename1)
            return None

        # PREPROCESSING
        cdo = Cdo()
        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        #1) select timeperiod and generate monthly mean file
        file_monthly = filename1[:-3] + '_' + s_start_time + '_' + s_stop_time + '_T63_monmean.nc'
        file_monthly = get_temporary_directory() + os.path.basename(file_monthly)
        cdo.monmean(options='-f nc', output=file_monthly, input='-remapcon,t63grid -seldate,' + s_start_time + ',' + s_stop_time + ' ' + filename1, force=force_calc)

        #2) calculate monthly or seasonal climatology
        if interval == 'monthly':
            sup_clim_file = file_monthly[:-3] + '_ymonmean.nc'
            sup_sum_file = file_monthly[:-3] + '_ymonsum.nc'
            sup_N_file = file_monthly[:-3] + '_ymonN.nc'
            sup_clim_std_file = file_monthly[:-3] + '_ymonstd.nc'
            cdo.ymonmean(options='-f nc -b 32', output=sup_clim_file, input=file_monthly, force=force_calc)
            cdo.ymonsum(options='-f nc -b 32', output=sup_sum_file, input=file_monthly, force=force_calc)
            cdo.ymonstd(options='-f nc -b 32', output=sup_clim_std_file, input=file_monthly, force=force_calc)
            cdo.div(options='-f nc', output=sup_N_file, input=sup_sum_file + ' ' + sup_clim_file, force=force_calc)  # number of samples
        elif interval == 'season':
            sup_clim_file = file_monthly[:-3] + '_yseasmean.nc'
            sup_sum_file = file_monthly[:-3] + '_yseassum.nc'
            sup_N_file = file_monthly[:-3] + '_yseasN.nc'
            sup_clim_std_file = file_monthly[:-3] + '_yseasstd.nc'
            cdo.yseasmean(options='-f nc -b 32', output=sup_clim_file, input=file_monthly, force=force_calc)
            cdo.yseassum(options='-f nc -b 32', output=sup_sum_file, input=file_monthly, force=force_calc)
            cdo.yseasstd(options='-f nc -b 32', output=sup_clim_std_file, input=file_monthly, force=force_calc)
            cdo.div(options='-f nc -b 32', output=sup_N_file, input=sup_sum_file + ' ' + sup_clim_file, force=force_calc)  # number of samples
        else:
            print interval
            raise ValueError('Unknown temporal interval. Can not perform preprocessing! ')

        if not os.path.exists(sup_clim_file):
            print 'File not existing (sup_clim_file): ' + sup_clim_file
            return None

        #3) read data
        sup = Data(sup_clim_file, 'rsus', read=True, label=self.model, unit='$W m^{-2}$', lat_name='lat', lon_name='lon', shift_lon=False)
        sup_std = Data(sup_clim_std_file, 'rsus', read=True, label=self.model + ' std', unit='-', lat_name='lat', lon_name='lon', shift_lon=False)
        sup.std = sup_std.data.copy()
        del sup_std
        sup_N = Data(sup_N_file, 'rsus', read=True, label=self.model + ' std', unit='-', lat_name='lat', lon_name='lon', shift_lon=False)
        sup.n = sup_N.data.copy()
        del sup_N

        # ensure that climatology always starts with January, therefore set date and then sort
        sup.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
        sup.timsort()

        #4) read monthly data
        supall = Data(file_monthly, 'rsus', read=True, label=self.model, unit='$W m^{-2}$', lat_name='lat', lon_name='lon', shift_lon=False)
        supall.adjust_time(day=15)
        if not supall._is_monthly():
            raise ValueError('Monthly timecycle expected here!')
        supmean = supall.fldmean()

        #/// return data as a tuple list
        retval = (supall.time, supmean, supall)
        del supall

        #/// mask areas without radiation (set to invalid): all data < 1 W/m**2
        #sup.data = np.ma.array(sis.data,mask=sis.data < 1.)

        return sup, retval

#-------------------------------------------------------------------------------------------------------------

    def get_albedo_data(self, interval='season', **kwargs):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        """

        force_calc = False

        # read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        # get fluxes
        Fu = self.get_surface_shortwave_radiation_up(interval=interval)
        if Fu is None:
            print 'File not existing for UPWARD flux!: ', self.name
            return None
        else:
            Fu_i = Fu[0]
        lab = Fu_i.label
        Fd = self.get_surface_shortwave_radiation_down(interval=interval, **{'CMIP5': {'valid_mask': 'land'}, 'CMIP5RAW': {'valid_mask': 'land'}})  # todo: take routine name from the configuration setup in JSON file !!!!
        if Fd is None:
            print 'File not existing for DOWNWARD flux!: ', self.name
            return None
        else:
            Fd_i = Fd[0]

        #albedo for chosen interval as caluclated as ratio of means of fluxes in that interval (e.g. season, months)
        Fu_i.div(Fd_i, copy=False)
        del Fd_i  # Fu contains now the albedo
        Fu_i._apply_mask(ls_mask.data)

        #albedo for monthly data (needed for global mean plots )
        Fu_m = Fu[1][2]
        del Fu
        Fd_m = Fd[1][2]
        del Fd

        Fu_m.div(Fd_m, copy=False)
        del Fd_m
        Fu_m._apply_mask(ls_mask.data)
        Fu_m._set_valid_range(0., 1.)
        Fu_m.label = lab + ' albedo'
        Fu_i.label = lab + ' albedo'
        Fu_m.unit = '-'
        Fu_i.unit = '-'

        # center dates of months
        Fu_m.adjust_time(day=15)
        Fu_i.adjust_time(day=15)

        # return data as a tuple list
        retval = (Fu_m.time, Fu_m.fldmean(), Fu_m)

        return Fu_i, retval


class CMIP5RAWData(CMIP5Data):
    """
    This class is supposed to use CMIP5 data in RAW format.
    This means that it builds on the CMORIZED CMIP5 data, but
    performs all necessary preprocessing step like e.g. calculation
    of ensemble means
    """
    def __init__(self, data_dir, model, experiment, dic_variables, name='', shift_lon=False, **kwargs):
        super(CMIP5RAWData, self).__init__(data_dir, model, experiment, dic_variables, name=model, shift_lon=shift_lon, **kwargs)
        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon
        self.type = 'CMIP5RAW'
        self._unique_name = self._get_unique_name()

    def _get_ensemble_filename(self, the_variable):
        """
        get filename of ensemble mean file
        if required, then all pre-processing steps are done

        Parameters
        ----------
        the_variable : str
            variable name to be processed

        Returns
        -------
        returns filename of file with multiensemble means
        """

        # use model parser to generate a list of available institutes and
        # models from data directory
        data_dir = self.data_dir
        if data_dir[-1] != os.sep:
            data_dir += os.sep

        CMP = preprocessor.CMIP5ModelParser(self.data_dir)
        model_list = CMP.get_all_models()

        # model name in configuration file is assumed to be INSTITUTE:MODEL
        institute = self.model.split(':')[0]
        model = self.model.split(':')[1]

        # TODO why is the institute not in the model output name ???
        output_file = get_temporary_directory() + the_variable + '_Amon_' + model + '_' + self.experiment + '_ensmean.nc'

        if institute not in model_list.keys():
            raise ValueError('Data for this institute is not existing: %s' % institute)

        # do preprocessing of data from multiple ensembles if file
        # already existing, then no processing is done
        C5PP = preprocessor.CMIP5Preprocessor(data_dir, output_file,
                                              the_variable, model,
                                              self.experiment,
                                              institute=institute)
        res_file = C5PP.ensemble_mean(delete=False,
                                      start_time=self.start_time,
                                      stop_time=self.stop_time)

        return res_file

class CMIP3Data(CMIP5Data):
    """
    Class for CMIP3 model simulations. This class is derived from C{Model}.
    """
    def __init__(self, data_dir, model, experiment, dic_variables, name='', shift_lon=False, **kwargs):
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
        super(CMIP3Data, self).__init__(data_dir, model, experiment, dic_variables, name=model, shift_lon=shift_lon, **kwargs)

        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon
        self.type = 'CMIP3'

        self._unique_name = self._get_unique_name()

