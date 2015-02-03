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
import ast
import numpy as np

from pycmbs.benchmarking import preprocessor
from pycmbs.benchmarking.utils import get_T63_landseamask, get_temporary_directory
from pycmbs.benchmarking.models.model_basic import *

from pycmbs.utils import print_log, WARNING


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
        if name == '':
            name = model
        super(CMIP5Data, self).__init__(data_dir, dic_variables, name=name, shift_lon=shift_lon, **kwargs)

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
        s = self.model.replace(' ', '') + '-' + self.experiment.replace(' ', '')
        s = s.replace('#', '-')
        if hasattr(self, 'ens_member'):
            s += '-' + str(self.ens_member)
        return s

    def get_rainfall_data(self, interval='season', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_wind(self, interval='season', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_evaporation(self, interval='season', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_latent_heat_flux(self, interval='season', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_model_data_generic(self, interval='season', **kwargs):
        """
        unique parameters are:
            filename - file basename
            variable - name of the variable as the short_name in the netcdf file

            kwargs is a dictionary with keys for each model. Then a dictionary with properties follows

        """

        if not self.type in kwargs.keys():
            print ''
            print 'WARNING: it is not possible to get data using generic function, as method missing: ', self.type, kwargs.keys()
            assert False


        locdict = kwargs[self.type]

        # read settings and details from the keyword arguments
        # no defaults; everything should be explicitely specified in either the config file or the dictionaries
        varname = locdict.pop('variable', None)
        #~ print self.type
        #~ print locdict.keys()
        assert varname is not None, 'ERROR: provide varname!'

        units = locdict.pop('unit', None)
        assert units is not None, 'ERROR: provide unit!'

        lat_name = locdict.pop('lat_name', 'lat')
        lon_name = locdict.pop('lon_name', 'lon')
        model_suffix = locdict.pop('model_suffix', None)
        model_prefix = locdict.pop('model_prefix', None)
        file_format = locdict.pop('file_format')
        scf = locdict.pop('scale_factor')
        valid_mask = locdict.pop('valid_mask')
        custom_path = locdict.pop('custom_path', None)
        thelevel = locdict.pop('level', None)

        target_grid = self._actplot_options['targetgrid']
        interpolation = self._actplot_options['interpolation']

        if custom_path is None:
            filename1 = self.get_raw_filename(varname, **kwargs)   # routine needs to be implemented by each subclass
        else:
            filename1 = custom_path + self.get_raw_filename(varname, **kwargs)

        if filename1 is None:
            print_log(WARNING, 'No valid model input data')
            return None

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
        mdata = Data(mdata_clim_file, varname, read=True, label=self._unique_name, unit=units, lat_name=lat_name, lon_name=lon_name, shift_lon=False, scale_factor=scf, level=thelevel, time_cycle=thetime_cylce)
        mdata_std = Data(mdata_clim_std_file, varname, read=True, label=self._unique_name + ' std', unit='-', lat_name=lat_name, lon_name=lon_name, shift_lon=False, level=thelevel, time_cycle=thetime_cylce)
        mdata.std = mdata_std.data.copy()
        del mdata_std
        mdata_N = Data(mdata_N_file, varname, read=True, label=self._unique_name + ' std', unit='-', lat_name=lat_name, lon_name=lon_name, shift_lon=False, scale_factor=scf, level=thelevel)
        mdata.n = mdata_N.data.copy()
        del mdata_N

        # ensure that climatology always starts with January, therefore set date and then sort
        mdata.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
        mdata.timsort()

        #4) read monthly data
        mdata_all = Data(file_monthly, varname, read=True, label=self._unique_name, unit=units, lat_name=lat_name, lon_name=lon_name, shift_lon=False, time_cycle=12, scale_factor=scf, level=thelevel)
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

        mdata._raw_filename = filename1
        mdata._monthly_filename = file_monthly
        mdata._clim_filename = mdata_clim_file
        mdata._varname = varname

        # return data as a tuple list
        retval = (mdata_all.time, mdata_mean, mdata_all)

        del mdata_all
        return mdata, retval


    def get_temperature_2m(self, interval='monthly', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_surface_shortwave_radiation_down(self, interval='monthly', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_surface_shortwave_radiation_up(self, interval='monthly', **kwargs):
        return self.get_model_data_generic(interval=interval, **kwargs)

    def get_albedo(self, interval='season', dic_up=None, dic_down=None):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,

        As the usage of different variables requires knowledge of the configuration of
        the input streams, these need to be provided in addition

        Parameters
        ----------
        dic_up : dict
            dictionary for get_surface_shortwave_radiation_up() as specified in model_data_routines.json
        dic_down : dict
            dictionary for get_surface_shortwave_radiation_down() as specified in model_data_routines.json
        """

        assert dic_up is not None, 'ERROR: dic_up needed'
        assert dic_down is not None, 'ERROR: dic_down needed'

        force_calc = False

        # read land-sea mask
        #~ ls_mask = get_T63_landseamask(self.shift_lon)

        #~ target grid  ??? valid mask ????

        def _extract_dict_from_routine_name(k, s):
            # extract dictionary name from routine name in model_data_routines.json
            res = ast.literal_eval(s[k].split('**')[1].rstrip()[:-1])
            #~ print res, type(res)
            return res

        # extract coniguration dictionaries for flues from model_data_routines
        kw_up = _extract_dict_from_routine_name('surface_upward_flux', dic_up)
        kw_down = _extract_dict_from_routine_name('sis', dic_down)

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        # get fluxes
        Fu = self.get_surface_shortwave_radiation_up(interval=interval, **kw_up)
        if Fu is None:
            print 'File not existing for UPWARD flux!: ', self.name
            return None
        else:
            Fu_i = Fu[0]
        lab = Fu_i._get_label()

        Fd = self.get_surface_shortwave_radiation_down(interval=interval, **kw_down)
        if Fd is None:
            print 'File not existing for DOWNWARD flux!: ', self.name
            return None
        else:
            Fd_i = Fd[0]

        # albedo for chosen interval as caluclated as ratio of means of fluxes in that interval (e.g. season, months)
        Fu_i.div(Fd_i, copy=False)
        del Fd_i  # Fu contains now the albedo
        #~ Fu_i._apply_mask(ls_mask.data)

        #albedo for monthly data (needed for global mean plots )
        Fu_m = Fu[1][2]
        del Fu
        Fd_m = Fd[1][2]
        del Fd

        Fu_m.div(Fd_m, copy=False)
        del Fd_m
        #~ Fu_m._apply_mask(ls_mask.data)
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
        super(CMIP5RAWData, self).__init__(data_dir, model, experiment, dic_variables, name=name, shift_lon=shift_lon, **kwargs)
        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon
        self.type = 'CMIP5RAW'
        self._unique_name = self._get_unique_name()

    def get_raw_filename(self, varname, **kwargs):
        mip = kwargs[self.type].pop('mip', None)
        assert mip is not None, 'ERROR: <mip> needs to be provided (CMIP5RAWSINGLE)'

        realm = kwargs[self.type].pop('realm')
        assert realm is not None, 'ERROR: <realm> needs to be provided (CMIP5RAWSINGLE)'

        return self._get_ensemble_filename(varname, mip, realm)

    def _get_ensemble_filename(self, the_variable, mip, realm):
        """
        get filename of ensemble mean file
        if required, then all pre-processing steps are done

        Parameters
        ----------
        the_variable : str
            variable name to be processed

        Returns
        -------
        returns filename of file with multi-ensemble means
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
        output_file = get_temporary_directory() + the_variable + '_' + mip + '_' + model + '_' + self.experiment + '_ensmean.nc'

        if institute not in model_list.keys():
            raise ValueError('Data for this institute is not existing: %s' % institute)

        # do preprocessing of data from multiple ensembles if file
        # already existing, then no processing is done
        C5PP = preprocessor.CMIP5Preprocessor(data_dir, output_file,
                                              the_variable, model,
                                              self.experiment,
                                              institute=institute, mip=mip, realm=realm)

        # calculate the ensemble mean and store as file
        # also the STDV is calculated on the fly calculated
        # resulting filenames are available by C5PP.outfile_ensmean and C5PP.outfile_ensstd
        C5PP.ensemble_mean(delete=False,
                                      start_time=self.start_time,
                                      stop_time=self.stop_time)

        return C5PP.outfile_ensmean


class CMIP5RAW_SINGLE(CMIP5RAWData):
    """
    This class is supposed to use CMIP5 data in RAW format.

    It is supposed to handle single emsemble members
    """

    def __init__(self, data_dir, model, experiment, dic_variables, name='', shift_lon=False, **kwargs):
        """
        Parameters
        ----------
        model_type : str
            model type like specified in the configuration file. It is
            supossed to be of format MPI-M:MPI-ESM-LR#1 etc.
            where after # there needs to be an integer number specifying
            the emsemble member number
        """

        if name == '':
            name = model

        # split between model type and ensemble member
        s = model.split('#')
        if len(s) != 2:
            print model, s
            raise ValueError('ERROR: invalid ensemble member specification')
        else:
            model = s[0]
            self.ens_member = int(s[1])

        self.institute = model.split(':')[0]

        super(CMIP5RAWData, self).__init__(data_dir, model, experiment, dic_variables, name=name, shift_lon=shift_lon, **kwargs)

        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon
        self.type = 'CMIP5RAWSINGLE'
        self._unique_name = self._get_unique_name()

    def get_raw_filename(self, variable, **kwargs):
        """
        return RAW filename for class CMIP5RAWSINGLE
        """

        # information comes from model_data_routines.json
        mip = kwargs[self.type].pop('mip', None)
        assert mip is not None, 'ERROR: <mip> needs to be provided (CMIP5RAWSINGLE)'

        realm = kwargs[self.type].pop('realm')
        assert realm is not None, 'ERROR: <realm> needs to be provided (CMIP5RAWSINGLE)'

        temporal_resolution = kwargs[self.type].pop('temporal_resolution')
        assert temporal_resolution is not None, 'ERROR: <temporal_resolution> needs to be provided (CMIP5RAWSINGLE)'

        data_dir = self.data_dir
        if data_dir[-1] != os.sep:
            data_dir += os.sep

        model = self.model.split(':')[1]
        fp = data_dir + self.institute + os.sep + model + os.sep + self.experiment + os.sep + temporal_resolution + os.sep + realm + os.sep + mip + os.sep + 'r' + str(self.ens_member) + 'i1p1' + os.sep + variable + os.sep + variable + '_' + mip + '_' + model + '_' + self.experiment + '_r' + str(self.ens_member) + 'i1p1_*.nc'

        files = glob.glob(fp)
        if len(files) == 0:
            return None
        if len(files) != 1:
            print files
            raise ValueError('More than one file found!')
        return files[0]


class CMIP3Data(CMIP5Data):
    """
    Class for CMIP3 model simulations. This class is derived from C{Model}.
    """
    def __init__(self, data_dir, model, experiment, dic_variables, name='', shift_lon=False, **kwargs):
        """

        Parameters
        ----------
        data_dir: directory that specifies the root directory where the data is located
        model: TBD tood
        experiment: specifies the ID of the experiment (str)
        dic_variables:
        name: TBD todo
        shift_lon: specifies if longitudes of data need to be shifted
        kwargs: other keyword arguments
        """
        super(CMIP3Data, self).__init__(data_dir, model, experiment, dic_variables, name=model, shift_lon=shift_lon, **kwargs)

        self.model = model
        self.experiment = experiment
        self.data_dir = data_dir
        self.shift_lon = shift_lon
        self.type = 'CMIP3'

        self._unique_name = self._get_unique_name()
