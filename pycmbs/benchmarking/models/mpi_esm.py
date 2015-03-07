# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
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


class JSBACH_BOT(Model):

    def __init__(self, filename, dic_variables, experiment, name='', shift_lon=False, **kwargs):
        super(JSBACH_BOT, self).__init__(filename, dic_variables, name=name, **kwargs)

        self.experiment = experiment
        self.shift_lon = shift_lon
        self.type = 'JSBACH_BOT'

        self._unique_name = self._get_unique_name()

    def _get_unique_name(self):
        """
        get unique name from model and experiment
        @return: string with unique combination of models and experiment
        """
        return self.name.replace(' ', '') + '-' + self.experiment.replace(' ', '')

    def get_albedo_data(self, interval='season'):
        """
        get albedo data for JSBACH

        returns Data object
        """

        if interval != 'season':
            raise ValueError('Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry')

        v = 'var176'

        filename = self.data_dir + 'data/model1/' + self.experiment + '_echam6_BOT_mm_1979-2006_albedo_yseasmean.nc'
        ls_mask = get_T63_landseamask(self.shift_lon)

        albedo = Data(filename, v, read=True,
                      label='MPI-ESM albedo ' + self.experiment, unit='-', lat_name='lat', lon_name='lon',
                      shift_lon=self.shift_lon,
                      mask=ls_mask.data.data)

        return albedo

    def get_tree_fraction(self, interval='season'):
        """
        todo implement this for data from a real run !!!
        """

        if interval != 'season':
            raise ValueError('Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry')

        ls_mask = get_T63_landseamask(self.shift_lon)

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_forest_shrub.nc'
        v = 'var12'
        tree = Data(filename, v, read=True,
                    label='MPI-ESM tree fraction ' + self.experiment, unit='-', lat_name='lat', lon_name='lon',
                    shift_lon=self.shift_lon,
                    mask=ls_mask.data.data, start_time=pl.num2date(pl.datestr2num('2001-01-01')), stop_time=pl.num2date(pl.datestr2num('2001-12-31')))

        return tree

    def get_grass_fraction(self, interval='season'):
        """
        todo implement this for data from a real run !!!
        """

        if interval != 'season':
            raise ValueError('Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry')

        ls_mask = get_T63_landseamask(self.shift_lon)

        filename = '/home/m300028/shared/dev/svn/trstools-0.0.1/lib/python/pyCMBS/framework/external/vegetation_benchmarking/VEGETATION_COVER_BENCHMARKING/example/historical_r1i1p1-LR_1850-2005_grass_crop_pasture_2001.nc'
        v = 'var12'
        grass = Data(filename, v, read=True,
                     label='MPI-ESM tree fraction ' + self.experiment, unit='-', lat_name='lat', lon_name='lon',
                     #shift_lon=shift_lon,
                     mask=ls_mask.data.data, start_time=pl.num2date(pl.datestr2num('2001-01-01')), stop_time=pl.num2date(pl.datestr2num('2001-12-31')), squeeze=True)

        return grass

    def get_surface_shortwave_radiation_down(self, interval='season'):
        """
        get surface shortwave incoming radiation data for JSBACH

        returns Data object
        """

        if interval != 'season':
            raise ValueError('Other temporal sampling than SEASON not supported yet for JSBACH BOT files, sorry')

        v = 'var176'

        y1 = '1979-01-01'
        y2 = '2006-12-31'
        rawfilename = self.data_dir + 'data/model/' + self.experiment + '_echam6_BOT_mm_1979-2006_srads.nc'

        if not os.path.exists(rawfilename):
            return None

        #--- read data
        cdo = pyCDO(rawfilename, y1, y2)
        if interval == 'season':
            seasfile = cdo.seasmean()
            del cdo
            print 'seasfile: ', seasfile
            cdo = pyCDO(seasfile, y1, y2)
            filename = cdo.yseasmean()
        else:
            raise ValueError('Invalid interval option %s ' % interval)

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SIS data
        sis = Data(filename, v, read=True,
                   label='MPI-ESM SIS ' + self.experiment, unit='-', lat_name='lat', lon_name='lon',
                   #shift_lon=shift_lon,
                   mask=ls_mask.data.data)

        return sis

    def get_rainfall_data(self, interval='season'):
        """
        get rainfall data for JSBACH
        returns Data object
        """

        if interval == 'season':
            pass
        else:
            raise ValueError('Invalid value for interval: %s' % interval)

        #/// PREPROCESSING: seasonal means ///
        s_start_time = str(self.start_time)[0:10]
        s_stop_time = str(self.stop_time)[0:10]

        filename1 = self.data_dir + self.experiment + '_echam6_BOT_mm_1980_sel.nc'
        tmp = pyCDO(filename1, s_start_time, s_stop_time).seldate()
        tmp1 = pyCDO(tmp, s_start_time, s_stop_time).seasmean()
        filename = pyCDO(tmp1, s_start_time, s_stop_time).yseasmean()

        #/// READ DATA ///

        #1) land / sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #2) precipitation data
        try:
            v = 'var4'
            rain = Data(filename, v, read=True, scale_factor=86400.,
                        label='MPI-ESM ' + self.experiment, unit='mm/day', lat_name='lat', lon_name='lon',
                        shift_lon=self.shift_lon,
                        mask=ls_mask.data.data)
        except:
            v = 'var142'
            rain = Data(filename, v, read=True, scale_factor=86400.,
                        label='MPI-ESM ' + self.experiment, unit='mm/day', lat_name='lat', lon_name='lon',
                        shift_lon=self.shift_lon,
                        mask=ls_mask.data.data)

        return rain


class JSBACH_RAW2(Model):
    """
    Class for RAW JSBACH model output
    works on the real raw output
    """

    #def __init__(self, filename, dic_variables, experiment, name='', shift_lon=False, model_dict=None, input_format='grb', raw_outdata='outdata/jsbach/', **kwargs):
    def __init__(self, filename, dic_variables, experiment, name='', shift_lon=False, input_format='grb', raw_outdata='outdata/jsbach/', **kwargs):
        """

        The assignment of certain variables to different input streams is done in the routine
        get_jsbach_data_generic()


        Parameters
        ----------
        input_format : str
            specifies file format of input data
            ['nc','grb']
        """

        super(JSBACH_RAW2, self).__init__(filename, dic_variables, name=name, **kwargs)

        self.experiment = experiment
        self.shift_lon = shift_lon
        #self.get_data()
        self.type = 'JSBACH_RAW2'

        self.input_format = input_format
        assert self.input_format in ['nc', 'grb']

        self.raw_outdata = raw_outdata

        self._unique_name = self._get_unique_name()

        # do preprocessing of streams (only needed once!) ---
        self.files = {}
        self._preproc_streams()
        #~ self.model_dict = copy.deepcopy(model_dict)

        self.model = 'JSBACH'

    def _get_filenames_jsbach_stream(self):
        return self.data_dir + self.raw_outdata + self.experiment + '_jsbach_main_mm_*.' + self.input_format

    def _get_filenames_veg_stream(self):
        return self.data_dir + self.raw_outdata + self.experiment + '_jsbach_veg_mm_*.' + self.input_format

    def _get_filenames_land_stream(self):
        return self.data_dir + self.raw_outdata + self.experiment + '_jsbach_land_mm_*.' + self.input_format

    def _get_filenames_surf_stream(self):
        return self.data_dir + self.raw_outdata + self.experiment + '_jsbach_surf_mm_*.' + self.input_format

    def _get_filenames_albedo_VIS(self):
        return self.data_dir + self.raw_outdata + self.experiment + '_jsbach_mm_*_VIS_albedo.' + self.input_format

    def _get_filenames_albedo_NIR(self):
        return self.data_dir + self.raw_outdata + self.experiment + '_jsbach_mm_*_NIR_albedo.' + self.input_format

    def _get_filenames_echam_BOT(self):
        return self.data_dir + self.raw_outdata + '../echam6/' + self.experiment + '_echam6_BOT_mm_*.sz'

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

        # jsbach stream
        print '   JSBACH stream ...'
        outfile = get_temporary_directory() + self.experiment + '_jsbach_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach.codes'
            tmp = tempfile.mktemp(suffix='.nc', prefix=self.experiment + '_jsbach_', dir=get_temporary_directory())  # temporary file
            #~ print self.data_dir
            #~ print self.raw_outdata
            #~ print 'Files: ', self._get_filenames_jsbach_stream()
            #~ stop
            if len(glob.glob(self._get_filenames_jsbach_stream())) > 0:  # check if input files existing at all
                print 'Mering the following files:', self._get_filenames_jsbach_stream()
                cdo.mergetime(options='-f nc', output=tmp, input=self._get_filenames_jsbach_stream())
                if os.path.exists(codetable):
                    cdo.monmean(options='-f nc', output=outfile, input='-setpartab,' + codetable + ' ' + tmp)  # monmean needed here, as otherwise interface does not work
                else:
                    cdo.monmean(options='-f nc', output=outfile, input=tmp)  # monmean needed here, as otherwise interface does not work
                print 'Outfile: ', outfile
                #~ os.remove(tmp)

                print 'Temporary name: ', tmp

        self.files.update({'jsbach': outfile})

        # veg stream
        print '   VEG stream ...'
        outfile = get_temporary_directory() + self.experiment + '_jsbach_veg_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach_veg.codes'
            tmp = tempfile.mktemp(suffix='.nc', prefix=self.experiment + '_jsbach_veg_', dir=get_temporary_directory())  # temporary file
            if len(glob.glob(self._get_filenames_veg_stream())) > 0:  # check if input files existing at all
                cdo.mergetime(options='-f nc', output=tmp, input=self._get_filenames_veg_stream())
                if os.path.exists(codetable):
                    cdo.monmean(options='-f nc', output=outfile, input='-setpartab,' + codetable + ' ' + tmp)  # monmean needed here, as otherwise interface does not work
                else:
                    cdo.monmean(options='-f nc', output=outfile, input=tmp)  # monmean needed here, as otherwise interface does not work
                os.remove(tmp)
        self.files.update({'veg': outfile})

        # veg land
        print '   LAND stream ...'
        outfile = get_temporary_directory() + self.experiment + '_jsbach_land_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach_land.codes'
            tmp = tempfile.mktemp(suffix='.nc', prefix=self.experiment + '_jsbach_land_', dir=get_temporary_directory())  # temporary file
            if len(glob.glob(self._get_filenames_land_stream())) > 0:  # check if input files existing at all
                cdo.mergetime(options='-f nc', output=tmp, input=self._get_filenames_land_stream())
                if os.path.exists(codetable):
                    cdo.monmean(options='-f nc', output=outfile, input='-setpartab,' + codetable + ' ' + tmp)  # monmean needed here, as otherwise interface does not work
                else:
                    cdo.monmean(options='-f nc', output=outfile, input=tmp)  # monmean needed here, as otherwise interface does not work
                os.remove(tmp)
        self.files.update({'land': outfile})

        # surf stream
        print '   SURF stream ...'
        outfile = get_temporary_directory() + self.experiment + '_jsbach_surf_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_jsbach_surf.codes'
            tmp = tempfile.mktemp(suffix='.nc', prefix=self.experiment + '_jsbach_surf_', dir=get_temporary_directory())  # temporary file
            if len(glob.glob(self._get_filenames_surf_stream())) > 0:  # check if input files existing at all
                print glob.glob(self._get_filenames_surf_stream())
                cdo.mergetime(options='-f nc', output=tmp, input=self._get_filenames_surf_stream())
                if os.path.exists(codetable):
                    cdo.monmean(options='-f nc', output=outfile, input='-setpartab,' + codetable + ' ' + tmp)  # monmean needed here, as otherwise interface does not work
                else:
                    cdo.monmean(options='-f nc', output=outfile, input=tmp)  # monmean needed here, as otherwise interface does not work
                os.remove(tmp)

        self.files.update({'surf': outfile})

        # ECHAM BOT stream
        print '   BOT stream ...'
        outfile = get_temporary_directory() + self.experiment + '_echam6_echam_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            codetable = self.data_dir + 'log/' + self.experiment + '_echam6_echam.codes'
            tmp = tempfile.mktemp(suffix='.nc', prefix=self.experiment + '_echam6_echam_', dir=get_temporary_directory())  # temporary file
            if len(glob.glob(self._get_filenames_echam_BOT())) > 0:  # check if input files existing at all
                cdo.mergetime(options='-f nc', output=tmp, input=self._get_filenames_echam_BOT())
                if os.path.exists(codetable):
                    cdo.monmean(options='-f nc', output=outfile, input='-setpartab,' + codetable + ' ' + tmp)  # monmean needed here, as otherwise interface does not work
                else:
                    cdo.monmean(options='-f nc', output=outfile, input=tmp)  # monmean needed here, as otherwise interface does not work
                os.remove(tmp)
        self.files.update({'echam': outfile})

        # ALBEDO file
        # albedo files as preprocessed by a script of Thomas
        print '   ALBEDO VIS stream ...'
        outfile = get_temporary_directory() + self.experiment + '_jsbach_VIS_albedo_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            if len(glob.glob(self._get_filenames_albedo_VIS())) > 0:  # check if input files existing at all
                cdo.mergetime(options='-f nc', output=outfile, input=self._get_filenames_albedo_VIS())
        self.files.update({'albedo_vis': outfile})

        print '   ALBEDO NIR stream ...'
        outfile = get_temporary_directory() + self.experiment + '_jsbach_NIR_albedo_mm_full.nc'
        if os.path.exists(outfile):
            pass
        else:
            if len(glob.glob(self._get_filenames_albedo_NIR())) > 0:  # check if input files existing at all
                cdo.mergetime(options='-f nc', output=outfile, input=self._get_filenames_albedo_NIR())
        self.files.update({'albedo_nir': outfile})

    def _get_unique_name(self):
        """
        get unique name from model and experiment
        @return: string with unique combination of models and experiment
        """
        return self.name.replace(' ', '') + '-' + self.experiment.replace(' ', '')

    def get_albedo_data(self, interval='season'):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,

        This routine uses the definitions of the routines how to
        read upward and downward fluxes
        """

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        #~ tmpdict = copy.deepcopy(kwargs)
        #~ print self.dic_vars

        routine_up = self.dic_vars['surface_upward_flux']
        routine_down = self.dic_vars['sis']

        #sw_down = self.get_surface_shortwave_radiation_down(interval=interval, **kwargs)
        cmd = 'sw_down = self.' + routine_down
        exec(cmd)

        #sw_up = self.get_surface_shortwave_radiation_up(interval=interval, **kwargs)
        cmd = 'sw_up = self.' + routine_up
        exec(cmd)

        # climatological mean
        alb = sw_up[0].div(sw_down[0])
        alb.label = self.experiment + ' albedo'
        alb.unit = '-'

        # original data
        alb_org = sw_up[1][2].div(sw_down[1][2])
        alb_org.label = self.experiment + ' albedo'
        alb_org.unit = '-'

        retval = (alb_org.time, alb_org.fldmean(), alb_org)

        return alb, retval

    def get_albedo_data_vis(self, interval='season', **kwargs):
        """
        This routine retrieves the JSBACH albedo information for VIS
        it requires a preprocessing with a script that aggregates from tile
        to box values!

        Parameters
        ----------
        interval : str
            ['season','monthly']
        """
        #~ tmpdict = copy.deepcopy(self.model_dict['albedo_vis'])
        return self.get_jsbach_data_generic(interval=interval, **kwargs)

    def get_albedo_data_nir(self, interval='season', **kwargs):
        """
        This routine retrieves the JSBACH albedo information for VIS
        it requires a preprocessing with a script that aggregates from tile
        to box values!

        Parameters
        ----------
        interval : str
            ['season','monthly']
        """
        #~ tmpdict = copy.deepcopy(self.model_dict['albedo_nir'])
        return self.get_jsbach_data_generic(interval=interval, **kwargs)

    def get_surface_shortwave_radiation_up(self, interval='season', **kwargs):
        return self.get_jsbach_data_generic(interval=interval, **kwargs)

    def get_surface_shortwave_radiation_down(self, interval='season', **kwargs):
        return self.get_jsbach_data_generic(interval=interval, **kwargs)

    def get_rainfall_data(self, interval='season', **kwargs):
        return self.get_jsbach_data_generic(interval=interval, **kwargs)

    def get_temperature_2m(self, interval='season', **kwargs):
        return self.get_jsbach_data_generic(interval=interval, **kwargs)

    def get_jsbach_data_generic(self, interval='season', **kwargs):
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

        varname = locdict.pop('variable')
        units = locdict.pop('unit', 'Unit not specified')

        lat_name = locdict.pop('lat_name', 'lat')
        lon_name = locdict.pop('lon_name', 'lon')
        #model_suffix = locdict.pop('model_suffix')
        #model_prefix = locdict.pop('model_prefix')
        file_format = locdict.pop('file_format')
        scf = locdict.pop('scale_factor')
        valid_mask = locdict.pop('valid_mask')
        custom_path = locdict.pop('custom_path', None)
        thelevel = locdict.pop('level', None)

        target_grid = self._actplot_options['targetgrid']
        interpolation = self._actplot_options['interpolation']

        if self.type != 'JSBACH_RAW2':
            print self.type
            raise ValueError('Invalid data format here!')

        # define from which stream of JSBACH data needs to be taken for specific variables
        if varname in ['swdown_acc', 'swdown_reflect_acc']:
            filename1 = self.files['jsbach']
        elif varname in ['precip_acc']:
            filename1 = self.files['land']
        elif varname in ['temp2']:
            filename1 = self.files['echam']
        elif varname in ['var14']:  # albedo vis
            filename1 = self.files['albedo_vis']
        elif varname in ['var15']:  # albedo NIR
            filename1 = self.files['albedo_nir']
        else:
            print varname
            raise ValueError('Unknown variable type for JSBACH_RAW2 processing!')

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
            raise ValueError('Unknown temporal interval. Can not perform preprocessing! ')

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

        #ensure that climatology always starts with J  anuary, therefore set date and then sort
        mdata.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
        mdata.timsort()

        #4) read monthly data
        mdata_all = Data(file_monthly, varname, read=True, label=self.model, unit=units, lat_name=lat_name, lon_name=lon_name, shift_lon=False, time_cycle=12, scale_factor=scf, level=thelevel)
        mdata_all.adjust_time(day=15)

        if target_grid == 't63grid':
            mdata._apply_mask(get_T63_landseamask(False, area=valid_mask))
            mdata_all._apply_mask(get_T63_landseamask(False, area=valid_mask))
        else:
            tmpmsk = get_generic_landseamask(False, area=valid_mask, target_grid=target_grid)
            mdata._apply_mask(tmpmsk)
            mdata_all._apply_mask(tmpmsk)
            del tmpmsk

        mdata_mean = mdata_all.fldmean()

        # return data as a tuple list
        retval = (mdata_all.time, mdata_mean, mdata_all)

        del mdata_all

        return mdata, retval


class JSBACH_SPECIAL(JSBACH_RAW2):
    """
    special class for more flexible reading of JSBACH input data
    it allows to specify the input format and the directory of the input data

    in case that you use a different setup, it is probably easiest to
    just copy this class and make the required adaptations.
    """
    def __init__(self, filename, dic_variables, experiment, name='', shift_lon=False, model_dict=None, input_format='nc', raw_outdata='', **kwargs):
        super(JSBACH_SPECIAL, self).__init__(filename, dic_variables, experiment, name=name, shift_lon=shift_lon, model_dict=model_dict, input_format=input_format, raw_outdata=raw_outdata, **kwargs)


class xxxxxxxxJSBACH_RAW(Model):
    """
    Class for RAW JSBACH model output
    works on manually preprocessed already concatenated data
    """

    def __init__(self, filename, dic_variables, experiment, name='', shift_lon=False, intervals='monthly', **kwargs):
        super(JSBACH_RAW, self).__init__(filename, dic_variables, name=name, intervals=intervals, **kwargs)

        print('WARNING: This model class should be depreciated as it contained a lot of hardcoded dependencies and is only intermediate')
        #TODO: depreciate this class
        stop

        self.experiment = experiment
        self.shift_lon = shift_lon
        self.type = 'JSBACH_RAW'
        self._unique_name = self._get_unique_name()

    def _get_unique_name(self):
        """
        get unique name from model and experiment
        """
        return self.name.replace(' ', '') + '-' + self.experiment.replace(' ', '')

    def get_temperature_2m(self, interval='monthly', **kwargs):
        """
        get surface temperature (2m) from JSBACH model results

        Parameters
        ----------
        interval : str
            specifies the aggregation interval. Possible options: ['season','monthly']
        """

        locdict = kwargs[self.type]

        y1 = '1980-01-01'  # TODO move this to the JSON dictionary or some parameter file
        y2 = '2010-12-31'
        variable = 'temp2'
        rawfile = self.data_dir + self.experiment + '_echam6_echam_' + variable + '_ALL.nc'
        files = glob.glob(rawfile)
        if len(files) != 1:
            print 'Inputfiles: ', files
            raise ValueError('Something went wrong: Invalid number of input files!')
        else:
            rawfile = files[0]
        mdata, retval = self._do_preprocessing(rawfile, variable, y1, y2, interval=interval, valid_mask=locdict['valid_mask'])
        return mdata, retval

    def get_albedo_data(self, interval='monthly', **kwargs):
        """
        calculate albedo as ratio of upward and downwelling fluxes
        first the monthly mean fluxes are used to calculate the albedo,
        """

        # read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)  # TODO make this more flexible

        if self.start_time is None:
            raise ValueError('Start time needs to be specified')
        if self.stop_time is None:
            raise ValueError('Stop time needs to be specified')

        Fd = self.get_surface_shortwave_radiation_down(**kwargs)
        Fu = self.get_surface_shortwave_radiation_up(**kwargs)

        if Fu is None:
            print 'File not existing for UPWARD flux!: ', self.name
            return None
        else:
            Fu_i = Fu[0]
            if Fu_i is None:
                return None

        if Fd is None:
            print 'File not existing for DOWNWARD flux!: ', self.name
            return None
        else:
            Fd_i = Fd[0]
            if Fd_i is None:
                return None
        lab = Fu_i.label

        # albedo for chosen interval as caluclated as ratio of means of fluxes in that interval (e.g. season, months)
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

#-----------------------------------------------------------------------

    def _do_preprocessing(self, rawfile, varname, s_start_time, s_stop_time, interval='monthly', force_calc=False, valid_mask='global', target_grid='t63grid'):
        """
        perform preprocessing
        * selection of variable
        * temporal subsetting
        """
        cdo = Cdo()

        if not os.path.exists(rawfile):
            print('File not existing! %s ' % rawfile)
            return None, None

        # calculate monthly means
        file_monthly = get_temporary_directory() + os.sep + os.path.basename(rawfile[:-3]) + '_' + varname + '_' + s_start_time + '_' + s_stop_time + '_mm.nc'
        if (force_calc) or (not os.path.exists(file_monthly)):
            cdo.monmean(options='-f nc', output=file_monthly, input='-seldate,' + s_start_time + ',' + s_stop_time + ' ' + '-selvar,' + varname + ' ' + rawfile, force=force_calc)
        else:
            pass
        if not os.path.exists(file_monthly):
            raise ValueError('Monthly preprocessing did not work! %s ' % file_monthly)

        # calculate monthly or seasonal climatology
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

        # read data
        if interval == 'monthly':
            thetime_cylce = 12
        elif interval == 'season':
            thetime_cylce = 4
        else:
            print interval
            raise ValueError('Unsupported interval!')

        mdata = Data(mdata_clim_file, varname, read=True, label=self.name, shift_lon=False, time_cycle=thetime_cylce, lat_name='lat', lon_name='lon')
        mdata_std = Data(mdata_clim_std_file, varname, read=True, label=self.name + ' std', unit='-', shift_lon=False, time_cycle=thetime_cylce, lat_name='lat', lon_name='lon')
        mdata.std = mdata_std.data.copy()
        del mdata_std
        mdata_N = Data(mdata_N_file, varname, read=True, label=self.name + ' std', shift_lon=False, lat_name='lat', lon_name='lon')
        mdata.n = mdata_N.data.copy()
        del mdata_N

        # ensure that climatology always starts with January, therefore set date and then sort
        mdata.adjust_time(year=1700, day=15)  # set arbitrary time for climatology
        mdata.timsort()

        #4) read monthly data
        mdata_all = Data(file_monthly, varname, read=True, label=self.name, shift_lon=False, time_cycle=12, lat_name='lat', lon_name='lon')
        mdata_all.adjust_time(day=15)

        #mask_antarctica masks everything below 60 degree S.
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

    def get_surface_shortwave_radiation_down(self, interval='monthly', **kwargs):
        """
        get surface shortwave incoming radiation data for JSBACH

        Parameters
        ----------
        interval : str
            specifies the aggregation interval. Possible options: ['season','monthly']
        """

        locdict = kwargs[self.type]

        y1 = '1980-01-01'  # TODO move this to the JSON dictionary or some parameter file
        y2 = '2010-12-31'
        rawfile = self.data_dir + self.experiment + '_jsbach_' + y1[0: 4] + '_' + y2[0: 4] + '.nc'
        mdata, retval = self._do_preprocessing(rawfile, 'swdown_acc', y1, y2, interval=interval, valid_mask=locdict['valid_mask'])
        return mdata, retval

#-----------------------------------------------------------------------

    def get_surface_shortwave_radiation_up(self, interval='monthly', **kwargs):
        """
        get surface shortwave upward radiation data for JSBACH

        Parameters
        ----------
        interval : str
            specifies the aggregation interval. Possible options: ['season','monthly']
        """

        locdict = kwargs[self.type]

        y1 = '1980-01-01'  # TODO: move this to the JSON dictionary or some parameter file
        y2 = '2010-12-31'
        rawfile = self.data_dir + self.experiment + '_jsbach_' + y1[0: 4] + '_' + y2[0: 4] + '.nc'
        mdata, retval = self._do_preprocessing(rawfile, 'swdown_reflect_acc', y1, y2, interval=interval, valid_mask=locdict['valid_mask'])
        return mdata, retval

#-----------------------------------------------------------------------

    def get_model_data_generic(self, interval='monthly', **kwargs):
        """
        This is only a wrapper to redirect to individual functions
        for the JSBACH_RAW class

        Currently only the usage for rainfall is supported!
        """
        # HACK: only a wrapper, should be depreciated
        raise ValueError('Rainfall analysis not working yet!')
        self.get_rainfall_data(interval=interval, **kwargs)

    def get_rainfall_data(self, interval='monthly', **kwargs):
        """
        get surface rainfall data for JSBACH
        uses already preprocessed data where the convective and
        advective rainfall has been merged

        Parameters
        ----------
        interval : str
            specifies the aggregation interval. Possible options: ['season','monthly']
        """

        locdict = kwargs[self.type]

        y1 = '1980-01-01'  # TODO : move this to the JSON dictionary or some parameter file
        y2 = '2010-12-31'
        variable = 'aprc'
        rawfile = self.data_dir + self.experiment + '_echam6_echam_*_precipitation.nc'
        files = glob.glob(rawfile)
        if len(files) != 1:
            print 'Inputfiles: ', files
            raise ValueError('Something went wrong: Invalid number of input files!')
        else:
            rawfile = files[0]

        mdata, retval = self._do_preprocessing(rawfile, variable, y1, y2, interval=interval, valid_mask=locdict['valid_mask'])
        return mdata, retval

#-----------------------------------------------------------------------

    def get_gpp_data(self, interval='season'):
        """
        get surface GPP data for JSBACH

        todo temporal aggregation of data --> or leave it to the user!
        """
        cdo = Cdo()
        v = 'var167'
        y1 = str(self.start_time)[0:10]
        y2 = str(self.stop_time)[0:10]
        rawfilename = self.data_dir + 'data/model/' + self.experiment + '_' + y1[0:4] + '-' + y2[0:4] + '.nc'
        times_in_file = int(''.join(cdo.ntime(input=rawfilename)))

        if interval == 'season':
            if times_in_file != 4:
                tmp_file = get_temporary_directory() + os.path.basename(rawfilename)
                cdo.yseasmean(options='-f nc -b 32 -r ', input='-selvar,' + v + ' ' + rawfilename, output=tmp_file[:-3] + '_yseasmean.nc')
                rawfilename = tmp_file[:-3] + '_yseasmean.nc'

        if interval == 'monthly':
            if times_in_file != 12:
                tmp_file = get_temporary_directory() + os.path.basename(rawfilename)
                cdo.ymonmean(options='-f nc -b 32 -r ', input='-selvar,' + v + ' ' + rawfilename, output=tmp_file[:-3] + '_ymonmean.nc')
                rawfilename = tmp_file[:-3] + '_ymonmean.nc'

        if not os.path.exists(rawfilename):
            return None

        filename = rawfilename

        #--- read land-sea mask
        ls_mask = get_T63_landseamask(self.shift_lon)

        #--- read SW up data
        gpp = Data4D(filename, v, read=True,
                     label=self.experiment + ' ' + v, unit='gC m-2 a-1', lat_name='lat', lon_name='lon',
                     shift_lon=self.shift_lon,
                     mask=ls_mask.data.data, scale_factor=3600. * 24. * 30. / 0.083
                     )

        return gpp.sum_data4D()

#-----------------------------------------------------------------------
