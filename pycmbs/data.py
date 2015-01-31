# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os
import sys

from pycmbs.statistic import get_significance, ttest_ind
from pycmbs.netcdf import NetCDFHandler
from pycmbs.polygon import Raster
from pycmbs.polygon import Polygon as pycmbsPolygon

import numpy as np
from matplotlib import pylab as plt
import matplotlib.pylab as pl
from scipy import stats
from netCDF4 import netcdftime

from calendar import monthrange
from cdo import Cdo
import datetime
import pytz
import pickle
import datetime
import calendar
import tempfile
import struct
import gzip


class Data(object):

    """
    Data class: main class
    """

    def __init__(self, filename, varname, lat_name=None, lon_name=None,
                 read=False, scale_factor=1., label=None, unit=None,
                 shift_lon=False, start_time=None, stop_time=None,
                 mask=None, time_cycle=None, squeeze=False, level=None,
                 verbose=False, cell_area=None, time_var='time',
                 checklat=True, weighting_type='valid', oldtime=False,
                 warnings=True, calc_cell_area=True):
        """
        Constructor for Data class

        Parameters
        ----------
        filename : str
            name of the file that contains the data
            (specify None if not data from file)

        varname : str
            name of the variable that contains the data
            (specify None if not data from file)

        lat_name : str
            name of latitude field

        lon_name : str
            name of longitude field

        read : bool
            specifies if the data should be read directly when
            creating the Data object

        scale_factor : float
            scale factor to rescale the data after reading. Be aware,
            that this IS NOT the scale factor that is used for data
            compression in e.g. netCDF variables
            (variable attribute scale_factor) but this variable has
            the purpose to perform a rescaling (e.g. change of units)
            of the data immediately after it was read (e.g. convert
            rainfall intensity [mm/h] to [mm/day], the scale_factor
            would be 1./24.

        label : str
            the label of the data. This is automatically used
            e.g. in plotting routines

        unit : str
            specify the unit of the data. Will be used for plotting.

        shift_lon : bool
            if this flag is True, then the longitudes are shifted to
            [-180 ... 180] if they are given in [0 ... 360]

        start_time : datetime object
            specifies the start time of the data to be read from
            file (if None=default, then beginning of file is used)

        stop_time : datetime object
            specifies the stop time of the data to be read from
            file (if None, then end of file is used)

        mask : array(ny,nx)
            mask that is applied to the data when it is read.
            Needs to have same geometry of data.

        time_cycle : int
            specifies size of the periodic cycle of the data.
            This is needed if deseasonalized anomalies should be
            calculated. If one works for instance with monthly
            data, time_cycle needs to be 12. Assumption is that the
            temporal sampling of the data is regular. When reading
            from netCDF file, the regularity is checked.

        squeeze : bool
            remove singletone dimensions in the data when it is read.

        level : int
            specify level to work with (needed for 4D data)

        verbose : bool
            verbose mode for printing

        cell_area : numpy array (ny,nx)
            area size [m**2] of each grid cell. These weights are uses
            as an INITIAL weight and are renormalized in accordance to
            the valid data only.

        warnings : bool
            log warnings in a separate logfile

        weighting_type : str
           specifies how normalization shall be done ['valid','all']
           'valid': the weights are calculated based on all VALID
           values, thus the sum of all these weights is one
           'all': contrary, weights are calculated based on ALL
           (valid and invalid) data. The latter option can be useful,
           if one is interested e.g. in a global area weighted mean,
           whereas the values in 'self' are only valid for e.g.
           land areas. If one wants to calculate e.e. the change in
           global mean surface fluxes, given only land fluxes, one can
           normalize using normtype='all' and then gets the impact on
           the global mean fluxes. In the other case (normtype='valid'),
           one would get the change in the global mean of the LAND
           fluxes only!

        oldtime : bool
            if True, then the old definition for time is used, which is
            compliant with the pylab time definition, as *x* is a float
            value which gives the number of days (fraction part
            represents hours, minutes, seconds) since
            0001-01-01 00:00:00 UTC *plus* *one*.

            NOTE that there is a *PLUS ONE*. The actual data object
            supports different calendars, while this is not possible
            using the pylab num2date/date2num functions. (Un)fortunately,
            the new routine, which is implemented in self.num2date(),
            self.date2num() *performs correctly* the calculations. Thus
            there is *NO* *PLUS ONE* needed. As a consequence, all files
            that have been written with older versions of pyCMBS have a
            wrong time variable included in the file. To allow for
            backwards compliance, the option oldtime=True can be used
            which will then mimic a similar behaviour as the pylab
            date functions.
        calc_cell_area : bool
                calculate cell area
        """

        self.weighting_type = weighting_type
        self.filename = filename
        self.varname = varname
        self.scale_factor = scale_factor
        self.lat_name = lat_name
        self.lon_name = lon_name
        self.squeeze = squeeze
        self.squeezed = False
        self.detrended = False
        self.verbose = verbose
        # assume that coordinates are always in 0 < lon < 360
        self._lon360 = True
        self._calc_cell_area = calc_cell_area

        self.inmask = mask
        self.level = level
        self.gridtype = None
        self._oldtime = oldtime
        # specifies if latitudes have been checked for increasing order
        #(required for zonal plot)
        self._latitudecheckok = False
        self.warnings = warnings

        if label is None:
            if self.filename is None:
                self.label = ''
            else:
                self.label = self.filename
        else:
            self.label = label

        if unit is None:
            self.unit = None
        else:
            self.unit = unit

        if time_cycle is not None:
            self.time_cycle = time_cycle

        self.cell_area = cell_area  # [m**2]

        self.lat = None
        self.lon = None

        #/// read data from file ///
        if read:
            self.read(shift_lon, start_time=start_time, stop_time=stop_time,
                      time_var=time_var, checklat=checklat)

        # check if longitudes are from 0 ... 360
        if self.lon is not None:
            if self._lon360:
                if self.lon.min() < 0.:
                    # self._shift_lon_360() #shift coordinates to [0 ... 360.]
                    # here we assume that the coordinates are between -180 ...
                    # 180
                    self._lon360 = False
                if self.lon.max() > 360.:
                    print self.lon.max()
                    raise ValueError('invalid longitudes needs shifting !!!')
            else:
                print '_lon360: ', self._lon360
                print self.lon.min(), self.lon.max()
                print 'WARNING: plotting etc not supported for longitudes which are not equal to 0 ... 360'

    def _get_shape(self):
        return self.data.shape
    shape = property(_get_shape)

    def _get_date(self):
        #--- convert to datetime objects ---
        # use this approach to ensure that a datetime.datetime array is available for further processing
        # set also timezone as UTC as otherwise comparisons of dates is not possible!
        # CAUTION: assumes that timezone is always UTC !!

        try:
            return np.asarray(
                [datetime.datetime(x.year, x.month, x.day, x.hour, x.minute, x.second, 0, pytz.UTC) for x in
                 self.num2date(self.time)])
        # if an exception occurs then write data on screen for bughandling
        except:
            if os.path.exists('dump.pkl'):
                os.remove('dump.pkl')
            pickle.dump(self, open('dump.pkl', 'w'))
            print 'An error occured in Data.date! Printing data for bugfixing'
            for i in range(len(self.time)):
                x = self.num2date(self.time[i])
                print self.time[i], x, datetime.datetime(
                    x.year, x.month, x.day, x.hour, x.minute, x.second, 0,
                    pytz.UTC)
            raise ValueError(
                'Some error in time conversion happened! Look in dump.pkl to fix it')

    date = property(_get_date)

    def _get_ndim(self):
        return self.data.ndim
    ndim = property(_get_ndim)

    def _get_nt(self):
        return len(self.time)
    nt = property(_get_nt)

    def _get_ny(self):
        if self.ndim == 2:
            ny, nx = self.data.shape
        elif self.ndim == 3:
            nt, ny, nx = self.data.shape
        else:
            raise ValueError('Invalid geometry!')
        return ny
    ny = property(_get_ny)

    def _get_nx(self):
        if self.ndim == 2:
            ny, nx = self.data.shape
        elif self.ndim == 3:
            nt, ny, nx = self.data.shape
        else:
            raise ValueError('Invalid geometry!')
        return nx
    nx = property(_get_nx)

    def _get_mindate(self, base=None):
        """
        Brief
        -----
        get minimum date

        Parameters
        ----------
        base : str
            ['day','month']; if given, then e.g. the first of the month
            or the first time of the day is returned instead of the
            actual minimum value. This allows to easily round the mindate.

        Test
        ----
        unittest implemented
        """
        dmin = self.date.min()

        if base is None:
            rval = dmin
        elif base == 'month':
            rval = datetime.datetime(
                dmin.year, dmin.month, 1, 0, 0, 0, 0, dmin.tzinfo)
        elif base == 'day':
            rval = datetime.datetime(
                dmin.year, dmin.month, dmin.day, 0, 0, 0, 0, dmin.tzinfo)
        elif base == 'year':
            rval = datetime.datetime(dmin.year, 1, 1, 0, 0, 0, 0, dmin.tzinfo)
        else:
            raise ValueError('Invalid base in _get_mindate()')
        return rval

    def _get_maxdate(self, base=None):
        """
        get maximum date

        Parameters
        ----------
        base : str
            ['day','month']; if given, then e.g. the last day of the month or the last time of the day is
            returned instead of the actual maximum value; this allows to easily round the maxdate

        Test
        ----
        unittest implemented
        """

        dmax = self.date.max()

        if base is None:
            rval = dmax
        elif base == 'month':
            rval = datetime.datetime(dmax.year, dmax.month,
                                     calendar.monthrange(dmax.year, dmax.month)[1], 23, 59, 59, 0, dmax.tzinfo)
        elif base == 'day':
            rval = datetime.datetime(
                dmax.year, dmax.month, dmax.day, 23, 59, 59, 0, dmax.tzinfo)
        elif base == 'year':
            rval = datetime.datetime(
                dmax.year, 12, 31, 23, 59, 59, 0, dmax.tzinfo)
        else:
            raise ValueError('Invalid base in _get_maxdate()')
        return rval

    def _days_per_month(self):
        """return the number of days per month in Data timeseries (unittest)"""
        return [float(calendar.monthrange(d.year, d.month)[1]) for d in self.date]

    def _log_warning(self, s, write_log=False):
        """
        log warnings for class in a logfile

        Parameters
        ----------
        s : str
            string with warning message
        write_log : bool
            do actual data loging (default = False to avoid unnecessary generation of log file)
        """

        if not write_log:
            return

        if 'DATA_WARNING_FILE' in os.environ.keys():
            file = os.environ['DATA_WARNING_FILE']
            # create output directory if necessary
            if not os.path.exists(os.path.dirname(file)):
                os.makedirs(os.path.dirname(file))
        else:
            file = 'data_warnings.log'

        if os.path.exists(file):
            mode = 'a'
        else:
            mode = 'w'

        print("   " + s)

        f = open(file, mode)
        if self.filename is None:
            filename = 'data object has not filename'
        else:
            filename = self.filename

        f.write(filename + '\t' + s + '\n')
        f.close()

    def _oldtimeoffset(self):
        """
        return offset to convert to old time
        offset is one day *PLUS ONE* following pylab documentation
        This routine takes care of different time units
        """

        print('This routine is depreciated!')  # TODO

        if not hasattr(self, 'time_str'):
            raise ValueError('ERROR: time offset can not be determined!')

        if 'hours' in self.time_str:
            return 1. * 24.
        elif 'seconds' in self.time_str:
            return 1. * 86400.
        elif 'days' in self.time_str:
            return 1.
        else:
            print self.time_str
            raise ValueError(
                'ERROR: Invalid timestring: conversion not possible!')

    def num2date(self, t):
        """
        convert a numeric time to a datetime object
        Routine is similar to pylab num2date, but allows to make use
        of different calendars. It encapsulates the corresponding
        netcdftime function

        Parameters
        ----------
        t : float
            numerical value describing time in accordance with calendar
            settings

        Returns
        -------
        d : datetime
        """

        if self._oldtime:  # see documentation in __init__ of self
            offset = self._oldtimeoffset()
        else:
            offset = 0.
        if not hasattr(self, 'time_str'):
            raise ValueError('num2date can not work without timestr!')
        if self.time_str is None:
            raise ValueError('num2date can not work without timestr!')
        else:
            return netcdftime.num2date(t + offset, self.time_str,
                                       calendar=self.calendar)

    def date2num(self, t):
        """
        convert a datetime object into a numeric time variable
        Routine is similar to pylab date2num, but allows to make use
        of different calendars. It encapsulates the corresponding
        netcdftime function

        Parameters
        ----------
        t : datetime
            datetime object with current date

        Returns
        -------
        d : float
            numerical value describing time in accordance with calendar
            settings
        """

        if self._oldtime:  # see documentation in __init__ of self
            offset = self._oldtimeoffset()
        else:
            offset = 0.
        if not hasattr(self, 'time_str'):
            raise ValueError('date2num can not work without timestr!')
        if self.time_str is None:
            raise ValueError('date2num can not work without timestr!')
        else:
            return netcdftime.date2num(t, self.time_str,
                                       calendar=self.calendar) - offset

    def save(self, filename, varname=None, format='nc',
             delete=False, mean=False, timmean=False, compress=True):
        """
        saves the data object to a file

        Parameters
        ----------
        filename : str
            filename the file should be saved too
        varname : str
            variable name in output file. If *None*, then
            the variables are just named like var001 ...
        format : str
            output format ['nc','ascii','nc3','nc4']
        delete : bool
            delete file if existing without asking. If *False*, and the
            file is existing already, then an error is raised
        mean : bool
            save spatial mean field only instead of the full field
        timmean : bool
            save temporal mean field
        compress : bool
            compress resulting file if supported by library
        """

        map_formats = {'nc': 'NETCDF4', 'nc3':
                       'NETCDF3_CLASSIC', 'nc4': 'NETCDF4'}

        if mean and timmean:
            raise ValueError(
                'Only the MEAN or the TIMMEAN option can be given, but not together!')

        # either store full field or just spatial mean field
        if mean:
            print('Saving MEAN FIELD of object in file %s' % filename)
            tmp = self.fldmean(return_data=True)
        elif timmean:
            tmp = self.timmean(return_object=True)
        else:
            print('Saving object in file %s' % filename)
            tmp = self

        # store data now ...
        if format in ['nc', 'nc3', 'nc4']:
            tmp._save_netcdf(
                filename, varname=varname, delete=delete, compress=compress,
                format=map_formats[format])
        elif format == 'ascii':
            tmp._save_ascii(filename, varname=varname, delete=delete)
        else:
            raise ValueError('This output format is not defined yet!')

    def _save_ascii(self, filename, varname=None, delete=False):
        """
        saves the data object to an ASCII file as follows

        time  lon lat value

        Parameters
        ----------
        filename : str
            filename of output file
        varname : str
            name of output variable; this explicitely overwrites
            self.varname, which is tried to be used in first place
        delete : bool
            delete file if existing without further asking

        Test
        ----
        unittest implemented
        """

        #/// check if output file already there
        if os.path.exists(filename):
            if delete:
                os.remove(filename)
            else:
                raise ValueError(
                    'File already existing. Please delete manually or use DELETE option: %s' % filename)

        if hasattr(self, 'time'):
            if self.time is None:
                notime = True
            else:
                notime = False
        else:
            notime = True

        F = open(filename, 'w')
        if notime:
            if self.ndim == 2:
                F.write(self._arr2string(self.data, prefix=''))
            else:
                raise ValueError(
                    'Saving ASCII not implemented for data without time yet!')
        else:

            if self.ndim == 2:  # temporal mean field assumed
                F.write(self._arr2string(self.data, prefix='timmean'))
            elif self.ndim == 3:
                for i in xrange(len(self.time)):
                    F.write(
                        self._arr2string(self.data[i, :, :], prefix=str(self.date[i])))
            else:
                raise ValueError('Invalid geometry!')

        F.close()

    def _arr2string(self, a, prefix='', sep='\t'):
        """
        convert a 2D numpy array to an ASCII list
        this routine is supposed to be used to convert the 2D field
        of a timestep and return a list as follows

        prefix  lon1  lat1  value1
        prefix  lon2  lat2  value2
        ...

        Parameters
        ----------
        a : ndarray (2D)
            numpy array with data; needs to have geometry ny, nx
        prefix : str
            prefix to be appended
        """

        ny, nx = a.shape

        assert (self.ny == ny)
        assert (self.nx == nx)

        s = ''
        sep = '\t'
        eol = '\n'

        if len(prefix) > 0:
            prefix += sep

        for i in xrange(ny):
            for j in xrange(nx):
                if a.mask[i, j]:  # in case of masked values
                    pass
                else:
                    s += prefix + \
                        str(self.lon[i, j]) + sep + \
                        str(self.lat[i, j]) + sep + str(a[i, j]) + eol

        return s

    def _save_netcdf(self, filename, varname=None, delete=False, compress=True, format='NETCDF4'):
        """
        saves the data object to a netCDF file

        Parameters
        ----------
        filename : str
            filename of output file
        varname : str
            name of output variable; this explicitely overwrites
            self.varname, which is tried to be used as a first order
        delete : bool
            delete file if existing without asking
        compress : bool
            compress resulting file if supported by backend
        format : str
            output file format specifier as used by netCDF4 library
        """

        # check if output file already there
        if os.path.exists(filename):
            if delete:
                os.remove(filename)
            else:
                raise ValueError('File already existing. Please delete \
                                  manually or use DELETE option: \
                                  %s' % filename)
        # variable name
        if varname is None:
            if self.varname is None:
                varname = 'var1'
            else:
                varname = self.varname

        # create new file
        File = NetCDFHandler()
        File.open_file(filename, 'w', format=format)

        # create dimensions
        if hasattr(self, 'data'):
            if self.data.ndim == 3:
                if self.time is None:
                    raise ValueError('No time variable existing! \
                                      Can not write 3D data!')
                nt, ny, nx = self.shape
            elif self.data.ndim == 2:
                ny, nx = self.shape
            else:
                raise ValueError('Current shape not supported here %s' %
                                 self.shape)
        else:
            ny, nx = np.shape(self.lat)

        # Create Dimensions
        File.create_dimension('ny', size=ny)
        File.create_dimension('nx', size=nx)

        # Create variables
        if hasattr(self, 'time'):
            if self.time is not None:
                File.create_dimension('time', size=len(self.time))
                File.create_variable('time', 'd', ('time',))
                File.F.variables['time'].units = self.time_str

        if hasattr(self, 'data'):
            if self.data.ndim == 3:
                File.create_variable(
                    varname, 'd', ('time', 'ny', 'nx'), zlib=compress)
            elif self.data.ndim == 2:
                File.create_variable(varname, 'd', ('ny', 'nx'), zlib=compress)

        if self.lat is not None:
            File.create_variable('lat', 'd', ('ny', 'nx'))
            File.F.variables['lat'].units = 'degrees_north'
            File.F.variables['lat'].axis = "Y"
            File.F.variables['lat'].long_name = "latitude"

        if self.lon is not None:
            File.create_variable('lon', 'd', ('ny', 'nx'))
            File.F.variables['lon'].units = 'degrees_east'
            File.F.variables['lon'].axis = "X"
            File.F.variables['lon'].long_name = "longitude"

        if hasattr(self, 'data'):
            if hasattr(self, 'cell_area'):
                File.create_variable('cell_area', 'd', ('ny', 'nx'))

        #/// write data
        if hasattr(self, 'time'):
            if self.time is not None:
                File.assign_value('time', self.time)
                File.F.variables['time'].calendar = self.calendar

        if hasattr(self, 'data'):
            File.assign_value(varname, self.data)

        if self.lat is not None:
            File.assign_value('lat', self.lat)
        if self.lon is not None:
            File.assign_value('lon', self.lon)

        if hasattr(self, 'data'):
            if hasattr(self, 'cell_area'):
                if self.cell_area is not None:
                    File.assign_value('cell_area', self.cell_area)

        if hasattr(self, 'data'):
            if hasattr(self, 'long_name'):
                if self.long_name is not None:
                    File.F.variables[varname].long_name = self.long_name
            if hasattr(self, 'unit'):
                if self.unit is not None:
                    File.F.variables[varname].units = self.unit

            File.F.variables[varname].scale_factor = 1.
            File.F.variables[varname].add_offset = 0.
            File.F.variables[varname].coordinates = "lon lat"

        File.close()

    def _equal_lon(self):
        """
        This routine identifies if all longitudes in the dataset
        are the same (except for numerical uncertainties)

        Returns
        -------
        bool

        Test
        ----
        unittest implemented
        """
        if self.lon.ndim == 1:
            return True
        elif self.lon.ndim == 2:
            lu = self.lon.mean(axis=0)
            # this corresponds to an accuracy of 1m
            if any(np.abs(lu - self.lon[0, :]) > 1.E-5):
                return False
            else:
                return True
        else:
            raise ValueError('Unsupported geometry for longitude')

    def _get_unique_lon(self):
        """
        estimate if the Data contains unique longitudes and if so,
        returns a vector with these longitudes

        Returns
        -------
        lon : ndarray
            unique longitudes

        Test
        ----
        unittest implemented
        """

        if self.lon is None:
            raise ValueError(
                'Can not estimate longitude, as no longitudes existing!')

        if self.lon.ndim == 1:
            return self.lon
        elif self.lon.ndim == 2:
            if self._equal_lon():  # ... check for unique lons
                return self.lon[0, :]
            else:
                print self.filename
                raise ValueError(
                    'The dataset does not contain unique LONGITUDES!')
        else:
            raise ValueError(
                'Data dimension for longitudes not supported yet!')

    def get_bounding_box(self):
        """
        estimates bounding box of valid data. It returns the indices
        of the bounding box which frames all valid data

        CAUTION
        note that the indices can not be used directly for array slicing!
        One typically needs to add '1' to the last index

        Returns
        -------
        r : list
            returns indices for bounding box [i1,i2,j1,j2]

        Test
        ----
        unittest implemented
        """

        # return mask where all timesteps are valid!
        msk = self.get_valid_mask()  # gives a 2D mask

        # estimate boundary box indices
        xb = msk.sum(axis=0)
        yb = msk.sum(axis=1)

        j1 = -99
        for i in xrange(len(xb)):
            if j1 > 0:
                continue
            if (xb[i] > 0) & (j1 == -99):
                j1 = i

        j2 = -99
        for i in xrange(len(xb) - 1, -1, -1):
            if j2 > 0:
                continue
            if (xb[i] > 0) & (j2 == -99):
                j2 = i

        i1 = -99
        for i in xrange(len(yb)):
            if i1 > 0:
                continue
            if (yb[i] > 0) & (i1 == -99):
                i1 = i

        i2 = -99
        for i in xrange(len(yb) - 1, -1, -1):
            if i2 > 0:
                continue
            if (yb[i] > 0) & (i2 == -99):
                i2 = i
        return i1, i2, j1, j2

    def _squeeze(self):
        """
        remove singletone dimensions in data variable
        """

        if self.verbose:
            print 'SQUEEZING data ... ', self.data.ndim, self.data.shape

        if self.data.ndim > 2:
            self.data = self.data.squeeze()
            self.squeezed = True

        if self.verbose:
            print 'AFTER SQUEEZING data ... ', self.data.ndim, self.data.shape

    def _set_cell_area(self):
        """
        set cell area size. If a cell area was already given (either by user or from file)
        nothing will happen. Otherwise it will be tried to calculate cell_area from
        coordinates using the CDO's

        The estimation of the cell area follows the following steps:
        1) try to estimate cellarea using cdo gridarea
        2) if this does not work:
            a) directory is write protected --> write to temporary directory
            b) unknown grid --> try to select another grid, using cdo selgrid
        3) it could be that the cdo's dont support the input filetype. In that case
           it is tried to read the lat/lon information using the netCDF4 library
           store these fields in a dummy nc3 file and then apply the cdo's
        4) if all this does not work, then cell_area is set to unity for all grid cells and a WARNING is raised
        """

        # TODO unittest implementation

        if hasattr(self, 'cell_area'):
            if self.cell_area is not None:
                return

        # in case that cell area shall explicitely NOT be caluclated
        if not self._calc_cell_area:
            if self.ndim == 2:
                self.cell_area = np.ones(self.data.shape)
            elif self.ndim == 3:
                self.cell_area = np.ones(self.data[0, :, :].shape)
            else:
                raise ValueError('Invalid geometry!')
            return

        if (self.lat is None) or (self.lon is None):
            self._log_warning(
                "WARNING: cell area can not be calculated (missing coordinates)!")
            if self.ndim == 2:
                self.cell_area = np.ones(self.data.shape)
            elif self.ndim == 3:
                self.cell_area = np.ones(self.data[0, :, :].shape)
            else:
                raise ValueError('Invalid geometry!')
            return

        # calculate cell area from coordinates
        cell_file = self.filename[:-3] + '_cell_area.nc'

        if not os.path.exists(cell_file):  # calculate grid area using CDO's
            cdo = Cdo()
            try:
                cdo.gridarea(options='-f nc', output=cell_file,
                             input=self.filename)
            except:
                # occurs if you dont have write permissions
                print '   Seems that cell_area file can not be generated, try to generate in temporary directory'
                # generate some temporary filename
                cell_file = tempfile.mktemp(prefix='cell_area_', suffix='.nc')
                try:
                    cdo.gridarea(options='-f nc', output=cell_file,
                                 input=self.filename)
                    print '   Cell area file generated sucessfully in temporary file: ' + cell_file
                except:
                    # not sucessfull so far ... last try here by selecting an
                    # alternative grid (if available)
                    print(
                        "   Try to calculate gridarea using alternative grid")
                    # generate some temporary filename
                    cell_file = tempfile.mktemp(
                        prefix='cell_area_', suffix='.nc')
                    try:
                        cdo.gridarea(options='-f nc', output=cell_file,
                                     input='-selgrid,2 ' + self.filename)
                        print '   Cell area file generated sucessfully in temporary file: ' + cell_file
                    except:
                        try:
                            # store lat/lon coordinates in nc3 file and then
                            # apply the cdo's (last try)
                            coord_file = tempfile.mktemp(
                                prefix='coordinates_', suffix='.nc')
                            tmpdat = Data(None, None)
                            tmpdat.lat = self.lat
                            tmpdat.lon = self.lon
                            tmpxxx = np.zeros(self.lat.shape)
                            tmpdat.data = np.ma.array(
                                tmpxxx, mask=tmpxxx != tmpxxx)
                            tmpdat.save(coord_file, format='nc3')
                            del tmpxxx, tmpdat
                            cell_file = tempfile.mktemp(
                                prefix='cell_area_', suffix='.nc')
                            cdo.gridarea(options='-f nc',
                                         output=cell_file, input=coord_file)
                        except:
                            print('WARNING: no cell area could be generated!')

        # read cell_area file ---
        if os.path.exists(cell_file):
            # read cell area from file
            File = NetCDFHandler()
            File.open_file(cell_file, 'r')
            self.cell_area = File.get_variable('cell_area')
            File.close()

            # check geometries
            if self.data.ndim == 2:
                if self.cell_area.shape != self.data.shape:
                    raise ValueError(
                        'Invalid cell_area file: delete it manually and check again!')
            elif self.data.ndim == 1:
                if self.cell_area.shape != self.data.shape:
                    raise ValueError(
                        'Invalid cell_area file: delete it manually and check again!')
            elif self.data.ndim == 3:
                if self.cell_area.shape != self.data[0, :, :].shape:
                    raise ValueError(
                        'Invalid cell_area file: delete it manually and check again!')
        else:
            # no cell area calculation possible!!!
            # logger.warning('Can not estimate cell area! (setting all equal) ' + cell_file)

            self._log_warning(
                '*** WARNING: Can not estimate cell area! ' + cell_file)
            self._log_warning(' setting cell_area all to equal')

            if self.data.ndim == 2:
                self.cell_area = np.ones(self.data.shape)
            elif self.data.ndim == 3:
                self.cell_area = np.ones(self.data[0, :, :].shape)
            else:
                print 'actual geometry:  ', self.data.ndim, self.data.shape
                raise ValueError('Invalid geometry!')

    def get_zonal_mean(self, return_object=False):
        """
        calculate zonal mean statistics of the data for each timestep
        returns zonal statistics [time,ny]

        uses area weighting of data
        gives exact same results as function 'zonmean' in cdo's

        Parameters
        ----------
        return_object : bool
            if True, then returns a Data object

        Returns
        -------
        r : ndarray, Data
            array with zonal statistics
        """

        # TODO implement check if latitudes in y-axis direction are all the
        # same! Otherwise the routine does not make sense

        if self.cell_area is None:
            self._log_warning(
                'WARNING: no cell area given, zonal means are based on equal weighting!')
            w = np.ones(self.data.shape)
        else:
            w = self._get_weighting_matrix()

        #/// weight data
        dat = self.data * w

        #/// calculate zonal mean
        if dat.ndim == 2:
            r = dat.sum(axis=1) / w.sum(axis=1)  # zonal mean
        elif dat.ndim == 3:
            nt, ny, nx = dat.shape
            r = np.ones((nt, ny)) * np.nan
            W = np.ones((nt, ny)) * np.nan

            for i in xrange(nt):
                # weighted sum, normalized by valid data why ???
                r[i] = dat[i, :, :].sum(axis=1) / w[i, :, :].sum(axis=1)
                W[i] = w[i, :, :].sum(axis=1)
            r = np.ma.array(r, mask=W == 0.)

        else:
            print dat.shape
            raise ValueError('Unsupported geometry')

        if return_object:
            res = self.copy()
            res.label = self.label + ' zonal mean'
            res.data = r.T  # [lat,time]
            res.lat = self.lat[:, 0]  # latitudes as a vector
        else:
            res = r
        return res

    def get_percentile(self, p, return_object=True):
        """
        calculate percentile

        Parameters
        ----------
        p : float
            percentile value to obtain, e.g. 0.05 corresponds to 5% percentil
        return_object : bool
            specifies of a C{Data} object shall be returned [True] or a numpy array [False]

        Returns
        -------
        r : ndarray, Data
            returns the percentiles as either C{Data} object or as numpy array
        """

        if self.data.ndim != 3:
            raise ValueError(
                'Percentile calculation only supported for 3D data!')

        nt = len(self.data)
        x = self.data.copy()
        x.shape = (nt, -1)

        # calculate percentile
        res = stats.mstats.scoreatpercentile(x, p * 100.)

        # reshape data array
        res.shape = np.shape(self.data[0, :, :])
        res = np.ma.array(res, mask=np.isnan(res))

        # return
        if return_object:
            r = self.copy()
            r.label = self.label + '\n percentile: ' + str(round(p, 2))
            r.data = res
            return r
        else:
            return res

    def _get_unit(self):
        """
        get a nice looking string for units like e.g. '[mm/d]'
        """
        if self.unit is None:
            u = ''
        else:
            u = '[' + self.unit + ']'
        return u

    def _shift_lon(self):
        """
        shift longitude coordinates. Coordinates given as [0...360] are
        converted to [-180...180]

        changes lon field of Data object and sets variable _lon360
        """
        self.lon[self.lon >= 180.] -= 360.
        self._lon360 = False

    def _shift_lon_360(self):
        """
        shift longitude coordinates. Coordinates given as [-180...180] are
        converted to [0 ...360]

        changes lon field of Data object and sets variable _lon360
        """
        self.lon[self.lon < 0.] += 360.
        self._lon360 = True
        print('Longitudes were shifted to 0 ... 360!')

    def _apply_temporal_mask(self, mask):
        """
        apply a temporal mask to data. All timesteps where the mask is
        True will be masked, but geometry will not be changed, thus no
        masking will be applied

        The Data.data is changed

        Parameters
        ----------
        mask : ndarray of type bool
            needs to be of size [time]

        Returns
        -------
        None
        """

        if self.data.ndim != 3:
            raise ValueError('temporal masking only possible for 3D data')

        if len(mask) != len(self.data):
            print len(mask), self.data.shape
            raise ValueError('Inconsistent length of data and mask')

        for i in xrange(len(mask)):
            if mask[i]:
                self.data.mask[i, :, :] = True

    def read(self, shift_lon, start_time=None, stop_time=None,
             time_var='time', checklat=True, fmt='nc'):
        """
        Read data from a file. This functions provides a wrapper for
        I/O of different file formats

        Parameters
        ----------
        shift_lon : bool
            if given, longitudes will be shifted
        start_time : datetime
            start time for reading the data
        stop_time : datetime
            stop time for reading the data
        time_var : str
            name of time variable field
        checklat : bool
            check if latitude is in decreasing order (N ... S)
        fmt : str
            format of data to read
            ['nc']
        """
        if not os.path.exists(self.filename):
            raise ValueError('Error: file not existing: %s' % self.filename)

        self.time_var = time_var

        netcdf_backend = 'netCDF4'
        # read data
        if fmt == 'nc':
            self.data = self.read_netcdf(
                self.varname, netcdf_backend=netcdf_backend)
        else:
            raise ValueError('ERROR: invalid input format')

        if self.data is None:
            print self.varname
            raise ValueError(
                'The data in the file %s is not existing. This must not happen!' % self.filename)
        if self.scale_factor is None:
            raise ValueError(
                'The scale_factor for file %s is NONE, this must not happen!' % self.filename)

        # ensure that no Nan values occur
        np.ma.masked_where(np.isnan(self.data), self.data, copy=False)

        # this scaling is related to unit conversion and NOT
        # due to data compression
        self.data *= self.scale_factor

        # squeeze data to singletone
        if self.squeeze:
            self._squeeze()

        # mask data when desired ---
        if self.inmask is not None:
            self._apply_mask(self.inmask)

        # read lat/lon (try default names)
        if self.lat_name is None:
            self.lat_name = 'lat'
        if self.lon_name is None:
            self.lon_name = 'lon'
        if self.lat_name is not None:
            self.lat = self.read_netcdf(
                self.lat_name, netcdf_backend=netcdf_backend)
            # ensure that lat has NOT dimension (1,nlat)
            if self.lat is not None:
                if self.lat.ndim == 2:
                    if self.lat.shape[0] == 1:
                        self.lat = self.lat[0, :]
        else:
            self.lat = None

        if not self.lon_name is None:
            self.lon = self.read_netcdf(
                self.lon_name, netcdf_backend=netcdf_backend)
            # ensure that lon has NOT dimension (1,nlon)
            if self.lon is not None:
                if self.lon.ndim == 2:
                    if self.lon.shape[0] == 1:
                        self.lon = self.lon[0, :]
                # shift longitudes such that -180 < lon < 180
            if shift_lon:
                self._shift_lon()
        else:
            self.lon = None

        if self.lat is None:
            print('*** WARNING!!! No coordinates available!')

        # read time
        if self.time_var is not None:
            # returns either None or a masked array
            self.time = self.read_netcdf(
                self.time_var, netcdf_backend=netcdf_backend)
            if hasattr(self.time, 'mask'):
                self.time = self.time.data
            else:
                self.time = None
            if self.time is not None:
                if self.time.ndim != 1:
                    self.time = self.time.flatten()  # remove singletone
        else:
            self.time = None

        # determine time
        if self.time is not None:
            self.set_time()

        # lat lon to 2D matrix
        try:
            self._mesh_lat_lon()
        except:
            if self.verbose:
                print '        WARNING: No lat/lon mesh was generated!'

        # check if cell_area is already existing. if not,
        # try to calculate from coordinates
        self._set_cell_area()

        #  check if latitude in decreasing order (N ... S)?
        if checklat:
            if hasattr(self, 'lat'):
                if self.lat is not None:
                    # increasing order!
                    if np.all(np.diff(self.lat[:, 0]) > 0.):
                        self._flipud()
                        self._latitudecheckok = True
                    # decreasing order!
                    elif np.all(np.diff(self.lat[:, 0]) < 0.):
                        self._latitudecheckok = True
                    else:
                        print 'WARNING: latitudes not in systematic order! Might cause trouble with zonal statistics!'
                        self._latitudecheckok = False

        # calculate climatology from ORIGINAL (full dataset)
        if hasattr(self, 'time_cycle'):
            if self.time is not None:
                self._climatology_raw = self.get_climatology()

        # perform temporal subsetting
        if self.time is not None:

            # no temporal subsetting for 2D data! --> results in invalid
            # results!
            if self.ndim == 3:
                #- now perform temporal subsetting
                # BEFORE the conversion to the right time is required!
                m1, m2 = self._get_time_indices(start_time, stop_time)
                self._temporal_subsetting(m1, m2)

        # calculate time_cycle automatically if not set already.
        if self.time is not None:
            if hasattr(self, 'time_cycle'):
                if self.time_cycle is None:
                    self._set_timecycle()
            else:
                self._set_timecycle()

    def _get_binary_filehandler(self, mode='r'):
        """
        get filehandler for binary file
        it allows to access also gzip files

        Parameters
        ----------
        mode : str
            ['w','r']

        Returns
        -------
        file handler
        """
        mode += 'b'
        if self.filename[-3:] == '.gz':
            f = gzip.open(self.filename, mode)
        else:
            f = open(self.filename, mode)
        return f

    def _read_binary_file(self, dtype=None, lat=None, lon=None,
                          lonmin=None, lonmax=None, latmin=None,
                          latmax=None, ny=None, nx=None, nt=None):
        """
        read data from binary file
        this routine also allows spatial subsetting during reading
        from the binary file

        It requires however that two vectors of lon/lat are provided
        in case that a subsetting shall be made

        Parameters
        ----------
        dtype : str
            datatype specification
        lat : ndarray
            vector of latitudes
        lon : ndarray
            vector of longitudes
        """
        if dtype is None:
            raise ValueError('ERROR: dtype not provided')
        if lat is not None:
            assert (lon is not None)
        if lat is None:
            assert (ny is not None)
            assert (nx is not None)
        else:
            assert (lat.ndim == 1)
            assert (lon.ndim == 1)

        # specify bytes as well as format string for struct for each datatype
        dtype_spec = {
            'int16': 'H',
            'double': 'd',
            'int32': 'i',
            'float': 'f'
        }

        if dtype not in dtype_spec.keys():
            raise ValueError('ERROR: invalid data type')

        #~ 'int8'   :'b',
        #~ 'uint8'  :'B',
        #~ 'int16'  :'h',
        #~ 'uint16' :'H',
        #~ 'int32'  :'i',
        #~ 'uint32' :'I',
        #~ 'int64'  :'q',
        #~ 'uint64' :'Q',
        #~ 'float'  :'f',
        #~ 'double' :'d',
        #~ 'char'   :'s'}

        # set boundaries
        if lon is not None:
            if lonmin is None:
                lonmin = lon.min()
            if lonmax is None:
                lonmax = lon.max()
        else:
            lonmin = None
            lonmax = None

        if lat is not None:
            if latmin is None:
                latmin = lat.min()
                latmax = lat.max()
        else:
            latmin = None
            latmax = None

        # get file handler
        f = self._get_binary_filehandler()

        # read actual data
        if lon is None:
            file_content = f.read()  # read entire file
            # TODO if specifie, then read lat/lon information from file
            self.lat = None
            self.lon = None
        else:
            # check if lat/lon is increasing
            assert np.all(np.diff(lat) > 0.)
            assert np.all(np.diff(lon) > 0.)

            lonminpos = np.abs(lon - lonmin).argmin()
            if lon[lonminpos] < lonmin:
                lonminpos += 1
            lonmaxpos = np.abs(lon - lonmax).argmin()
            if lon[lonmaxpos] > lonmax:
                lonmaxpos -= 1

            latminpos = np.abs(lat - latmin).argmin()
            if lat[latminpos] < latmin:
                latminpos += 1

            latmaxpos = np.abs(lat - latmax).argmin()
            if lat[latmaxpos] > latmax:
                latmaxpos -= 1

            olon = lon[lonminpos:lonmaxpos + 1]
            olat = lat[latminpos:latmaxpos + 1]

            ny = len(olat)
            nx = len(olon)

            self.lon, self.lat = np.meshgrid(olon, olat)

            print 'coordinates: ', lonmin, lonmax, latmin, latmax
            print 'Positions: ', lonminpos, lonmaxpos, latminpos, latmaxpos

            file_content = self._read_binary_subset2D(f, struct.calcsize(dtype_spec[dtype]), ny=len(
                lat), nx=len(lon), xbeg=lonminpos, xend=lonmaxpos + 1, ybeg=latminpos, yend=latmaxpos + 1)

        # close file
        f.close()

        # put data into final matrix by decoding in accordance to the filetype

        self.data = np.reshape(
            np.asarray(struct.unpack(dtype_spec[dtype] * ny * nx * nt, file_content)), (ny, nx))

        del file_content

    def _read_binary_subset2D(self, f, nbytes, xbeg=None, xend=None, ybeg=None, yend=None, ny=None, nx=None):
        """
        read subset from binary file
        deconding takes place outside of this routine!

        Parameters
        ----------
        f : file handler
            file handler to read binary data
        nbytes : int
            size of a single value on the file [bytes]
        ybeg : int
            index of y coordinate start
        yend : int
            index of y coordinate end
        xbeg : int
            index of x coordinate start
        xend : int
            index of x coordinate end
        ny : int
            y-size of the original array on file
        nx : int
            x-size of the original array on file

        Returns
        -------
        string which needs to be decoded using struct
        """
        if ybeg is None:
            raise ValueError('ERROR: Need to specify YBEG')
        if xbeg is None:
            raise ValueError('ERROR: Need to specify XBEG')
        if yend is None:
            raise ValueError('ERROR: Need to specify YEND')
        if xend is None:
            raise ValueError('ERROR: Need to specify XEND')

        file_content = ''
        for i in xrange(ny):
            if i >= yend:
                break
            if i < ybeg:
                continue
            # position
            #~ print i*nbytes*nx, i*nx
            #~ stop
            pos = i * nbytes * nx + xbeg * nbytes
            f.seek(pos)
            # read content
            bytes_to_read = (xend - xbeg) * nbytes
            r = f.read(bytes_to_read)
            file_content += r
        return file_content

    def get_yearmean(self, mask=None, return_data=False):
        """
        This routine calculate the yearly mean of the data field
        A vector with a mask can be provided for further filtering

        e.g. if all the months from JAN-March are masked as TRUE, the
        result will correspnd to the JFM mean for each year

        Parameters
        ----------
        mask : ndarray (bool)
            temporal mask [time]

        return_data : bool
            specifies if results should be returned as C{Data} object
        """

        if mask is None:
            # if no mask is provided, take everything
            mask = np.ones(len(self.time)).astype('bool')
        else:
            if mask.ndim != 1:
                raise ValueError('Mask needs to be 1-D of length of time!')
            if len(mask) != len(self.time):
                raise ValueError('Mask needs to be 1-D of length of time!')

        #/// get data
        ye = np.asarray(self._get_years())
        years = np.unique(ye)
        dat = self.data

        # calculate mean
        if self.data.ndim == 1:
            res = np.zeros(len(years)) * np.nan
            su = np.zeros(len(years)) * np.nan
        elif self.data.ndim == 3:
            nt, ny, nx = self.data.shape
            res = np.zeros((len(years), ny, nx))
            su = np.zeros((len(years), ny, nx))
        else:
            raise ValueError('Unsupported dimension!')

        for i in xrange(len(years)):
            y = years[i]
            hlp = (ye == y) & mask
            if self.data.ndim == 1:
                res[i] = dat[hlp].mean()
                # calculate sum also (needed for masking in the end)
                su[i] = dat[hlp].sum()
            else:
                res[i, :, :] = dat[hlp, :].mean(axis=0)
                # calculate sum also (needed for masking in the end)
                su[i, :, :] = dat[hlp].sum(axis=0)

        # this is still not the best solution, but works
        res = np.ma.array(res, mask=(su == 0.))

        if return_data:
            # generate data object
            r = self.copy()
            r.data = res
            r.time = self.date2num(
                np.asarray([datetime.datetime(year, 1, 1) for year in years]))
            r.time_cycle = 1
            return r
        else:
            return years, res

    def get_yearsum(self, mask=None, return_data=False):
        """
        This routine calculates the yearly sum of the data field
        A vector with a mask can be provided for further filtering

        e.g. if all the months from JAN-March are masked as TRUE, the
        result will correspnd to the JFM sum for each year

        Parameters
        ----------
        mask : ndarray
            mask [time]
        return_data : bool
            specifies if a Data object shall be returned
        """

        if mask is None:
            # if no mask is provided, take everything
            mask = np.ones(len(self.time)).astype('bool')
        else:
            if mask.ndim != 1:
                raise ValueError('Mask needs to be 1-D of length of time!')
            if len(mask) != len(self.time):
                raise ValueError('Mask needs to be 1-D of length of time!')

        #/// get data
        ye = pl.asarray(self._get_years())
        years = pl.unique(ye)
        dat = self.data

        #/// calculate mean
        res = []
        for y in years:
            hlp = (ye == y) & mask
            if self.data.ndim == 1:
                res.append(dat[hlp].sum())
            else:
                res.append(dat[hlp, :].sum(axis=0))

        res = pl.asarray(res)
        msk = dat.count(0) == 0

        for i in xrange(len(res)):
            res[i, msk] = np.nan
        # mask all data that contained no single valid value!
        res = np.ma.array(res, mask=np.isnan(res))

        if return_data:
            r = self.copy()
            r.data = res
            r.time = self.date2num(
                np.asarray([datetime.datetime(year, 1, 1) for year in years]))
            r.time_cycle = 1
            return r
        else:
            return years, res

#-----------------------------------------------------------------------

    def partial_correlation(self, Y, Z, ZY=None, pthres=1.01, return_object=True):
        """
        perform partial correlation analysis.

        This function calculates the partial correlation between
        variables (self) and Y, removing the effect of variable Z before
        (condition). The partial correlation represents the correlation
        between X and Y, when the common effect, related to Z has been
        removed.

        The function allows to have two datasets used as a condition
        (Z,ZY). Lets say, you have two datasets which were generated
        with a two different forcings which you want to remove from
        X/Y before analyzing their relationship, then this is the
        right choice to specify a second independent variable ZY

        REFERENCES
        ----------
        [1] http://en.wikipedia.org/wiki/Partial_correlation#Using_linear_regression


        Parameters
        ----------

        Y : Data
            variable to calculate with
        Z : Data
            condition for either both variables or if ZY is given,
            then Z is used for SELF only
        pthres : float
            threshold to flag insignificant correlations
        return_object : bool
            specifies if a C{Data} object shall be returned

        Returns
        -------
        r : Data
            returns C{Data} objects with partial correlation parameters
        """

        assert isinstance(Y, Data)
        assert isinstance(Z, Data)

        # if a second condition is given, use it ...
        if ZY is not None:
            assert isinstance(ZY, Data)
        else:
            ZY = Z

        #--- calculate correlations
        rxy, pxy = self.correlate(Y, pthres=pthres)
        rxz, pxz = self.correlate(Z, pthres=pthres)
        rzy, pzy = ZY.correlate(Y, pthres=pthres)

        #--- calculate partial correlation coefficients
        res = (rxy.data - (rxz.data * rzy.data)) / (
            np.sqrt(1. - rxz.data * rxz.data) * np.sqrt(1. - rzy.data * rzy.data))

        if return_object:
            r = self.copy()
            r.time = None
            r.unit = ''
            r.data = res
            r.label = 'partial correlation coefficient'
            return r
        else:
            return res

#-----------------------------------------------------------------------

    def correlate(self, Y, pthres=1.01, spearman=False, detrend=False):
        """
        correlate present data on a grid cell basis with another dataset

        The routine currently supports to calculate either the Pearson
        product-moment correlation coefficient (default) or to calculate
        the Spearman Rank correlation coefficient

        Parameters
        ----------
        Y : Data
            dataset to correlate the present one with. The data set of
            self will be used as X in the calculation
        phres : float
            threshold for masking insignificant pixels. For a
            significance level of 95% pthres needs to be e.g. 0.05
            Then all results with p>pthres will be mased.
        spearman : bool
            option that specifies if spearman correlation should be
            calculated
        detrend : bool
            perform linear detrending before analysis

        Returns
        -------
        R : Data
            correlation coefficient
        P : significance level (p-value)

        Todo
        ----
        * more efficient implementation
        * slope calculation as well?
        * todo: significance correct ??? -- not if stats.mstats.linregress would be used!!!!

        """

        if not Y.data.shape == self.data.shape:
            print Y.data.shape, self.data.shape
            raise ValueError('unequal shapes: correlation not possible!')

        # generate a mask of all samples that are valid in BOTH datasets
        vmask = self.data.mask | Y.data.mask
        vmask = ~vmask

        # copy original data
        xv = self.data.copy()
        yv = Y.data.copy()
        sdim = self.data.shape

        # detrend data if required
        if detrend:
            xv = self.detrend(return_object=True).data.copy()
            yv = Y.detrend(return_object=True).data.copy()

        # ... and reshape it
        nt = len(self.data)
        xv.shape = (nt, -1)
        yv.shape = (nt, -1)
        vmask.shape = (nt, -1)

        # generate new mask for data
        xv.data[xv.mask] = np.nan
        yv.data[yv.mask] = np.nan
        xv[~vmask] = np.nan
        yv[~vmask] = np.nan

        xv = np.ma.array(xv, mask=np.isnan(xv))
        yv = np.ma.array(yv, mask=np.isnan(yv))

        # number of valid data sets where x and y are valid
        nvalid = vmask.sum(axis=0)

        # calculate correlation only for grid cells with at least 3 valid
        # samples
        mskvalid = nvalid > 2
        r = np.ones(sum(mskvalid)) * np.nan
        p = np.ones(sum(mskvalid)) * np.nan

        xn = xv[:, mskvalid]
        yn = yv[:, mskvalid]  # copy data
        nv = nvalid[mskvalid]

        # do correlation calculation; currently using np.ma.corrcoef as this
        # supports masked arrays, while stats.linregress doesn't!

        if mskvalid.sum() > 0:
            if spearman:
                res = [stats.mstats.spearmanr(xn[:, i], yn[:, i])
                       for i in xrange(sum(mskvalid))]
                res = np.asarray(res)
                r = res[:, 0]
                p = res[:, 1]
                r[p > pthres] = np.nan

            else:  # Pearson product-moment correlation
                res = [np.ma.corrcoef(xn[:, i], yn[:, i]) for i in xrange(sum(
                    mskvalid))]  # <<<< as an alternative one could use stats.mstats.linregress ; results are however equal for R-VALUE, but NOT for P-value, here mstats.linregress seems to be buggy!, see unittests
                res = np.asarray(res)
                r = res[:, 0, 1]  # correlation coefficient
                p = get_significance(r, nv)
                r[p > pthres] = np.nan
        else:
            r = None
            p = None

        # remap to original geometry
        R = np.ones(xv.shape[1]) * np.nan  # matrix for results
        P = np.ones(xv.shape[1]) * np.nan  # matrix for results
        R[mskvalid] = r
        P[mskvalid] = p
        orgshape = (sdim[1], sdim[2])
        R = R.reshape(orgshape)  # generate a map again
        P = P.reshape(orgshape)

        R = np.ma.array(R, mask=np.isnan(R))
        P = np.ma.array(P, mask=np.isnan(P))

        RO = self.copy()
        RO.data = R
        if spearman:
            RO.label = '$r_{spear}$: ' + self.label + ' vs. ' + Y.label
        else:
            RO.label = '$r_{pear}$: ' + self.label + ' vs. ' + Y.label
        RO.unit = ''

        PO = self.copy()
        PO.data = P
        PO.label = 'p-value'  # + self.label + ' ' + y.label
        PO.unit = ''

        return RO, PO

    def get_temporal_mask(self, v, mtype='monthly'):
        """
        returns a temporal mask which marks specific months or years
        that match the desired mask

        v : list
            list of values to be analyzed, e.g. [1,2,3] for JAN/FEB/MAR

        mtype : str
            specifies which mask should be applied
            valid values: ['monthly','yearly']

        Example
        -------
        >>> self.get_temporal_mask([1,2,3], mtype='monthly')
        will return a mask, where the months of Jan-Mar are set to True
        this can be used e.g. further with the routine get_yearmean()
        """

        valid_types = ['monthly', 'yearly']
        if mtype in valid_types:
            pass
        else:
            raise ValueError('Invalid type for mask generation %s' % mtype)

        # get months or years
        if mtype == 'monthly':
            vals = pl.asarray(self._get_months())
        elif mtype == 'yearly':
            vals = pl.asarray(self._get_years())
        else:
            raise ValueError('Invalid type for mask generation %s ' % mtype)

        # generate mask with all months
        mask = pl.zeros(self.nt).astype('bool')

        for m in v:
            hlp = vals == m
            mask[hlp] = True
        return pl.asarray(mask)

    def get_climatology(self, return_object=False, nmin=1, ensure_start_first=True):
        """
        calculate climatological mean for a time increment
        specified by self.time_cycle

        *Note*: one can not assume that the climatology starts from
        January if you use a time_cycle = 12
        Instead, the climatology simply starts with the value which
        corresponds to the first value of the data.

        Parameters
        ----------
        return_object : bool
            specifies if a Data object shall be returned
        nmin : int
            specifies the minimum number of datasets used for
            climatology; else the result is masked
        ensure_start_first : bool
            ensure that the timeseries of the resulting climatology
            always starts with the first date. If you have e.g.
            a datasets that starts with dates in March, then also the
            climatology will start in March. Using this option will
            ensure that the climatology will then start in January
        """
        if hasattr(self, 'time_cycle'):
            pass
        else:
            raise ValueError(
                'Climatology can not be calculated without a valid time_cycle')

        # generate output fields
        if self.data.ndim > 1:
            clim = np.ones(np.shape(self.data[0:self.time_cycle, :])) * np.nan
            slim = np.ones(np.shape(self.data[0:self.time_cycle, :])) * np.nan
        else:
            clim = np.ones(np.shape(self.data[0:self.time_cycle])) * np.nan
            slim = np.ones(np.shape(self.data[0:self.time_cycle])) * np.nan

        if clim.ndim == 1:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle] = self.data[
                    i::self.time_cycle].mean(axis=0)
                slim[i::self.time_cycle] = self.data[
                    i::self.time_cycle].sum(axis=0)
        elif clim.ndim == 2:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle, :] = self.data[
                    i::self.time_cycle, :].mean(axis=0)
                slim[i::self.time_cycle, :] = self.data[
                    i::self.time_cycle, :].sum(axis=0)
        elif clim.ndim == 3:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle, :, :] = self.data[
                    i::self.time_cycle, :, :].mean(axis=0)
                slim[i::self.time_cycle, :, :] = self.data[
                    i::self.time_cycle, :, :].sum(axis=0)
        else:
            raise ValueError('Invalid dimension when calculating climatology')

        n = slim / clim
        del slim  # number of data taken into account for climatology
        clim = np.ma.array(
            clim, mask=(np.isnan(clim) | (n < nmin) | np.isnan(n)))
        del n

        # create a data object
        r = self.copy()
        r.label += ' - climatology'
        r.data = clim
        r.time = []
        for i in xrange(self.time_cycle):
            r.time.append(self.time[i])  # use data for the first timesteps
        r.time = np.asarray(r.time)
        r.adjust_time(year=1200)  # set some arbitrary time

        if len(r.time) != len(r.data):
            print len(r.time)
            print len(r.data)
            raise ValueError(
                'Data and time are inconsistent in get_climatology()')

        if ensure_start_first:
            r._shift_time_start_firstdate()

        if return_object:
            return r
        else:
            return r.data

    def _shift_time_start_firstdate(self):
        """
        shift dataset that the timeseries is ensured to be in ascending order
        usefull e.g. if you have a climatology and this does not start with January
        and you want to shift it automatically
        """

        # search for point in timeseries where break occurs
        di = np.diff(self.time)
        m = di < 0.
        if m.sum() == 0:
            n = 0
        elif m.sum() == 1:  # a single breakpoint
            n = m.argmax() + 1  # position where the break in timeseries occurs
        else:
            raise ValueError(
                'More than a single breakpoint found. Can not process this data as it is not in cyclic ascending order')

        # shift data now
        self.timeshift(n, shift_time=True)

    def get_deseasonalized_anomaly(self, base=None, ensure_start_first=True):
        """
        calculate deseasonalized anomalies

        The anomalies are calculated by removing the mean seasonal cycle
        The seasonality is specified by the time_cycle variable of the
        Data object

        Parameters
        ----------
        base : str
            specifies the base to be used for the climatology
            'all': use the WHOLE original dataset as a reference
            'current' : use current data as a reference

        Returns
        -------
        return a Data object of deseasonalized anomalies
        """

        if base == 'current':
            clim = self.get_climatology(ensure_start_first=ensure_start_first)
        elif base == 'all':
            # if raw climatology not available so far, try to calculate it
            if hasattr(self, '_climatology_raw'):
                clim = self._climatology_raw
            else:
                if hasattr(self, 'time_cycle'):
                    self._climatology_raw = self.get_climatology(
                        ensure_start_first=ensure_start_first)
                    clim = self._climatology_raw
                else:
                    raise ValueError(
                        'Climatology can not be calculated because of missing time_cycle!')
        else:
            raise ValueError('Anomalies can not be calculated, invalid BASE')

        if hasattr(self, 'time_cycle'):
            pass
        else:
            raise ValueError(
                'Anomalies can not be calculated without a valid time_cycle')
        ret = np.ones_like(self.data) * np.nan

        if ret.ndim == 1:
            for i in xrange(self.time_cycle):
                ret[i::self.time_cycle] = self.data[i::self.time_cycle] - \
                    clim[i]
        elif ret.ndim == 2:
            for i in xrange(self.time_cycle):
                ret[i::self.time_cycle, :] = self.data[
                    i::self.time_cycle, :] - clim[i, :]
        elif ret.ndim == 3:
            for i in xrange(self.time_cycle):
                ret[i::self.time_cycle, :, :] = self.data[i::self.time_cycle, :, :] - clim[i, :, :]
        else:
            raise ValueError('Invalid dimension when calculating anomalies')
        ret = np.ma.array(ret, mask=(np.isnan(ret) | self.data.mask))

        # return a data object
        res = self.copy()
        res.data = ret.copy()
        res.label = self.label + ' anomaly'

        return res

    def condstat(self, M):
        """
        Conditional statistics of data

        This routine calculates conditions statistics over the current data. Given a mask M, the routine calculates for
        each unique value in M the mean, stdv, min and max from the current data

        Parameters
        ----------
        M : ndarray or Data
            mask to be used. Needs to be a 2D array of dimension ny x nx

        Returns
        -------
        res : dict
            dictionary with results where each entry has shape (nt,nvals) with nvals beeing the number of
            unique ID values in the mask
            res = {'id': vals, 'mean': means, 'sum': sums, 'min': mins, 'max': maxs, 'std': stds}

        Example
        -------
        > Let us assume you have a data object D and we assign some sample data to it and generate a mask with a few pixels
        > D.data = pl.randn(100,3,1) #some sample data
        > msk = np.asarray([[1,1,3],]).T #(3,1) mask
        > res = D.condstat(msk)  # calculate conditional statistics
        > This returns a dictionary with the following keys ['max', 'sum', 'min', 'id', 'mean']

        """

        if isinstance(M, Data):
            m = M.data
        else:
            m = M

        if self.data.ndim == 2:
            if self.data.shape != m.shape:
                print self.shape
                print m.shape
                raise ValueError('Invalid geometry!')
        elif self.data.ndim == 3:
            if self.data[0, :, :].shape != m.shape:
                print self.shape
                print m.shape
                raise ValueError('Invalid geometry!')
        else:
            raise ValueError('Unsupported Data geometry!')

        # calculate conditional statistics
        vals = np.unique(m).astype(int)
        if isinstance(vals, np.ma.core.MaskedArray):
            # for masked arrays the unique() returns also a placeholder
            # for the invalid data. To avoid that, we ensure here
            # that only the real mask values are stored
            vals = vals.data[~vals.mask]

        def _get_stat(a, msk, v):
            # get statistics of a single 2D field and a specific value v
            # a: masked array
            # msk: mask to use for analysis
            # v: float

            x = a[msk == v].flatten()
            m = np.nan
            s = np.nan
            mi = np.nan
            ma = np.nan
            su = np.nan

            if len(x) > 0:
                m = x.mean()
                mi = x.min()
                ma = x.max()
                su = x.sum()
            if len(x) > 2:
                s = x.std()

            return m, s, su, mi, ma

        if self.data.ndim == 2:
            means = np.ones((1, len(vals))) * np.nan
            sums = np.ones((1, len(vals))) * np.nan
            stds = np.ones((1, len(vals))) * np.nan
            mins = np.ones((1, len(vals))) * np.nan
            maxs = np.ones((1, len(vals))) * np.nan

            for i in xrange(len(vals)):
                means[0, i], stds[0, i], sums[0, i], mins[0,
                                                          i], maxs[0, i] = _get_stat(self.data, m, vals[i])

        elif self.data.ndim == 3:
            nt = len(self.data)
            means = np.ones((nt, len(vals))) * np.nan
            sums = np.ones((nt, len(vals))) * np.nan
            stds = np.ones((nt, len(vals))) * np.nan
            mins = np.ones((nt, len(vals))) * np.nan
            maxs = np.ones((nt, len(vals))) * np.nan

            # calculate for each timestep and value the conditional statistic
            for t in xrange(nt):
                for i in xrange(len(vals)):
                    means[t, i], stds[t, i], sums[t, i], mins[t, i], maxs[t, i] = _get_stat(
                        self.data[t, :, :],
                        m, vals[i])
        else:
            raise ValueError('Invalid geometry!')

        # output arrays are all of shape (nt,nvals)
        # now we reformat output as such that the ID is the key for the
        # dictionary
        res = {}
        try:
            thedate = self.date
        except:
            thedate = None

        for i in xrange(len(vals)):
            id = vals[i]
            res.update(
                {id: {'mean': means[:, i], 'std': stds[:, i], 'sum': sums[:, i], 'min': mins[:, i],
                      'max': maxs[:, i], 'time': thedate}})
        if len(res) == 0:
            return None
        else:
            return res

    def set_time(self):
        """
        convert times that are in a specific format
        If the time string is already known to be handled by the
        netcdf4time module, then nothing is happening
        """
        if self.time_str is None:
            raise ValueError(
                'ERROR: time can not be determined, as units for time not available!')
        if not hasattr(self, 'calendar'):
            raise ValueError('ERROR: no calendar specified!')
        if not hasattr(self, 'time'):
            raise ValueError('ERROR: no time specified!')

        if self.time_str == 'day as %Y%m%d.%f':
            # in case of YYYYMMDD, convert to other time with basedate
            # 0001-01-01 00:00:00
            self._convert_time()
        if self.time_str == 'day as YYYYMMDD':
            self._convert_time_YYYYMMDD()
        elif self.time_str == 'month as %Y%m.%f':
            self._convert_timeYYYYMM()
        elif self.time_str == 'year as %Y.%f':
            self._convert_timeYYYY()
        elif 'months since' in self.time_str:
            # months since is not supported by netCDF4 library at the moment.
            # Therefore implementation here.
            self._convert_monthly_timeseries()

            #--- time conversion using netCDF4 library routine ---
            # actually nothing needs to be done, as everything shall
            # be handled by self.num2date() in all subsequent subroutines
            # to properly handle difference in different calendars.

    def _get_date_from_month(self, nmonths):
        """
        calculate a datetime object for a time given in 'months since'
        a basedate. The routine increments itteratively the number of
        months and returns a datetime object.

        This is done for a *single* timestep!

        Parameters
        ----------
        nmonths : int, float
            time as numeric value (number of months since basedate)

        Returns
        -------
        d : datetime
            datetime object with actual date
        """

        if not 'months since' in self.time_str:
            print(self.time_str)
            raise ValueError(
                'This routine is only for conversion of monthly data!')

        basedate = self.time_str.split('since')[1].lstrip()

        # start date
        start_date = pl.datestr2num(basedate)
        act_date = start_date * 1.

        for i in xrange(int(nmonths)):  # increment months
            d = pl.num2date(act_date)
            # number of days in current month
            ndays = monthrange(d.year, d.month)[1]
            act_date += ndays

        return pl.num2date(act_date)

    def _convert_monthly_timeseries(self):
        """
        convert monthly timeseries to a daily timeseries
        """
        if self.calendar not in ['standard', 'gregorian', None]:
            print self.calendar
            raise ValueError('Not sure if monthly timeseries conversion \
                                works with this calendar!')

        newtime = [self._get_date_from_month(t)
                   for t in self.time]  # ... estimate new time
        self.calendar = 'standard'
        self.time_str = 'days since 0001-01-01 00:00:00'
        # plus one because of the num2date() basedate definition
        self.time = pl.date2num(newtime) + 1.

    def apply_temporal_subsetting(self, start_date, stop_date):
        """
        perform temporal subsetting of data

        Parameters
        ----------
        start_date : datetime
            start time
        stop_date : datetime
            stop time

        Returns
        -------
        Returns nothing, but modifies the actual data
        """
        i1, i2 = self._get_time_indices(start_date, stop_date)
        self._temporal_subsetting(i1, i2)

    def _temporal_subsetting(self, i1, i2):
        """
        perform temporal subsetting of the data based on given
        time indices i1,i2

        Parameters
        ----------
        i1 : int
            start time index
        i2 : int
            stop time index

        Returns
        -------
        Returns nothing, but modifies the current Data object
        """

        if i2 < i1:
            raise ValueError('Invalid indices _temporal_subsetting')
        # incremet last index, as otherwise the last dataset is missing!
        i2 += 1
        if i2 > len(self.time):
            i2 = len(self.time)
        self.time = self.time[i1:i2]

        if self.data.ndim == 3:
            self.data = self.data[i1:i2, :, :]
        elif self.data.ndim == 2:
            # data has already been squeezed and result was 2D (thus without
            # time!)
            if self.squeezed:
                print(
                    'Data was already squeezed: no temporal subsetting is performed!')
            else:
                self.data = self.data[i1:i2, :]
        elif self.data.ndim == 1:  # single temporal vector assumed
            self.data = self.data[i1:i2]
        else:
            raise ValueError('Error temporal subsetting: invalid dimension!')

    def align(self, y, base=None):
        """
        Temporal alignment of two Data objects.
        The datasets need to have a similar time stepping. The desired timestepping is explicitely specified
        by the user in the *base* argument. It is obligatory to provide this argument.write

        Parameters
        ----------
        y : Data
            Data object that should be aligned with current data
        base : str
            specifies the temporal basis for the alignment. Data needs to have been preprocessed already with such
            a time stepping. Currently supported values: ['month','day']

        Returns
        -------
        x, y : Data
            returns two dataobjects which are aligned to eacht other
        """

        assert (isinstance(y, Data))
        if base is None:
            raise ValueError(
                'You need to specify the base for the alignment [month,day] !')

        # ensure ascending time order
        x = self.copy()

        if not x._is_sorted():
            raise ValueError('Time series in X is not sorted ascending!')
        if not y._is_sorted():
            raise ValueError('Time series in Y is not sorted ascending!')

        if base == 'month':
            if not x._is_monthly():
                print x.date
                raise ValueError('Dataset X is not monthly data!')
            if not y._is_monthly():
                print y.date
                raise ValueError('Dataset Y is not monthly data!')
        elif base == 'day':
            if not x._is_daily():
                print x.date
                raise ValueError('Dataset X is not daily data!')
            if not y._is_daily():
                print y.date
                raise ValueError('Dataset Y is not daily data!')
        else:
            raise ValueError('Unsupported base for alignment!')

        # min/max dates
        ymin = y._get_mindate(base=base)
        ymax = y._get_maxdate(base=base)
        xmin = x._get_mindate(base=base)
        xmax = x._get_maxdate(base=base)

        start = None
        stop = None

        # search first for the dataset which starts first
        if xmin <= ymin:
            xfirst = True
        else:
            xfirst = False

        # check if overlap at all
        err = 0
        if xfirst:
            if xmax < ymin:
                err += 1  # no overlap
        else:
            if ymax < xmin:
                err += 1
        if err > 0:
            return None, None

        # from here onwards we know that there is an overlap
        start = max(xmin, ymin)
        stop = min(xmax, ymax)

        x1, x2 = x._get_time_indices(start, stop)
        y1, y2 = y._get_time_indices(start, stop)

        x._temporal_subsetting(x1, x2)
        y._temporal_subsetting(y1, y2)

        return x, y

    def interp_time(self, d, method='linear'):
        """
        interpolate data matrix in time. The existing data is
        interpolated to a new temporal spaceing that is specified by
        the time vector argument 't'. The interpolation is done,
        by constructing a weighting matrix which basically performs
        a linear interpolation as y = w*x(1) + (1-w)*x(2)


        Parameters
        ----------
        d : ndarray of datetime objects
            vector of time where the data should be interpolated to.

        method : str
            option to specify interpolation method. At the moment,
            only linear interpolation is supported! ['linear']

        Returns
        -------
        returns a new C{Data} object that contains the interpolated values

        # TODO : still some boundary effects for last timestep
        """

        # check if timezone information available. If not, then
        # set to UTC as default
        d = np.asarray([datetime.datetime(x.year, x.month, x.day, x.hour, x.minute, x.second, 0, pytz.UTC)
                       for x in d])

        if method not in ['linear']:
            raise ValueError(
                'Only linear interpolation supported at the moment so far!')

        # checks
        if self.data.ndim != 3:
            raise ValueError(
                'Interpolation currently only supported for 3D arrays!')
        if not np.all(np.diff(self.date2num(d)) > 0):
            raise ValueError(
                'Input time array is not in ascending order! This must not happen! Please ensure ascending order')
        if not np.all(np.diff(self.time) > 0):
            raise ValueError(
                'Time array of data is not in ascending order! This must not happen! Please ensure ascending order')

        nt0, ny, nx = self.shape  # original dimensions
        nt = len(d)  # target length of reference time

        # copy data
        X = self.data.copy()
        X.shape = (nt0, -1)  # [time,npix]
        nt0, npix = X.shape

        # preliminary checks
        f_err = False

        # A) all data is BEFORE desired period
        if self.date.max() < d.min():
            print self.date.max(), d.min()
            print(
                'WARNING: specified time period is BEFORE any data availability. NO INTERPOLATION CAN BE DONE!')
            f_err = True

        # B) all data is AFTER desired period
        if self.date.min() > d.max():
            print self.date.min(), d.max()
            print(
                'WARNING: specified time period is AFTER any data availability. NO INTERPOLATION CAN BE DONE!')
            f_err = True

        if f_err:
            tmp = np.zeros((nt, ny, nx))
            r = self.copy()
            r.data = np.ma.array(tmp, mask=tmp > 0.)
            r.time = self.date2num(d)  # return some dummy result
            return r

        # construct weighting matrix
        W = np.zeros((nt, nt0))
        i1 = 0
        i2 = 1  # indices in original data
        f_init = True
        for i in xrange(nt - 1):
            # 1) find start of interpolation period
            if f_init:
                # do nothing while data coverage not reached yet
                while self.date[i2] <= d[0]:
                    i2 += 1
                    continue
            f_init = False
            i1 = i2 - 1
            if i1 < 0:
                raise ValueError('Invalid index i1')

            # increment
            if i2 < nt0:
                if self.date[i2] < d[i]:
                    if i2 <= nt0 - 1:
                        i2 += 1
            if i1 < nt0 - 1:
                if self.date[i1 + 1] < d[i]:
                    i1 += 1
            else:
                continue

            # 2) check consistency
            if self.date[i1] > d[i]:
                # the first timeperiod with valid data has not been reached yet
                # ... loop
                continue
            if i2 >= nt0:
                continue

            if self.date[i2] < d[i]:
                print self.date[i1], d[i], self.date[i2]
                print i1, i, i2
                raise ValueError('interp_time: this should not happen!')

            if i2 > nt0 - 1:
                break

            #... here we have valid data
            t1 = self.date2num(self.date[i1])
            t2 = self.date2num(self.date[i2])

            W[i, i1] = (t2 - self.date2num(d[i])) / (t2 - t1)
            W[i, i2] = 1. - W[i, i1]

            #... now increment if needed
            if i < (nt0 - 1):
                if t2 < self.date2num(d[i + 1]):
                    i1 += 1
                    i2 += 1

        #/// generate interpolation Matrix and perform interpolation
        # could become a problem for really large matrices!
        N = np.ma.dot(W, X)
        # avoid boundary problem (todo: where is the problem coming from ??)
        N[nt - 1, :] = np.nan
        # mask all data that is outside of valid time period
        msk = (d < self.date.min()) | (d > self.date.max())
        N[msk, :] = np.nan
        N.shape = (nt, ny, nx)

        res = self.copy()

        res.time = self.date2num(d)
        res.time_str = self.time_str
        res.calendar = self.calendar
        res.data = np.ma.array(N, mask=np.isnan(N))
        del N

        return res

    def _get_time_indices(self, start, stop):
        """
        determine time indices start/stop based on data timestamps
        and desired start/stop dates

        Parameters
        ----------
        start : datetime
            start time
        stop : datetime
            stop time

        Returns
        -------
        returns start/stop indices (int)
        """

        def _check_timezone(d):
            # if no timezone, then set it
            if d.tzinfo is None:
                t = datetime.datetime(
                    d.year, d.month, d.day, d.hour, d.minute, d.second, d.microsecond, pytz.UTC)
                return t
            else:
                return d

        # no subsetting
        if start is None and stop is None:
            return 0, len(self.time) - 1
        if start is None:
            start = self.date[0]
        else:
            assert (isinstance(start, datetime.datetime))
            start = _check_timezone(start)
        if stop is None:
            stop = self.date[-1]
        else:
            assert (isinstance(stop, datetime.datetime))
            stop = _check_timezone(stop)
        if stop < start:
            raise ValueError('Error: startdate > stopdate')

        s1 = self.date2num(start)
        s2 = self.date2num(stop)

        #- check that time is increasing only
        if any(np.diff(self.time)) < 0.:
            raise ValueError('Error _get_time_indices: Time is not increasing')

        # determine indices
        m1 = abs(self.time - s1).argmin()
        m2 = abs(self.time - s2).argmin()

        if self.time[m1] < s1:
            m1 += 1
        if self.time[m2] > s2:
            m2 -= 1
        if m2 < m1:
            raise ValueError('Something went wrong _get_time_indices')
        return m1, m2

    def _get_years(self):
        """
        get years from timestamp
        """
        return [x.year for x in self.date]

    def _get_months(self):
        """
        get months from timestamp
        """
        return [x.month for x in self.date]

    def _get_days_per_month(self):
        """ get number of days for each month """
        return np.asarray([calendar.monthrange(x.year, x.month)[1] for x in self.date])

    def _mesh_lat_lon(self):
        """
        In case that the Data object lat/lon is given in vectors
        the coordinates are mapped to a 2D field. This routine
        resets the Data object coordinates
        """
        if (plt.isvector(self.lat)) & (plt.isvector(self.lon)):
            LON, LAT = np.meshgrid(self.lon, self.lat)
            self.lon = LON
            self.lat = LAT
        else:
            pass

    def read_netcdf(self, varname, netcdf_backend='netCDF4'):
        """
        read data from netCDF file

        varname : str
            name of variable to be read
        """
        File = NetCDFHandler(netcdf_backend=netcdf_backend)
        File.open_file(self.filename, 'r')

        if self.verbose:
            print 'Reading file ', self.filename
        if not varname in File.get_variable_keys():
            if self.verbose:
                self._log_warning(
                    'WARNING: data can not be read. Variable not existing! ', varname)
            File.close()
            return None

        try:
            data = File.get_variable(varname)
        except:
            print('ERROR when reading variable %s' % varname)
            return None
        var = File.get_variable_handler(varname)

        if data.ndim > 3:
            if self.level is None:
                print data.shape
                raise ValueError(
                    '4-dimensional variables not supported yet! Either remove a dimension or specify a level!')
            else:
                # [time,level,ny,nx ] --> [time,ny,nx]
                data = data[:, self.level, :, :]

        # in case of vector data, generate a dummy dimension
        if data.ndim == 1:
            tmp = np.zeros((1, len(data)))
            tmp[:] = data[:] * 1.
            data = tmp * 1.
            del tmp

        self.fill_value = None
        self.fill_value = File._get_fill_value(varname)
        if self.fill_value is not None:
            msk = data == self.fill_value
            # set to nan, as otherwise problems with masked and scaled data
            data[msk] = np.nan
            data = np.ma.array(data, mask=np.zeros(data.shape).astype(
                'bool'))  # generate an empty mask first to ensure that the mask has the same geometry as the data!
            data.mask[np.isnan(data)] = True
        else:
            data = np.ma.array(data, mask=np.zeros(data.shape).astype('bool'))
            self.fill_value = -99999.

        #--- scale factor
        scal = File._get_scale_factor(varname)
        self._scale_factor_netcdf = scal * 1.

        offset = File._get_add_offset(varname)
        self._add_offset_netcdf = offset * 1.

        #data = data * scal + offset
        data *= scal
        data += offset

        # set longname of variable
        self.long_name = File._get_long_name(varname)

        # check if file has cell_area attribute and only use it if it has not
        # been set by the user
        if 'cell_area' in File.get_variable_keys() and self.cell_area is None:
            self.cell_area = File.get_variable('cell_area')

        # set units if possible; if given by user, this is taken
        # otherwise unit information from file is used if available
        if self.unit is None:
            self.unit = File._get_unit(varname)

        if self.time_var in File.get_variable_keys():
            tvar = File.get_variable_handler(self.time_var)
            if hasattr(tvar, 'units'):
                self.time_str = tvar.units
            else:
                self.time_str = None

            if hasattr(tvar, 'calendar'):
                self.calendar = tvar.calendar
                # when climatology means, reset calendar to standard
                if self.calendar == 'climatology_bounds':
                    self.calendar = 'standard'
            else:
                print 'WARNING: no calendar specified!'
                self.calendar = 'standard'
        else:
            self.time = None
            self.time_str = None
        File.close()

        return data

    def temporal_trend(self, return_object=False, pthres=1.01):
        """
        calculate temporal trend of the data over time
        the slope of the temporal trend has unit [dataunit/day]

        Parameters
        ----------

        return_object : bool
            specifies if a C{Data} object shall be returned [True]
            or if a numpy array shall be returned [False]
        pthres : float
            specifies significance threshold; all values above this threshold will be masked

        Returns
        -------
        The following variables are returned:
        correlation, slope, intercept, p-value
        """
        # time difference in days
        dt = np.asarray([x.days for x in self.date - self.date[0]]
                        ).astype('float')
        # ensure that a masked array is used
        x = np.ma.array(dt, mask=dt != dt)

        R, S, I, P, C = self.corr_single(x, pthres=pthres)

        R.label = self.label + '(correlation)'
        S.label = self.label + '($\partial x / \partial t$)'
        I.label = self.label + '(offset)'
        P.label = self.label + '(p-value)'

        if return_object:
            S.unit += ' / day'
            R.unit = '-'
            I.unit = self.unit
            P.unit = '-'
            return R, S, I, P
        else:
            return R.data, S.data, I.data, P.data

    def timmean(self, return_object=False):
        """
        calculate temporal mean of data field

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True]; else a numpy array is returned
        """
        if self.data.ndim == 3:
            res = self.data.mean(axis=0)
        elif self.data.ndim == 2:
            res = self.data.copy()  # no temporal averaging
        else:
            print self.data.ndim
            raise ValueError(
                'Temporal mean can not be calculated as dimensions do not match!')

        if return_object:
            tmp = self.copy()
            tmp.data = res
            if hasattr(tmp, 'time'):
                del tmp.time
            return tmp
        else:
            return res

    def timmin(self, return_object=False):
        """
        calculate temporal minimum of data field

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True]; else a numpy array is returned
        """

        if self.data.ndim == 3:
            res = self.data.min(axis=0)
        elif self.data.ndim == 2:
            # no temporal averaging
            res = self.data.copy()
        else:
            print self.data.ndim
            raise ValueError(
                'Temporal minimum can not be calculated as dimensions do not match!')

        if return_object:
            tmp = self.copy()
            tmp.data = res
            if hasattr(tmp, 'time'):
                del tmp.time
            return tmp
        else:
            return res

    def timmax(self, return_object=False):
        """
        calculate temporal maximum of data field

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True]; else a numpy array is returned
        """
        if self.data.ndim == 3:
            res = self.data.max(axis=0)
        elif self.data.ndim == 2:
            # no temporal averaging
            res = self.data.copy()
        else:
            print self.data.ndim
            raise ValueError(
                'Temporal maximum can not be calculated as dimensions do not match!')

        if return_object:
            tmp = self.copy()
            tmp.data = res
            if hasattr(tmp, 'time'):
                del tmp.time
            return tmp
        else:
            return res

    def timcv(self, return_object=True):
        """
        calculate temporal coefficient of variation

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True]; else a numpy array is returned

        Test
        ----
        unittest implemented
        """
        res = self.timstd(return_object=False) / \
            self.timmean(return_object=False)
        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy()
                tmp.data = res
                tmp.label = self.label + ' (CV)'
                tmp.unit = '-'
                return tmp
        else:
            return res

#-----------------------------------------------------------------------

    def normalize(self, return_object=True):
        """
        normalize data by removing the mean and dividing by the standard deviation
        normalization is done for each grid cell

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned
        """

        if self.data.ndim != 3:
            raise ValueError('Normalization only possible for 3D data!')

        if return_object:
            d = self.copy()
        else:
            d = self

        d.sub(d.timmean(return_object=True), copy=False)
        d.div(d.timstd(return_object=True), copy=False)

        if return_object:
            return d
        else:
            return None

    def timstd(self, return_object=False):
        """
        calculate temporal standard deviation of data field

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True];
            else a numpy array is returned
        """
        if self.data.ndim == 3:
            res = self.data.std(axis=0)
        elif self.data.ndim == 2:
            # no temporal averaging
            res = None
        else:
            raise ValueError(
                'Temporal standard deviation can not be calculated as dimensions do not match!')

        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy()
                tmp.data = res
            if hasattr(tmp, 'time'):
                del tmp.time
                return tmp
        else:
            return res

#-----------------------------------------------------------------------

    def timvar(self, return_object=False):
        """
        calculate temporal variance of data field

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True];
            else a numpy array is returned
        """
        if self.data.ndim == 3:
            res = self.data.var(axis=0)
        elif self.data.ndim == 2:
            # no temporal averaging
            res = None

        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy()
                tmp.data = res
                if hasattr(tmp, 'time'):
                    del tmp.time
                return tmp
        else:
            return res

    def timsum(self, return_object=False):
        """
        calculate temporal sum of data field

        Parameters
        ----------
        return_object : bool
            specifies if a C{Data} object shall be returned [True];
            else a numpy array is returned
        """
        if self.data.ndim == 3:
            pass
        elif self.data.ndim == 1:
            pass
        else:
            raise ValueError(
                'Temporal sum can not be calculated as dimensions do not match!')
        res = self.data.sum(axis=0)
        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy()
                tmp.data = res
                if hasattr(tmp, 'time'):
                    del tmp.time
                return tmp
        else:
            return res

#-----------------------------------------------------------------------

    def timn(self, return_object=False):
        """
        calculate number of valid samples per time
        The implementation is by using timmean() and timsum()

        Parameters
        ----------
        return_object : bool
            return Data object
        """
        res = self.timsum() / self.timmean()

        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy()
                tmp.data = res
                return tmp
        else:
            return res

    def hp_filter(self, lam, return_object=True):
        """
        implements the Hodrick-Prescott filter

        Todo
        ----
        - use more efficient implementation from statsmodels
        - support HP filter for multidimensional data
        - implement unittests
        - how to handle gaps ???

        Parameters
        ----------
        lam : float
            lambda parameter of HP filter. The larger it is, the smoother
            the resulting trend timeseries will be
        return_object : bool
            return a Data object if True

        References
        ----------
        http://python4econ.blogspot.de/2012/05/hodrick-prescott-filter.html
        """
        from scipy import linalg as la
        from scipy import sparse
        import scipy as sp

        def _hp_filter(y, w):
            # make sure the inputs are the right shape
            m, n = y.shape
            if m < n:
                y = y.T
                m = n
            a = sp.array([w, -4 * w, ((6 * w + 1) / 2.)])
            d = sp.tile(a, (m, 1))

            d[0, 1] = -2. * w
            d[m - 2, 1] = -2. * w
            d[0, 2] = (1 + w) / 2.
            d[m - 1, 2] = (1 + w) / 2.
            d[1, 2] = (5 * w + 1) / 2.
            d[m - 2, 2] = (5 * w + 1) / 2.

            B = sparse.spdiags(d.T, [-2, -1, 0], m, m)
            B = B + B.T
            # report the filtered series, s
            return sp.dot(la.inv(B.todense()), y)

        if self.ndim != 1:
            if self.ndim == 3:
                if (self.shape[1] == 1) and (self.shape[2] == 1):
                    pass
                else:
                    print self.shape
                    raise ValueError(
                        'HP filter currently only implemented for 1D data! (A)')
            else:
                print self.shape
                raise ValueError(
                    'HP filter currently only implemented for 1D data! (B)')

        if lam < 0.:
            raise ValueError('HP filter needs lambda>0. as input!')

        # the HP filter is based on the log() of the data
        # avoid therefore negative numbers
        dmin = self.data.min()
        x = self.data.flatten() - dmin + 1.
        # work only on valid data; note that this approech is probably not the
        # best one!
        msk = x.mask
        hp = _hp_filter(np.log(np.asarray([x[~msk]])), lam)  # 2D input needed
        y = np.ones_like(x) * np.nan
        y[~msk] = np.exp(hp)
        y += dmin - 1.
        y[msk] = np.nan

        if return_object:
            r = self.copy()
            tmp = np.ones((self.nt, 1, 1)) * np.nan
            tmp[:, 0, 0] = y[:]
            r.data = np.ma.array(tmp, mask=np.isnan(tmp))
            return r
        else:
            return y

    def _get_weighting_matrix(self):
        """
        get matrix for area weighting of grid cells. For each timestep
        the weights are calculated as a function of either the  number
        of valid grid cells or all grid cells. Which one of the two
        approaches is used depends on self.weighting_type

        The returned array contains weights for each timestep. The sum
        of these weights is equal to one for each timestep.

        Returns
        -------
        w : ndarray
            weighting matrix in same geometry as original data
        """

        normtype = self.weighting_type

        if normtype in ['valid', 'all']:
            pass
        else:
            raise ValueError('Invalid option for normtype: %s' % normtype)

        w = np.zeros(self.data.shape)

        if self.data.ndim == 2:
            if normtype == 'valid':
                m = ~self.data.mask
                self.totalarea = self.cell_area[m].sum()
                w[m] = self.cell_area[m] / self.totalarea
                w = np.ma.array(w, mask=~m)
            else:
                self.totalarea = self.cell_area.sum()
                w = self.cell_area / self.totalarea
                w = np.ma.array(w, mask=w != w)
            return w

        elif self.data.ndim == 3:
            nt = len(self.data)

            # 1) repeat cell area nt-times
            cell_area = self.cell_area.copy()
            s = np.shape(cell_area)
            cell_area.shape = (-1)
            if len(s) == 2:
                w = cell_area.repeat(nt).reshape((s[0] * s[1], nt)).T
            elif len(s) == 1:
                w = cell_area.repeat(nt).reshape((1, nt)).T
            else:
                print 'nt: ', nt
                print 's: ', s
                print 'len(s): ', len(s)
                raise ValueError('Invalid geometry!')

            w.shape = self.data.shape  # geometry is the same now as data

            if normtype == 'valid':
                # 2) mask areas that do not contain valid data
                w = np.ma.array(w, mask=self.data.mask)

                # 3) calculate for each time the sum of all VALID grid cells
                # --> normalization factor
                no = w.reshape(nt, -1).sum(axis=1)  # ... has size nt
                self.totalarea = no * 1.

                # 4) itterate over all timesteps and calculate weighting matrix
                for i in xrange(nt):
                    w[i, :, :] /= no[i]
            else:
                # 2) mask areas that do not contain valid data
                w = np.ma.array(w, mask=(w != w))
                self.totalarea = self.cell_area.sum()
                # normalization by total area. This does NOT result in sum(w)
                # == 1 for each timestep!
                w /= self.totalarea
            return w
        else:  # dimension
            raise ValueError(
                'weighting matrix not supported for this data shape')

    def areasum(self, return_data=False, apply_weights=True):
        """
        calculate area weighted sum of the spatial field for each time using area weights
        NOTE, that results must not be the same as from cdo fldsum(), as fldsum() DOES NOT
        perform an area weighting!

        Parameters
        ----------
        return_data : bool
            if True, then a C{Data} object is returned
        apply_weights : bool
            apply weights when calculating area weights

        Returns
        -------
        vector of spatial mean array[time]
        """

        if self.data.ndim == 3:
            pass
        elif self.data.ndim == 2:
            pass
        else:
            raise ValueError('Areasum currently only supported for 2D/3D data')

        if apply_weights:
            # area weighting
            # get weighting matrix for each timestep (taking care of invalid
            # data)
            w = self._get_weighting_matrix()
            # multiply the data with the weighting matrix in memory efficient
            # way
            w *= self.data
            if self.data.ndim == 3:
                w.shape = (len(self.data), -1)
                tmp = w.sum(axis=1)  # ... gives weighted sum
            elif self.data.ndim == 2:
                tmp = np.asarray([np.asarray(w.sum())])
            else:
                raise ValueError('Undefined!')

            # mean = sum { w * x } = sum { area * x / totalarea } ==> mean *
            # totalarea = sum {area * x}
            # this is the difference to fldmean() !; Here we rescale the result
            # by the total area used for calculating the weights
            tmp *= self.totalarea

        else:
            # no area weighting
            if self.data.ndim == 3:
                tmp = np.reshape(self.data, (len(self.data), -1)).sum(axis=1)
            elif self.data.ndim == 2:
                tmp = np.asarray([self.data.sum()])
            else:
                raise ValueError('Undefined')

        tmp = np.ma.array(tmp, mask=tmp != tmp)

        #////
        if return_data:  # return data object
            if self.data.ndim == 3:
                x = np.zeros((len(tmp), 1, 1))
                x[:, 0, 0] = tmp
            elif self.data.ndim == 2:
                x = np.zeros((1, 1))
                x[:, :] = tmp[0]
            else:
                raise ValueError('Undefined')

            assert (isinstance(tmp, np.ma.masked_array))
            r = self.copy()
            # use mask of array tmp (important if all values are invalid!)
            r.data = np.ma.array(x.copy(), mask=tmp.mask)

            # return cell area array with same size of data
            r.cell_area = np.array([1.])

            return r
        else:  # return numpy array
            return tmp

    def fldmean(self, return_data=False, apply_weights=True):
        """
        calculate mean of the spatial field for each time using weighted averaging
        results are exactly the same as one would obtain with the similar
        cdo function

        Parameters
        ----------
        return_data : bool
            if True, then a C{Data} object is returned

        apply_weights : bool
            apply weights when calculating area weights

        Returns
        -------
        r : ndarray or Data
            vector of spatial mean array[time]
        """

        if self.data.ndim == 3:
            pass
        elif self.data.ndim == 2:
            pass
        else:
            raise ValueError('fldmean currently only supported for 2D/3D data')

        if apply_weights:
            # area weighting
            # get weighting matrix for each timestep (taking care of invalid
            # data)
            w = self._get_weighting_matrix()
            # multiply the data with the weighting matrix in memory efficient
            # way
            w *= self.data
            if self.data.ndim == 3:
                w.shape = (len(self.data), -1)
                tmp = w.sum(axis=1)  # ... gives weighted sum
            elif self.data.ndim == 2:
                tmp = np.asarray([np.asarray(w.sum())])
            else:
                raise ValueError('Undefined!')
        else:
            # no area weighting
            if self.data.ndim == 3:
                tmp = np.reshape(self.data, (len(self.data), -1)).mean(axis=1)
            elif self.data.ndim == 2:
                tmp = np.asarray([self.data.mean()])
            else:
                raise ValueError('Undefined')

        if return_data:  # return data object
            if self.data.ndim == 3:
                x = np.zeros((len(tmp), 1, 1))
                x[:, 0, 0] = tmp
            elif self.data.ndim == 2:
                x = np.zeros((1, 1))
                x[:, :] = tmp[0]
            else:
                raise ValueError('Undefined')
            assert (isinstance(tmp, np.ma.masked_array))
            r = self.copy()
            r.data = np.ma.array(x.copy(),
                                 mask=tmp.mask)  # use mask of array tmp (important if all values are invalid!)

            # return cell area array with same size of data
            r.cell_area = np.array([1.])
            return r
        else:  # return numpy array
            return tmp

    def fldstd(self, return_data=False, apply_weights=True, ddof=0):
        """
        calculate stdv of the spatial field using area weighting
        returns exactly same results as the same CDO function

        Parameters
        ----------
        return_data : bool
            if True, then a C{Data} object is returned
        apply_weights : bool
            specifies if area weighting should be applied
        ddof : int, optional
            Means Delta Degrees of Freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
            By default `ddof` is zero.

        Returns
        -------
        res : ndarray, Data
            vector of spatial std array[time]

        Test
        ----
        unittest implemented

        Todo
        ----
        * proper unittest reference for ddof=1
        * unittest together with cdo's
        """

        if self.data.ndim == 3:
            pass
        elif self.data.ndim == 2:
            pass
        else:
            raise ValueError('fldstd currently only supported for 3D data')

        if ddof not in [0, 1]:
            raise ValueError('ddof only supported for [0,1] so far!')
        if ddof == 1:
            raise ValueError(
                'Sorry, but for DDOF=1 there are still problems for weighted samples. Please check unittests first.')

        if apply_weights:
            # calculate weighted standard deviation.
            # http://en.wikipedia.org/wiki/Mean_square_weighted_deviation
            #(adapted from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy)

            # in general it is assumed that the weights are normalized,
            # thus that sum(w) = 1., but the routine below is coded
            # in a way that this is not obligatory

            # calculate weighting matrix
            # get weighting matrix for each timestep (taking care of invalid
            # data)
            w = self._get_weighting_matrix()

            if self.data.ndim == 2:
                mu = (self.data * w).sum() / w.sum()
                if ddof == 0:  # unbiased estimator
                    V1 = w.sum()
                    s = (w * (self.data - mu) ** 2.).sum() / V1
                    tmp = [np.sqrt(s)]
                elif ddof == 1:  # biased estimator
                    V1 = w.sum()
                    V2 = (w * w).sum()
                    print 'V1, V2: ', V1, V2
                    print 'mu: ', mu
                    s = V1 * np.sum(w * (self.data - mu) ** 2.) / \
                        (V1 * V1 - V2)
                    tmp = [np.sqrt(s)]
                else:
                    raise ValueError('DDOF /= 1 not implemented yet!')

            elif self.data.ndim == 3:
                nt, ny, nx = self.shape
                s = np.ones(nt) * np.nan
                if w.ndim != 3:
                    raise ValueError('need mask for each timestep!')

                if ddof == 0:
                    for i in xrange(nt):
                        V1 = w[i, :, :].sum()
                        mu = (self.data[i, :, :] * w[i, :, :]).sum() / V1
                        s[i] = (w[i, :, :] * (self.data[i, :, :] - mu)
                                ** 2.).sum() / V1
                    tmp = np.sqrt(s)
                elif ddof == 1:
                    for i in xrange(nt):
                        V1 = w[i, :, :].sum()
                        mu = (self.data[i, :, :] * w[i, :, :]).sum() / V1
                        V2 = np.sum(w[i, :, :] ** 2.)
                        s[i] = V1 * (w[i, :, :] * (self.data[i, :, :] - mu)
                                     ** 2.).sum() / (V1 * V1 - V2)
                    tmp = np.sqrt(s)
                else:
                    raise ValueError('DDOF /= 1 not implemented yet!')
            else:
                raise ValueError('Undefined')

        else:
            # no area weighting
            if self.data.ndim == 2:
                tmp = self.data.std(ddof=ddof)
            elif self.data.ndim == 3:
                tmp = np.reshape(self.data, (len(self.data), -1)
                                 ).std(axis=1, ddof=ddof)
            else:
                raise ValueError('Undefined in fldstd()')

        tmp = np.ma.array(tmp, mask=tmp != tmp)

        if return_data:  # return data object
            if self.data.ndim == 3:
                x = np.zeros((len(tmp), 1, 1))
                x[:, 0, 0] = tmp
            elif self.data.ndim == 2:
                x = np.zeros((1, 1, 1))
                x[0, 0, 0] = tmp
            else:
                raise ValueError('Undefined')

            if not (isinstance(tmp, np.ma.masked_array)):
                print type(tmp)
                print self.data.ndim
                raise ValueError('Invalid data type!')

            r = self.copy()
            r.data = np.ma.array(x.copy(),
                                 mask=tmp.mask)  # use mask of array tmp (important if all values are invalid!)

            # return cell area array with same size of data
            r.cell_area = np.array([1.])

            return r
        else:  # return numpy array
            return tmp

    def _get_label(self):
        """
        return a nice looking label
        """
        if hasattr(self, 'label'):
            pass
        else:
            self.label = ''
        return self.label

#-----------------------------------------------------------------------

    def _convert_time(self):
        """
        convert time that was given as YYYYMMDD.f
        and set time variable of Data object
        """
        s = map(str, self.time)
        T = []
        for t in s:
            y = t[0:4]
            m = t[4:6]
            d = t[6:8]
            h = t[8:]
            h = str(int(float(h) * 24.))
            tn = y + '-' + m + '-' + d + ' ' + h + ':00'
            T.append(tn)
        T = np.asarray(T)
        self.calendar = 'gregorian'
        self.time_str = 'days since 0001-01-01 00:00:00'
        # convert first to datetime object and then use own function !!!
        self.time = self.date2num(plt.num2date(plt.datestr2num(T)))

    def _convert_time_YYYYMMDD(self):
        """
        convert time that was given as YYYYMMDD
        and set time variable of Data object
        """
        s = map(str, self.time)
        T = []
        for t in s:
            y = t[0:4]
            m = t[4:6]
            d = t[6:8]  # always the first day is used as default
            h = '00'
            tn = y + '-' + m + '-' + d + ' ' + h + ':00'
            T.append(tn)
        T = np.asarray(T)
        self.calendar = 'gregorian'
        self.time_str = 'days since 0001-01-01 00:00:00'
        # convert first to datetime object and then use own function !!!
        self.time = self.date2num(plt.num2date(plt.datestr2num(T)))

    def _convert_timeYYYYMM(self):
        """
        convert time that was given as YYYYMM.f
        and set time variable of Data object
        """
        s = map(str, self.time)
        T = []
        for t in s:
            y = t[0:4]
            m = t[4:6]
            d = '01'  # always the first day is used as default
            h = '00'
            tn = y + '-' + m + '-' + d + ' ' + h + ':00'
            T.append(tn)
        T = np.asarray(T)
        self.calendar = 'gregorian'
        self.time_str = 'days since 0001-01-01 00:00:00'
        # convert first to datetime object and then use own function !!!
        self.time = self.date2num(plt.num2date(plt.datestr2num(T)))

    def _convert_timeYYYY(self):
        """
        convert time that was given as YYYY.f
        and set time variable of Data object

        The date is set to the first of January for each year
        """

        s = map(str, self.time)
        T = []
        for t in s:
            y = t[0:4]
            m = '01'
            d = '01'  # always the first day is used as default
            h = '00'
            tn = y + '-' + m + '-' + d + ' ' + h + ':00'
            T.append(tn)
        T = np.asarray(T)
        self.calendar = 'gregorian'
        self.time_str = 'days since 0001-01-01 00:00:00'
        # convert first to datetime object and then use own function !!!
        self.time = self.date2num(plt.num2date(plt.datestr2num(T)))

    def adjust_time(self, day=None, month=None, year=None):
        """
        correct all timestamps and assign same day and/or month
        for all timesteps

        Parameters
        ----------
        day : int
            if specified then the argument will be used as day for
            all timestamps
        month : int
            if specified then the argument will be used as month for
            all timestamps
        year : int
            if specified then the argument will be used as year for
            all timestamps

        Test
        ----
        unittest implemented
        """

        o = []
        for t in self.time:
            d = self.num2date(t)
            s = str(d)  # convert to a string
            if day is not None:
                s = s[0:8] + str(day).zfill(2) + s[10:]  # replace day
            if month is not None:
                s = s[0:5] + str(month).zfill(2) + s[7:]  # replace month
            if year is not None:
                s = str(year).zfill(4) + s[4:]

            # convert str. a number and then again to a datetime object
            # to allow to employ specific time conversion of data object
            o.append(self.date2num(plt.num2date(plt.datestr2num(s))))

        o = np.asarray(o)
        self.time = o.copy()

#-----------------------------------------------------------------------

    def timsort(self, return_object=False):
        """
        sorts a C{Data} object in accordance with its time axis.

        A typical application of this function would be for climatological mean
        values. If one calculates a climatology using the cdo ymonmean command,
        it is *not* guarantued that the data is actually in ascending order, namely, that
        January is the first dataset. The reason is, that the cdo's start with the dataset
        of the first month!

        by using timsort() one can ensure a proper sequence of the data.
        In case of a climatology, it is recommended that you first set the day and year to a common
        date to get in the end a sorting of the months. An example would look like

        self.adjust_time(day=15,year=2000) #sets year to a dummy = 2000 and day = 15
        self.timsort() #results in a sorted climatology

        Parameters
        ----------
        return_object : bool
            specifies if a Data object shall be returned

        Test
        ----
        unittest implemented
        """

        # checks
        if self.time is None:
            raise ValueError('Time array needed for timsort()')
        if self.data.ndim != 3:
            raise ValueError('3D array needed for timsort()')

        # specify object ot work on
        if return_object:
            x = self.copy()
        else:
            x = self

        # do the sorting
        s = np.argsort(x.time)
        x.data = x.data[s, :, :]
        x.time = x.time[s]
        if hasattr(x, 'std'):  # standard deviation
            x.std = x.std[s, :, :]
        if hasattr(x, 'n'):  # number of datasets
            x.n = x.n[s, :, :]

        # result
        if return_object:
            return x

    def get_aoi(self, region):
        """ region of class Region """

        # copy self
        d = self.copy()
        d.data = region.get_subset(d.data)
        d.cell_area = region.get_subset(d.cell_area)

        if hasattr(d, '_climatology_raw'):
            d._climatology_raw = region.get_subset(d._climatology_raw)

        if plt.isvector(d.lat):
            d.lat = d.lat[region.y1:region.y2]
        else:
            d.lat = region.get_subset(d.lat)

        if plt.isvector(d.lon):
            d.lon = d.lon[region.x1:region.x2]
        else:
            d.lon = region.get_subset(d.lon)

        d.label = d.label + ' (' + region.label + ')'
        return d

    def get_aoi_lat_lon(self, R, apply_mask=True):
        """
        get area of interest (AOI) given lat/lon coordinates

        the routine masks all area which is NOT in the given area
        coordinates of region are assumed to be in -180 < lon < 180

        CAUTION: the current object will be changed

        Parameters
        ----------
        R : Region
            region object that specifies region
        apply_mask : bool
            apply former data mask (default)
        """

        LON = self.lon.copy()
        if self._lon360:
            tmsk = LON > 180.
            LON[tmsk] -= 360.
            del tmsk
        msk_lat = (self.lat >= R.latmin) & (self.lat <= R.latmax)
        msk_lon = (LON >= R.lonmin) & (LON <= R.lonmax)

        # additional mask in Region object
        if (R.mask is None) | (apply_mask is False):
            msk_region = np.ones(np.shape(msk_lat)).astype('bool')
        else:
            if np.shape(msk_lat) != np.shape(R.mask):
                print np.shape(msk_lat), np.shape(R.mask)
                raise ValueError('Invalid geometries for mask')
            else:
                msk_region = R.mask
        msk = msk_lat & msk_lon & msk_region  # valid area
        self._apply_mask(msk)

    def cut_bounding_box(self, return_object=False):
        """
        estimate bounding box of data and subset dataset such that
        only valid data is contained in the bounding box.

        In the current setup, only the part of the grid is returned
        where *all* timesteps are valid

        Parameters
        ----------
        return_object : bool
            return data object, otherwise the modifications are
            applied to current object
        """

        # get bounding box
        # note that the indices can not be used directly for array
        # slicing. One tyipically needs to add '1' to the last index
        i1, i2, j1, j2 = self.get_bounding_box()

        if return_object:
            D = self.copy()
        else:
            D = self

        if self.ndim == 3:
            D.data = D.data[:, i1:i2 + 1, j1:j2 + 1]
        elif self.ndim == 2:
            D.data = D.data[i1:i2 + 1, j1:j2 + 1]
        else:
            raise ValueError(
                'Cutting of bounding box not implemented for data other than 2D/3D!')
        if hasattr(self, 'lat'):
            if D.lat is not None:
                D.lat = D.lat[i1:i2 + 1, j1:j2 + 1]
        if hasattr(self, 'lon'):
            if D.lon is not None:
                D.lon = D.lon[i1:i2 + 1, j1:j2 + 1]
        if hasattr(D, 'cell_area'):
            if D.cell_area is not None:
                D.cell_area = D.cell_area[i1:i2 + 1, j1:j2 + 1]

        if return_object:
            return D
        else:
            return None

    def get_valid_mask(self, frac=1., return_frac=False):
        """
        calculate a mask which is True, when a certain fraction of
        all timestamps of the field are valid

        Parameters
        ----------
        frac : float
            minimum fraction of valid data in all timesteps needed [0..1]
        return_frac : bool
            return also fraction and not only final mask

        Returns
        -------
        msk : ndarray
            array with mask where valid data is available for at least
            *frac* timesteps
        frac : ndarray
            optional: fraction

        Test
        ----
        unittest implemented
        """

        if (frac < 0.) or (frac > 1.):
            raise ValueError('Fraction needs to be between 0 ... 1!')

        if self.data.ndim == 1:
            thefrac = 1. - float(sum(self.data.mask)) / \
                float(len(self.data.mask))  # valid fraction
            if thefrac >= frac:
                msk = np.ones((1, 1)).astype('bool')
            else:
                msk = np.zeros((1, 1)).astype('bool')
            if return_frac:
                return msk, thefrac
            else:
                return msk
        elif self.data.ndim == 2:
            # keeps original data mask
            res = np.ones(self.data.shape).astype('bool')
            res[self.data.mask] = False
            if return_frac:
                thefrac = np.ones(self.data.shape)
                return res, thefrac
            else:
                return res
        elif self.data.ndim == 3:
            n = len(self.data)  # number of timesteps
            hlp = self.data.copy()
            if hasattr(hlp, 'mask'):
                hlp1 = hlp.data.copy()
                hlp1[hlp.mask] = np.nan
            else:
                hlp1 = hlp.data.copy()
            thefrac = (np.sum(~np.isnan(hlp1), axis=0) / float(n))
            msk = thefrac >= frac
            if return_frac:
                return msk, thefrac
            else:
                return msk
        else:
            raise ValueError('Unsupported dimension!')

    def get_valid_data(self, return_mask=False, mode='all'):
        """
        this routine calculates from the masked array
        only the valid data and returns it together with its
        coordinate as vector

        valid means that ALL timestamps need to be valid!

        Parameters
        ----------
        return_mask : bool
            specifies if the mask applied to the original data should
            be returned as well
        mode : str
            analysis mode ['all','one']
            'all': all timestamps need to be valid
            'one': at least a single dataset needs to be valid
        """

        if hasattr(self, 'lon'):
            if self.lon is not None:
                lon = self.lon.reshape(-1)
            else:
                lon = None
        else:
            lon = None
        if hasattr(self, 'lat'):
            if self.lat is not None:
                lat = self.lat.reshape(-1)
            else:
                lat = None
        else:
            lat = None

        if self.ndim == 3:

            n = len(self.time)

            # vectorize the data

            data = self.data.reshape(n, -1)
            # set pixels with NaN to invalid
            data.mask[np.isnan(data.data)] = True

            # extract only valid (not masked data)
            if mode == 'all':
                # identify all ngrid cells where all timesteps are valid
                msk = np.sum(~data.mask, axis=0) == n
            elif mode == 'one':
                # identify ONE grid cell where all timesteps are valid
                msk = np.sum(~data.mask, axis=0) > 0
            else:
                raise ValueError('Invalid option in get_valid_data() %s' %
                                 mode)

            data = data[:, msk]

        elif self.ndim == 2:
            data = self.data.reshape(-1)
            msk = ~data.mask
            data = data[msk]
        else:
            raise ValueError('Unsupported dimension!')

        if lon is not None:
            lon = lon[msk]
        if lat is not None:
            lat = lat[msk]
        if return_mask:
            return lon, lat, data, msk
        else:
            return lon, lat, data

    def _apply_mask(self, msk1, keep_mask=True):
        """
        apply a mask to C{Data}. All data where mask==True
        will be masked. Former data and mask will be stored.

        When a Data object is provided as a mask, the mask
        attribute of the data field will be used for mask identification

        Parameters
        ----------
        msk1 : ndarray or Data object
            mask to be applied to data. Needs to have same geometry as
            data.
        keep_mask : bool
            keep old mask
        """

        if isinstance(msk1, Data):
            msk = msk1.data.mask
        else:
            msk = msk1

        self.__oldmask = self.data.mask.copy()
        self.__olddata = self.data.data.copy()
        if hasattr(self, 'std'):
            if self.data.shape != self.std.shape:
                raise ValueError(
                    'Standard deviation has different geometry than data!')
            self.__oldstd = self.std.data.copy()

        if self.data.ndim == 2:
            # convert to float to allow for nan support
            tmp1 = self.data.copy().astype('float')
            tmp1[~msk] = np.nan

            if hasattr(self, 'std'):
                tmps = self.std.copy().astype('float')

            if keep_mask:
                if self.__oldmask.ndim > 0:
                    tmp1[self.__oldmask] = np.nan
                    if hasattr(self, 'std'):
                        tmps[self.__oldmask] = np.nan

            self.data = np.ma.array(tmp1, mask=np.isnan(tmp1))
            if hasattr(self, 'std'):
                self.std = np.ma.array(tmps, mask=np.isnan(tmps))
                del tmps
            del tmp1

        elif self.data.ndim == 3:
            for i in xrange(self.nt):
                tmp = self.data[i, :, :].copy()
                tmp[~msk] = np.nan
                self.data[i, :, :] = tmp[:, :] * 1.
                del tmp

                if hasattr(self, 'std'):
                    tmps = self.std[i, :, :].copy()
                    tmps[~msk] = np.nan
                    self.std[i, :, :] = tmps
                    del tmps

            if keep_mask:
                if self.__oldmask.ndim > 0:
                    self.data.data[self.__oldmask] = np.nan
                    if hasattr(self, 'std'):
                        self.std.data[self.__oldmask] = np.nan

            self.data = np.ma.array(
                self.data.data, mask=np.isnan(self.data.data))
            if hasattr(self, 'std'):
                self.std = np.ma.array(
                    self.std.data, mask=np.isnan(self.std.data))
        else:
            print np.shape(self.data)
            raise ValueError('Unsupported geometry _apply_mask')

        if hasattr(self, '_climatology_raw'):
            for i in range(len(self._climatology_raw)):
                tmp = self._climatology_raw[i, :, :].copy()
                tmp[~msk] = np.nan
                self._climatology_raw[i, :, :] = tmp[:, :]
                del tmp

    def shift_x(self, nx):
        """
        shift data array in x direction by nx steps

        Parameters
        ----------
        nx : int
            shift by nx steps
        """
        if self.data.ndim == 3:
            self.data = self.__shift3D(self.data, nx)
        else:
            self.data = self.__shift2D(self.data, nx)
        self.lat = self.__shift2D(self.lat, nx)
        self.lon = self.__shift2D(self.lon, nx)

    def __shift3D(self, x, n):
        """
        shift 3D data

        Parameters
        ----------
        x : ndarray
            data to be shifted
        n : int
            shifting step
        """
        tmp = x.copy()
        y = x.copy()
        y[:, :, :] = np.nan
        y[:, :, 0:n] = tmp[:, :, -n:]
        y[:, :, n:] = tmp[:, :, 0:-n]
        return y

    def timeshift(self, n, return_data=False, shift_time=False):
        """
        shift data in time by n-steps
        positive numbers mean, that the
        data is shifted leftwards (thus towards earlier times)

        e.g. timeseries 1980,1981,1982, a shift of n=1 would generate
        a dataset with the data ordered as follows
        [1981,1982,1980]

        Note that the timevector remains UNCHANGED!

        Parameters
        ----------
        n : int
            lag to shift data (n>=0)
        shift_time : bool
            shift also timevector

        return_data : bool
            specifies if a NEW C{Data} object shall be returned

        Tests
        -----
        unittests implemented
        """

        if n == 0:
            return
        if n < 0:
            raise ValueError('Negative timeshifts are not supported yet!')
        if self.data.ndim != 3:
            raise ValueError(
                'array of size [time,ny,nx] is needed for temporal shifting!')

        # copy data
        tmp = self.data.copy()

        if return_data:
            res = self.copy()
        else:
            res = self

        res.data[:, :, :] = np.nan
        res.data[:-n:, :, :] = tmp[n:, :, :]
        res.data[-n:, :, :] = tmp[0:n, :, :]
        res.data = np.ma.array(res.data, mask=np.isnan(res.data))
        del tmp

        if shift_time:  # shift also timevector
            tmp = res.time * 1.
            res.time[:-n:] = tmp[n:]
            res.time[-n:] = tmp[0:n]

        if return_data:
            return res
        else:
            return None

    def _set_valid_range(self, vmin, vmax):
        """
        sets the valid range of the data
        only data with vmin <= data <= vmax will be kept as valid

        Parameters
        ----------
        vmin : float
            minimum valid value
        vmax : float
            maximum valid value

        Tests
        -----
        unittest implemented
        """
        self.data = np.ma.array(
            self.data, mask=((self.data < vmin) | (self.data > vmax)))

    def __shift2D(self, x, n):
        """
        shift 2D data

        Parameters
        ----------
        x : Data
            data to be shifted
        n : int
            shifting step
        """
        tmp = x.copy()
        y = x.copy()
        y[:, :] = np.nan
        y[:, 0:n] = tmp[:, -n:]
        y[:, n:] = tmp[:, 0:-n]
        return y

    def copy(self):
        """
        copy complete C{Data} object including all attributes
        """
        d = Data(None, None)

        for attr, value in self.__dict__.iteritems():
            try:
                # copy (needed for arrays)
                cmd = "d." + attr + " = self." + attr + '.copy()'
                exec cmd
            except:
                # copy
                cmd = "d." + attr + " = self." + attr
                exec cmd
        return d

#-----------------------------------------------------------------------

    def add(self, x, copy=True):
        """
        Add a Data object to the current object field

        Parameters
        ----------
        x : Data
            object to be added to the current object
        copy : bool
            if True, then a copy is returned. Otherwise the actual data
            is modified

        Test
        ----
        unittest implemented
        """

        if np.shape(self.data) != np.shape(x.data):
            raise ValueError('Inconsistent geometry (add): can not calculate!')

        if copy:
            d = self.copy()
        else:
            d = self
        d.data = d.data + x.data
        d.label = self.label + ' + ' + x.label
        return d

    def sub(self, x, copy=True):
        """
        Substract a C{Data} object from the current object field

        Parameters
        ----------
        x : Data
            object to be substracted from the current object
        copy : bool
            if True, then a copy is returned. Otherwise the actual data
            is modified

        Test
        ----
        unittest implemented
        """

        f_elementwise = False
        if np.shape(self.data) != np.shape(x.data):
            s1 = np.shape(self.data)
            if (s1[1] == x.data.shape[0]) & (s1[2] == x.data.shape[1]):
                f_elementwise = True
            else:
                sys.stdout.write(str(self.shape))
                sys.stdout.write(str(x.shape))
                raise ValueError(
                    'Inconsistent geometry (sub): can not calculate!')

        if copy:
            d = self.copy()
        else:
            d = self
        if f_elementwise:
            for i in xrange(len(d.data)):
                d.data[i, :, :] = d.data[i, :, :] - x.data[:, :]
        else:
            d.data = d.data - x.data
        d.label = self.label + ' - ' + x.label
        return d

    def diff(self, x, axis=0, equal_var=True, mask_data=False,
             pthres=0.05):
        """
        Difference between two C{Data} objects

        This routine calculates the difference between the data of two
        datasets. It calculates the significance and returns
        the mean differences and their corresponding significances

        The significance is calculated using a two-tailored t-test or a welch test
        in case of different variances. Independent samples are assumed!
        (unittest test_diff)

        Parameters
        ----------
        x : Data
            Data object which will be substracted from self
        axis : int
            axis along which the data will be aggregated
            (typically axis=0 corresponds to time)
        equal_var : bool
            specifies if the two input datasets (self,x) are expected
            to have same variance. Dependent on this parameter.
            If the variance is equal, then a t-test is applied,
            if not, then a welch test is used.
        mask_data : bool
            specifies if the data field which is returned should be
            masked for all areas that do *not* show significant changes!
        pthres : float
            threshold for significant p-value; a value of e.g. 0.05 corresponds to the 95% significance level.
            this threshold will be used for the generation of a mask that might be used e.g. as an overlay in map_plot()

        Returns
        -------
        returns a C{Data} object that includes a) the difference map,
        b) the p-value, c) a mask that can be used e.g. as an overlay
        for map_plot()

        @todo: implementation of welch test. This should be actually already be implemented in stats.ttest_ind, but is not available in my python installation!
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html#scipy.stats.ttest_ind
        but the code would be available here: https://github.com/scipy/scipy/blob/v0.11.0/scipy/stats/stats.py#L2957
        """

        #/// check consistency
        if np.shape(self.data) != np.shape(x.data):
            raise ValueError('Inconsistent geometry (sub): can not calculate!')
        if axis > len(self.data.shape) - 1:
            raise ValueError('Invalid axis parameter: %s' % str(axis))
        if axis < 0:
            raise ValueError('Invalid axis parameter: %s' % str(axis))

        #/// create new data object
        d = self.copy()
        d.label = self.label + ' - ' + x.label

        #/// calculate statistical significance of the difference
        if isinstance(d.data, np.ma.masked_array):
            sta = stats.mstats
            # use routine in pyCMBS.statistic.py
            t, p = ttest_ind(d.data, x.data, axis=axis)
        else:
            t, p = stats.ttest_ind(d.data, x.data,
                                   axis=axis)  # todo equal var for welch test not part of my psthon installation!

        # invert p-value, as a p-value of 1. would correspond to the same data
        p = 1. - p

        #/// mean difference masked if p-value too low
        mask = p <= pthres
        if mask_data:
            # mean difference as masked array
            d.data = np.ma.array(self.timmean() - x.timmean(), mask=~mask)
        else:
            d.data = np.ma.array(self.timmean() - x.timmean(),
                                 mask=np.zeros(self.timmean().shape).astype('bool'))  # mean difference as masked array
        d.p_value = p
        # masks the grid cells that show significant changes (todo check this
        # again!) needs additional validation
        d.p_mask = mask
        d.t_value = t

        return d

    def subc(self, x, copy=True):
        """
        Substract a constant value from the current object field

        Parameters
        ----------
        x : float
            constant (can be either a scalar or a field that has
            the same geometry as the second and third dimension of self.data)
        copy : bool
            if True, then a new data object is returned
            else, the data of the present object is changed
        """
        if copy:
            d = self.copy()
        else:
            d = self
        if np.isscalar(x):
            d.data -= x
        elif x.ndim == 2:  # x is an array
            for i in xrange(len(self.time)):
                d.data[i, :, :] -= x[:, :]
        else:
            raise ValueError('Invalid geometry in detrend()')
        return d

    def addc(self, x, copy=True):
        """
        Add a constant value to the current object field

        Parameters
        ----------
        x : float
            constant value to add to current object
        copy : bool
            if True, then a new data object is returned
            else, the data of the present object is changed

        Test
        ----
        unittest implemented
        """

        if copy:
            d = self.copy()
        else:
            d = self
        d.data += x
        return d

    def mulc(self, x, copy=True):
        """
        Multiply current data by a constant

        Parameters
        ----------
        x : float
            constant value to multiply with current object
        copy : bool
            if True, then a new data object is returned
            else, the data of the present object is changed

        Test
        ----
        unittest implemented
        """

        if copy:
            d = self.copy()
        else:
            d = self
        d.data *= x
        return d

    def divc(self, x, copy=True):
        """
        Divide current data by a constant

        Parameters
        ----------
        x : float
            constant value to divide current object by
        copy : bool
            if True, then a new data object is returned
            else, the data of the present object is changed

        Test
        ----
        unittest implemented
        """

        if copy:
            d = self.copy()
        else:
            d = self
        d.data /= x
        return d

#-----------------------------------------------------------------------

    def div(self, x, copy=True):
        """
        Divide current object field by field of a C{Data} object

        Parameters
        ----------
        x : Data
            object in the denominator, (data needs to have either same
            geometry as self.data or second and third dimension need to match)
        copy : bool
            if True, then a new data object is returned
            else, the data of the present object is changed

        Test
        ----
        unittest implemented
        """

        if np.shape(self.data) != np.shape(x.data):
            if self.data.ndim == 3:
                if x.data.ndim == 2:
                    if np.shape(self.data[0, :, :]) == np.shape(x.data):
                        # second and third dimension match
                        pass
                    else:
                        print np.shape(self.data)
                        print np.shape(x.data)
                        raise ValueError(
                            'Inconsistent geometry (div): can not calculate!')
                else:
                    print np.shape(self.data)
                    print np.shape(x.data)
                    raise ValueError(
                        'Inconsistent geometry (div): can not calculate!')
            else:
                print np.shape(self.data)
                print np.shape(x.data)
                raise ValueError(
                    'Inconsistent geometry (div): can not calculate!')

        if copy:
            d = self.copy()
        else:
            d = self
        if np.shape(d.data) == np.shape(x.data):
            d.data /= x.data
        elif np.shape(d.data[0, :, :]) == np.shape(x.data):
            for i in xrange(len(self.time)):
                d.data[i, :, :] /= x.data
        else:
            raise ValueError('Can not handle this geometry in div()')

        d.label = self.label + ' / ' + x.label

        return d

    def mul_tvec(self, x, copy=True):
        """
        multiply the data with a time vector.
        Each timestep will be multiplied with the corresponding
        element from the time vector

        out[i,:,:] = in[i,:,:] * x[i]

        Parameters
        ----------
        x : ndarray
            data vector with scaling constants for each timestep
            len(x) = self.nt
        copy : bool
            if True, a Data object is returned. Otherwise changes are
            applied to the input data object (self)
        """

        if x.ndim != 1:
            raise ValueError('Only 1D arrays allowed')
        if len(x) != self.nt:
            print len(x), self.nt
            raise ValueError('Inconsistent geometries for temporal vector!')

        if copy:
            o = self.copy()
        else:
            o = self

        if False:
            # slow implementation ...
            for i in xrange(self.nt):
                o.data[i, :, :] = self.data[i, :, :] * x[i]
        else:
            # fast implementation using broadcasting ...
            # http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html
            o.data[:, :, :] = self.data[:, :, :] * x[:, np.newaxis, np.newaxis]

        if copy:
            return o
        else:
            return None

    def mul(self, x, copy=True):
        """
        Multiply current object field by field by a C{Data} object

        Parameters
        ----------
        x : Data
            object in the denominator, (data needs to have either same
            geometry as self.data or second and third dimension need to match)
        copy : bool
            if True, then a new data object is returned
            else, the data of the present object is changed

        Test
        ----
        unittest implemented
        """

        if np.shape(self.data) != np.shape(x.data):
            if self.data.ndim == 3:
                if x.data.ndim == 2:
                    if np.shape(self.data[0, :, :]) == np.shape(x.data):
                        # second and third dimension match
                        pass
                    else:
                        print np.shape(self.data)
                        print np.shape(x.data)
                        raise ValueError(
                            'Inconsistent geometry (div): can not calculate!')
                else:
                    print np.shape(self.data)
                    print np.shape(x.data)
                    raise ValueError(
                        'Inconsistent geometry (div): can not calculate!')
            else:
                print np.shape(self.data)
                print np.shape(x.data)
                raise ValueError(
                    'Inconsistent geometry (div): can not calculate!')

        if copy:
            d = self.copy()
        else:
            d = self
        if np.shape(d.data) == np.shape(x.data):
            d.data = d.data * x.data
        elif np.shape(d.data[0, :, :]) == np.shape(x.data):
            for i in xrange(len(self.time)):
                d.data[i, :, :] = d.data[i, :, :] * x.data
        else:
            raise ValueError('Can not handle this geometry in div()')

        d.label = self.label + ' * ' + x.label
        return d

    def _sub_sample(self, step):
        """
        perform spatial subsampling of data

        Parameters
        ----------
        step : int
            stepsize for subsampling
        """
        if self.data.ndim == 3:
            self.data = self.data[:, ::step, ::step]
        elif self.data.ndim == 2:
            self.data = self.data[::step, ::step]
        else:
            raise ValueError('Data Dimension not supported!')
        if hasattr(self, 'lat'):
            if self.lat is not None:
                self.lat = self.lat[::step, ::step]
        if hasattr(self, 'lon'):
            if self.lon is not None:
                self.lon = self.lon[::step, ::step]

    def corr_single(self, x, pthres=1.01, mask=None, method='pearson'):
        """
        The routine correlates a data vector with all data of the
        current object. It returns several correlation measures

        Example
        -------
        >> d = Data(None, None)
        >> x = np.random(100)
        >> rpears,slope,intercept,p,covar= d.corr_single(x, pthres=0.05)

        Parameters
        ----------
        x : ndarray
            the data vector correlations should be calculated with,
            numpy array [time]
        method : str
            correlation method to be used ['spearman','pearson']
        pthres : float
            significance threshold. All values below this
            threshold will be returned as valid
        mask : ndarray
            mask to flag invalid data

        Test
        ----
        unittest implemented
        """

        if method not in ['pearson', 'spearman']:
            raise ValueError(
                'Only pearson or spearman rank correlation supported so far.')
        if self.ndim != 3:
            raise ValueError('Invalid geometry!')

        nt, ny, nx = sz = self.shape

        if nt != len(x):
            raise ValueError('Inconsistent geometries')

        # check if 'x' is a masked array where the mask has the same
        # size as the data. If this is not the case, then set it
        if isinstance(x, np.ma.core.MaskedArray):
            if isinstance(x.mask, np.ndarray):
                if x.mask.shape != x.data.shape:
                    raise ValueError('Invalid mask geometry!')
            else:
                raise ValueError(
                    'The mask of the dataset needs to be an array!')
        else:
            raise ValueError('Expect masked array as input in corr_single')

        # get data with at least one valid value
        lo, la, dat, msk = self.get_valid_data(return_mask=True, mode='one')
        xx, n = dat.shape
        if self.verbose:
            print('   Number of grid points: ', n)

        R = np.ones((ny, nx)) * np.nan  # output matrix for correlation
        P = np.ones((ny, nx)) * np.nan  # output matrix for p-value
        S = np.ones((ny, nx)) * np.nan  # output matrix for slope
        I = np.ones((ny, nx)) * np.nan  # output matrix for intercept
        CO = np.ones((ny, nx)) * np.nan  # output matrix for covariance

        R.shape = (-1)
        S.shape = (-1)
        P.shape = (-1)
        I.shape = (-1)
        CO.shape = (-1)

        print 'Calculating correlation ...'
        if method == 'pearson':
            res = [stats.mstats.linregress(x, dat[:, i]) for i in xrange(n)]
            res = np.asarray(res)
            slope = res[:, 0]
            intercept = res[:, 1]
            r_value = res[:, 2]
            p_value = res[:, 3]
            std_err = res[:, 4]

        elif method == 'spearman':
            res = np.ones((n, 5)) * np.nan

            # better implementation
            #~ if n < 3: not so easy, as 'n' is the total number of points ???
            #~ set results as invalid
            #~ else:
            #~ res = [stats.mstats.spearmanr(x, dat[:, i]) for i in xrange(n)]
            #~ ...

            # this is implemented like this at the moment, as the number of
            # valid data points needs to be > 3
            for i in xrange(n):
                invalid = False
                if (~x.mask).sum() < 3:
                    invalid = True
                if (~dat[:, i].mask).sum() < 3:
                    invalid = True
                    # do processing only for at least 3 samples!
                if invalid:
                    res[i, :] = np.nan  # set all to nan
                    continue
                else:
                    # todo: implement it more efficiently
                    rho, prob = stats.mstats.spearmanr(x, dat[:, i])
                    res[i, 0] = np.nan  # slope
                    res[i, 1] = np.nan  # intercept
                    res[i, 2] = rho  # r_value
                    res[i, 3] = prob  # p-value
                    res[i, 4] = np.nan  # std_err
            slope = res[:, 0]
            intercept = res[:, 1]
            r_value = res[:, 2]
            p_value = res[:, 3]
            std_err = res[:, 4]

        else:
            raise ValueError('Invalid method!')

        R[msk] = r_value
        P[msk] = p_value
        I[msk] = intercept
        S[msk] = slope
        R.shape = (ny, nx)
        P.shape = (ny, nx)
        I.shape = (ny, nx)
        S.shape = (ny, nx)

        #--- prepare output data objects
        Rout = self.copy()  # copy object to get coordinates
        Rout.label = 'correlation'
        msk = (P > pthres) | (np.isnan(R))
        #msk = np.zeros_like(R).astype('bool')
        Rout.data = np.ma.array(R, mask=msk).copy()
        Rout.unit = '-'

        Sout = self.copy()  # copy object to get coordinates
        Sout.label = 'slope'
        Sout.data = np.ma.array(S, mask=msk).copy()
        Sout.unit = self.unit

        Iout = self.copy()  # copy object to get coordinates
        Iout.label = 'intercept'
        Iout.data = np.ma.array(I, mask=msk).copy()
        Iout.unit = self.unit

        Pout = self.copy()  # copy object to get coordinates
        Pout.label = 'p-value'
        Pout.data = np.ma.array(P, mask=msk).copy()
        Pout.unit = '-'

        Cout = self.copy()  # copy object to get coordinates
        Cout.label = 'covariance'
        # currently not supported: covariance!
        Cout.data = np.ma.array(np.ones(P.shape) * np.nan, mask=msk).copy()
        Cout.unit = '-'

        if mask is not None:
            # apply a mask
            Rout._apply_mask(mask)
            Sout._apply_mask(mask)
            Iout._apply_mask(mask)
            Pout._apply_mask(mask)
            Cout._apply_mask(mask)

        return Rout, Sout, Iout, Pout, Cout

    def detrend(self, return_object=True):
        """
        detrend data timeseries by removing linear trend over time.
        It is assumed that the timesamples have equidistant spacing.
        This assumption is important to consider, as regression is calculated
        only using the length of the data vector!

        Parameters
        ----------
        return_object : bool
            specifies if C{Data} object will be returned (default=True)

        Returns
        -------
        if specified the routine returns a Data object. Otherwise, the
        current Data object is modified.

        Tests
        -----
        unittests implemented
        """

        print('Detrending data ...')

        if self.data.ndim != 3:
            raise ValueError('Can not detrend data other than 3D!')

        # generate dummy vector for linear correlation (assumes equally spaced
        # data!!!!!) todo: generate unittest for this
        # @todo: replace this by using actual timestamp for regression calcuclation
        x = np.arange(len(self.time))
        x = np.ma.array(x, mask=x != x)

        # correlate and get slope and intercept
        Rout, Sout, Iout, Pout, Cout = self.corr_single(x)

        # calculate regression field
        reg = Data(None, None)
        reg.data = np.zeros(self.data.shape) * np.nan
        reg.label = 'trend line'
        nt = self.data.shape[0]
        for i in range(nt):
            reg.data[i, :, :] = Sout.data * i + Iout.data

        # substract regression line
        res = self.sub(reg)
        res.label = self.label + '(det.)'

        if return_object:
            res.detrended = True
            return res
        else:
            self.data = res.data
            self.detrended = True
            return None

    def _is_daily(self):
        """
        check if the timeseries is daily

        Test
        ----
        unittest implemented
        """
        if hasattr(self, 'time'):
            d = np.ceil(pl.date2num(self.date))
            return np.all(np.diff(d) == 1.)
        else:
            return False

    def _is_monthly(self):
        """
        check if the data is based on a sequence of increasing monthly
        values. The routine simply checks of the months of the
        timeseries is increasing. Days are not considered!

        Returns
        -------
        bool

        Test
        ----
        unittest implemented
        """
        if hasattr(self, 'time'):
            # get list of all months
            mo = self._get_months()

            # get value of unique differences between months; only values 1 and
            # -11 are allowed
            di = np.unique(np.diff(mo))
            if len(di) > 2:
                return False
            for d in di:
                if d not in [1, -11]:
                    return False
                # if we have reached this point, then the months are in
            # ascending monthly order. Now check if the years are as well
            di = np.unique(np.diff(self._get_years()))
            for d in di:
                if d not in [0, 1]:
                    return False
                #... everything is o.k., we have an increasing monthly and yearly timeseries
            return True
        else:
            return False

    def _pad_timeseries(self, fill_value=-99.):

        import numpy as np
        from matplotlib import dates
        from dateutil.rrule import rrule, MONTHLY

        self._log_warning('Trying to pad timeseries')

        data = self.data
        time = self.time

        dummy_data = np.ma.array(
            np.ones(data[0, :, :].shape) * data.fill_value,
            fill_value=data.fill_value,
            mask=np.ones(data[0, :, :].shape) * True)

        months = self._get_months()
        mondif = np.diff(months)

        gaps = np.where(mondif > 1)
        new_time = self.num2date(time)
        new_data = data.copy()

        idx_shift = 0

        for i in gaps[0]:
            start_time = new_time[i + idx_shift]
            stop_time = new_time[i + 1 + idx_shift]
            gap_months = rrule(MONTHLY, dtstart=start_time).between(
                start_time, stop_time, inc=True)[1:-1]

            new_time = np.insert(new_time, i + 1 + idx_shift, gap_months)
            new_data = np.insert(
                new_data, i + 1 + idx_shift, dummy_data, axis=0)
            idx_shift = idx_shift + len(gap_months)

        data_masked = np.ma.array(new_data.data,
                                  mask=new_data == fill_value,
                                  fill_value=fill_value)

        self.time = self.date2num(new_time)
        self.data = data_masked.copy()
        self._set_timecycle()

    def _set_timecycle(self):
        """
        determine automatically the timecycle of the data and
        set the appropriate variable if possible

        Test
        ----
        unittest implemented
        """
        if self._is_monthly():
            self.time_cycle = 12
        else:
            self._log_warning(
                'WARNING: timecycle can not be set automatically!')

    def _flipud(self):
        """
        flip dataset up down
        """
        if self.data.ndim == 3:
            self.data = self.data[:, ::-1, :]
        elif self.data.ndim == 2:
            self.data = self.data[:: -1, :]
        else:
            raise ValueError('Unsupported geometry for _flipud()')
        if hasattr(self, 'cell_area'):
            if self.cell_area is not None:
                self.cell_area = self.cell_area[::-1, :]
        if hasattr(self, 'lat'):
            if self.lat is not None:
                self.lat = self.lat[::-1, :]

    def _is_sorted(self):
        """
        Checks if timeseries is increasingly sorted

        Test
        ----
        unittest implemented
        """
        return np.all(np.diff(self.time) >= 0.)

    def temporal_smooth(self, N, return_object=True, frac=1.):
        """
        Temporal smoothing of datasets. The routine applies a fast
        approach based on convolution to calculate the temporal smoothing
        of the timeseries.

        Note that the routine does not take into account any time units,
        nor does it check if the data is without gaps. It simply applies
        the smoothing window on the whole timeseries. It is subject
        of the user to ensure that the data is prepared in a way that
        the results make sense.

        Parameters
        ----------
        N : int
            window size for smoothing (needs to be an odd number)
        return_object : bool
            True: return Data object
            False: return numpy array
        frac : float
            minimum fraction of valid timesteps required for calculation

        Test
        ----
        unittest implemented
        """

        # http://stackoverflow.com/questions/13728392/moving-average-or-running-mean
        def _runningMeanFast(x, N):
            N = float(N)
            return np.convolve(x, np.ones((N,)) / N)[
                (N - 1.):]  # returns average. The result is always on the first item!

        N = int(N)
        if N < 3:
            raise ValueError('window need to be >= 3')
        if N % 2 != 1:
            raise ValueError('window size needs to be an odd number!')

        msk = self.get_valid_mask(frac=frac)

        if self.data.ndim == 1:  # single timeseries
            tmp = np.ones_like(self.data) * np.nan
            r = _runningMeanFast(self.data.flatten(), N)
            tmp[N / 2:len(tmp) - N / 2] = r[:-(N - 1)]  # center data!
        elif self.data.ndim == 2:
            raise ValueError(
                'Invalid data geometry for temporal smoothing! (2D)')
        elif self.data.ndim == 3:
            r = np.ones_like(self.data) * np.nan
            nt, ny, nx = self.data.shape
            for i in xrange(ny):  # todo: more efficient implementation
                for j in xrange(nx):
                    if msk[i, j]:
                        r[:, i, j] = _runningMeanFast(self.data[:, i, j], N)
                # ensure that smoothed data is centred in time!
            tmp = np.ones_like(self.data) * np.nan
            # center data!
            tmp[N / 2:len(tmp) - N / 2, :, :] = r[:-(N - 1), :, :]

        # results
        tmp = np.ma.array(tmp, mask=np.isnan(tmp))
        if return_object:
            res = self.copy()
            res.data = tmp
            return res
        else:
            return tmp

    def distance(self, lon_deg, lat_deg, earth_radius=6371.):
        """
        calculate distance of all grid points to a given coordinate
        Note, that calculations are only approximate, as earth is approximated
        as sphere!

        Parameters
        ----------
        lon : float
            longitude [deg]
        lat : float
            latitude [deg]
        earth_radius : float
            earth radius [km]

        Returns
        -------
        returns distance [m]
        """
        from pycmbs.grid import Grid
        assert hasattr(self, 'lat')
        assert hasattr(self, 'lon')
        if not isinstance(self.lat, np.ndarray):
            raise ValueError('Numpy array required!')
        if not isinstance(self.lon, np.ndarray):
            raise ValueError('Numpy array required!')
        G = Grid(np.deg2rad(self.lat), np.deg2rad(self.lon),
                 sphere_radius=earth_radius * 1000.)
        d = G.orthodrome(np.deg2rad(self.lon), np.deg2rad(self.lat),
                         np.deg2rad(lon_deg), np.deg2rad(lat_deg))
        return d

    def _get_center_position(self):
        """
        returns indices of center position in data array
        in case of odd array sizes, the actual center position is
        returned. Otherwise (equal numbers), the center position - 1 is
        returned. Now interpolation is performed. Thus it is always
        ensured that some real data is returned which is close to the
        center

        Returns
        -------
        indices of center position [i,j]
        """

        if (self.ny % 2) == 0:  # equal numbers
            ipos = (self.ny - 1) / 2
        else:  # odd numbers
            ipos = (self.ny - 1) / 2

        if (self.nx % 2) == 0:
            jpos = (self.nx - 1) / 2
        else:
            jpos = (self.nx - 1) / 2

        return ipos, jpos

    def get_center_data(self, return_object=False, flatten=False):
        """
        returns data for center position

        Parameters
        ----------
        return_object : bool
            return the results as a Data object
        flatten : bool
            if True, then the resulting array will be flatteneds
        """
        i, j = self._get_center_position()
        if i is None:
            return None
        if j is None:
            return None

        if self.ndim == 2:
            res = self.data[i, j]
        elif self.ndim == 3:
            res = self.data[:, i, j]
        else:
            assert False

        if return_object:
            r = self.copy()
            if self.ndim == 2:
                res = np.asarray([[res]])
            elif self.ndim == 3:
                res = res.reshape((len(res), 1, 1))
            else:
                assert False
            if flatten:
                res = res.flatten()
            r.data = res
            r.cell_area = np.ones((1, 1))
            return r
        else:
            if flatten:
                return res.flatten()
            else:
                return res

    def _init_sample_object(self, nt=None, ny=20, nx=10):
        """
        initialize the current object as a samle object
        this is in particular usefull for testing

        use this e.g. as
        x = Data(None, None)
        x._init_sample_object(nt=100, ny=500, nx=200)

        Parameters
        ----------
        nt : int
            number of timesteps
        ny : int
            number of rows
        nx : int
            number of cols
        """

        if nt is None:
            data = np.random.random((ny, nx))
        else:
            if ny is None:
                if nx is not None:
                    raise ValueError(
                        'When only timeseries is provided, then nx and ny need to be None!')
                    data
                else:
                    data = np.random.random(nt)
            else:
                data = np.random.random((nt, ny, nx))

        self.data = np.ma.array(data, mask=data != data)
        self.verbose = True
        self.unit = 'myunit'
        self.label = 'testlabel'
        self.filename = 'testinputfilename.nc'
        self.varname = 'testvarname'
        self.long_name = 'This is the longname'
        if nt is not None:
            self.time = np.arange(nt) + plt.datestr2num('2001-01-01')
        self.time_str = "days since 0001-01-01 00:00:00"
        self.calendar = 'gregorian'
        self.oldtime = False
        if ny is None:
            self.cell_area = 1.
            self.lon = None
            self.lat = None
        else:
            self.cell_area = np.ones((ny, nx))
            lat = np.linspace(-90., 90., ny)
            lon = np.linspace(-180., 180., nx)
            self.lon, self.lat = np.meshgrid(lon, lat)

    def _rasterize(self, lon, lat, radius=None, return_object=True):
        """
        rasterize data to a target grid specified by the input arguments

        CAUTION: this is a rather slow function!s

        Parameters
        ----------
        lat : ndarray
            latitude [deg]
        lon : ndarray
            longitude [deg]
        radius : float
            threshold radius
        return_object : bool
            return a Data object

        Returns
        -------
        returns data object with gridded results
        """

        if lon.shape != lat.shape:
            raise ValueError('Inconsistent geometry!')
        if radius is None:
            raise ValueError('Search radius obligatory')

        # flatten data
        dlon = self.lon.flatten()
        dlat = self.lat.flatten()
        data = self.data.flatten()
        res = np.ones_like(lon) * np.nan

        for i in xrange(len(dlon)):
            # distance
            d = np.sqrt((lon - dlon[i]) ** 2. + (lat - dlat[i]) ** 2.)
            #d = np.ma.array(d, mask= d <= radius)
            dmin = d.min()
            m = d == dmin
            if dmin <= radius:  # threshold
                res[m] = data[i]

        if return_object:
            x = self.copy()
            x.lon = lon * 1.
            x.lat = lat * 1.
            x.data = np.ma.array(res, mask=np.isnan(res))
        else:
            raise ValueError('Not implemented yet!')

        return x

    def mask_region(self, r, return_object=False, method='full', maskfile=None, force=False):
        """
        Given a Region object, mask all the data which is outside of the region

        Parameters
        ----------
        r : Region
            Region which contains valid data
        return_object : bool
            if True then a new data object is returned with the masked data
            otherwise the original data is modified
        method : str
            ['full','fast'] two methods for rasterization are supported.
            full: calculate full raster = default
            fast: some faster method, but this might not be totally correct
        maskfile : str
            filename of maskfile
            if provided, then the generated mask is stored in a file specified
            by maskfile. In case that this file is already existing, no raster
            will be generated, but the mask will be read from file. Only exception is if
            force=True
            The filename needs to have the '.nc' extension!
        force : bool
            force always calculation of mask based on polygon information
        """

        print 'Masking by region ...'

        if maskfile is not None:
            if maskfile[-3:] != '.nc':
                print maskfile
                raise ValueError('Maskfile needs to eb a netcdf file!')
            if os.path.exists(maskfile):
                f_rasterize = False
                if force:
                    os.remove(maskfile)
                    f_rasterize = True
            else:
                f_rasterize = True
        else:
            f_rasterize = True

        # perform rasterization
        if f_rasterize:
            polylist = []
            polylist.append(pycmbsPolygon(r.id, zip(r.lon, r.lat)))

            print '   ... rasterizing'
            M = Raster(self.lon, self.lat)
            M.rasterize_polygons(polylist, method=method)
            themask = M.mask
        else:
            # load mask from file!
            MD = Data(maskfile, 'mask', read=True)
            assert MD.ny == self.ny
            assert MD.nx == self.nx
            themask = MD.data

        # save maskfile
        if maskfile is not None:
            if not os.path.exists(maskfile):
                MD = Data(None, None)
                MD._init_sample_object(ny=self.ny, nx=self.nx)
                MD.lon = self.lon
                MD.lat = self.lat
                MD.data = np.ma.array(M.mask, mask=M.mask != M.mask)
                MD.save(maskfile, varname='mask', delete=True)

        if return_object:
            x = self.copy()
        else:
            x = self
        x._apply_mask(themask > 0.)

        if return_object:
            return x
        else:
            return None

    def lomb_scargle_periodogram(self, P, return_object=True, frac=1., corr=True):
        """
        Calculate LOMB-SCARGLE periodogram
        This routine provides a wrapper to the function
        in statistic.py

        Parameters
        ----------
        P : ndarray
            periods [days]
        return_object : bool
            if True -> return Data objects, otherwise
            return ndarrays
        frac : float
            minimum fraction of valid data needed for timesteps to perform calculation
            This is done also for performance improvement!
        """
        from pycmbs.statistic import lomb_scargle_periodogram

        if self.ndim != 3:
            raise ValueError('Only 3D geometry supported!')

        n = len(P)
        A = np.ones((n, self.ny, self.nx)) * np.nan
        B = np.ones((n, self.ny, self.nx)) * np.nan
        if corr:
            R = np.ones((n, self.ny, self.nx)) * np.nan
            PV = np.ones((n, self.ny, self.nx)) * np.nan

        t = self.time
        if 'days since' not in self.time_str:
            raise ValueError('only time units in days currently supported!')

        # get mask where at least
        vmask = self.get_valid_mask(frac=frac)

        for i in xrange(self.ny):
            print i, self.ny
            for j in xrange(self.nx):
                if vmask[i, j]:
                    if corr:
                        A[:, i, j], B[:, i, j], R[:, i, j], PV[:, i, j] = lomb_scargle_periodogram(
                            t, P, self.data[:, i, j], corr=corr)
                    else:
                        A[:, i, j], B[:, i, j] = lomb_scargle_periodogram(
                            t, P, self.data[:, i, j], corr=corr)

        if return_object:
            Aout = Data(None, None)
            Aout._init_sample_object(nt=n, ny=self.ny, nx=self.nx)
            Aout.data = np.ma.array(A, mask=A != A)
            Aout.label = 'amplitude'
            Aout.unit = self.unit

            Bout = Data(None, None)
            Bout._init_sample_object(nt=n, ny=self.ny, nx=self.nx)
            Bout.data = np.ma.array(B, mask=B != B)
            Bout.label = 'phase'
            Bout.unit = 'days'

            if corr:
                Rout = Data(None, None)
                Rout._init_sample_object(nt=n, ny=self.ny, nx=self.nx)
                Rout.data = np.ma.array(R, mask=R != R)
                Rout.label = 'correlation'
                Rout.unit = '-'

                Pout = Data(None, None)
                Pout._init_sample_object(nt=n, ny=self.ny, nx=self.nx)
                Pout.data = np.ma.array(PV, mask=PV != PV)
                Pout.label = 'p-value'
                Pout.unit = '-'

                return Aout, Bout, Rout, Pout
            else:
                return Aout, Bout
        else:
            if corr:
                return A, B, R, PV
            else:
                return A, B

    def get_area(self, valid=True, frac=1.):
        """
        calculate area

        Parameters
        ----------
        valid : bool
            return area for valid pixels only; if false, the area from
            all grid cells is returned
        frac : float
            in case that temporal varying data is used, this parameter
            can define the fraction of timesteps that need to be
            valid to define a cell as valid
        """
        assert hasattr(self, 'cell_area')

        if valid:
            return self.cell_area[self.get_valid_mask(frac=frac)].sum()
        else:
            assert type(
                self.cell_area) == np.ndarray, 'Only numpy arrays for cell_area supported at the moment for this function'
            return self.cell_area.sum()
