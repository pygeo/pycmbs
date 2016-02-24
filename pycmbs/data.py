# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import os
import sys

#~ from geoval.statistic import get_significance, ttest_ind
from pycmbs.netcdf import NetCDFHandler
from pycmbs.polygon import Raster
from pycmbs.polygon import Polygon as pycmbsPolygon

import numpy as np
from matplotlib import pylab as plt
import matplotlib.pylab as pl
from scipy import stats

from netCDF4 import netcdftime
# define module functions here as they are defined in different ways in different versions of netcdftime
# in newer versions of netcdftime, the date2num and num2date functions are part of utime. It is tried
# here to handle newer and older versions
try:
    xxx = netcdftime.date2num
    old_netcdftime = True
    del xxx
except:
    old_netcdftime = False

try:
    xxx = netcdftime.num2date
    old_netcdftime = True
    del xxx
except:
    old_netcdftime = False

from calendar import monthrange
from cdo import Cdo
import datetime
import pytz
import pickle
import datetime

import tempfile
import struct
import gzip

from geoval.core.data import GeoData


class Data(GeoData):

    """
    Data class: main class
    """

    def __init__(self, filename, varname, **kwargs):
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
        geometry_file : str
            name of individual file with coordinates. This can be usefull in the case
            that no coordinate information is stored within the actual data files
        """

        super(Data, self).__init__(filename, varname, **kwargs)

        read=kwargs.pop('read', False)

        self.weighting_type = kwargs.pop('weighting_type', 'valid')
        self.scale_factor = kwargs.pop('scale_factor', 1.)
        self.lat_name = kwargs.pop('lat_name', None)
        self.lon_name = kwargs.pop('lon_name', None)
        label = kwargs.pop('label', None)
        unit = kwargs.pop('unit', None)

        self.squeeze = kwargs.pop('squeeze', False)
        self.squeezed = False
        self.detrended = False

        self.verbose = kwargs.pop('verbose', False)

        time_cycle=kwargs.pop('time_cycle', None)

        shift_lon=kwargs.pop('shift_lon', False)
        start_time=kwargs.pop('start_time', None)
        stop_time=kwargs.pop('stop_time', None)

        time_var=kwargs.pop('time_var', 'time')
        checklat=kwargs.pop('checklat', True)

        self.geometry_file = kwargs.pop('geometry_file', None)

        # assume that coordinates are always in 0 < lon < 360
        self._lon360 = True
        self._calc_cell_area = kwargs.pop('calc_cell_area', True)

        self.inmask = kwargs.pop('mask', None)
        self.level = kwargs.pop('level', None)
        self.gridtype = None
        self._oldtime = kwargs.pop('oldtime', False)
        if self._oldtime:
            print('WARNING: the option _oldtime is depreciated and will be removed in future versions')
        # specifies if latitudes have been checked for increasing order
        #(required for zonal plot)
        self._latitudecheckok = False
        self.warnings = kwargs.pop('warnings', True)

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

        self.cell_area = kwargs.pop('cell_area', None)  # [m**2]

        self.lat = None
        self.lon = None

        if self.geometry_file is not None:
            assert os.path.exists(
                self.geometry_file), 'ERROR: geometry filename provided, but file not existing! ' + self.geometry_file

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


    def _oldtimeoffset(self):
        """
        return offset to convert to old time
        offset is one day *PLUS ONE* following pylab documentation
        This routine takes care of different time units
        """

        assert False, 'This routine is depreciated (oldtimeoffset)!'

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

        # calculate correlations
        rxy, pxy = self.correlate(Y, pthres=pthres)
        rxz, pxz = self.correlate(Z, pthres=pthres)
        rzy, pzy = ZY.correlate(Y, pthres=pthres)

        # calculate partial correlation coefficients
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
        elif self.time_str == 'day as YYYYMMDDhhmm':
            self._convert_time_YYYYMMDDhhmm()
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

    def copy(self):
        """
        copy complete C{Data} object including all attributes
        """
        d = Data(None, None)
        return self._copy_all_attributes(d)




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
