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


from cdo import Cdo
import datetime
import pytz
#~ import pickle
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


        calc_cell_area : bool
                calculate cell area
        geometry_file : str
            name of individual file with coordinates. This can be usefull in the case
            that no coordinate information is stored within the actual data files
        """
        self.lat = None
        self.lon = None

        super(Data, self).__init__(filename, varname, **kwargs)

        self.detrended = False



        # assume that coordinates are always in 0 < lon < 360
        self._lon360 = True


        self.level = kwargs.pop('level', None)
        self.gridtype = None

        # specifies if latitudes have been checked for increasing order
        #(required for zonal plot)
        #~ self._latitudecheckok = False
        #~ self.warnings = kwargs.pop('warnings', True)
#~

#~
        #~ if self.geometry_file is not None:
            #~ assert os.path.exists(
                #~ self.geometry_file), 'ERROR: geometry filename provided, but file not existing! ' + self.geometry_file

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
        from geoval.statistic import lomb_scargle_periodogram

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
