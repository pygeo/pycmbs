#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"
__email__ = "alexander.loew@zmaw.de"

import os,sys

import Nio

import numpy as np

from matplotlib import pylab as plt

from statistic import get_significance

import matplotlib.pylab as pl

'''
@todo: data access via opendap
'''

class Data():
    '''
    Data class: main class
    '''
    def __init__(self,filename,varname,lat_name=None,lon_name=None,read=False,scale_factor = 1.,label=None,unit=None,shift_lon=False,start_time=None,stop_time=None,mask=None,time_cycle=None,squeeze=False):
        '''
        Constructor for Data class

        @param filename: name of the file that contains the data  (specify None if not data from file)
        @type filename: str

        @param varname: name of the variable that contains the data (specify None if not data from file)
        @type varname: str

        @param lat_name: name of latitude field
        @type lat_name: str

        @param lon_name: name of longitude field
        @type lon_name: str

        @param read: specifies if the data should be read directly when creating the Data object
        @type read: bool

        @param scale_factor: scale factor to rescale the data after reading. Be aware, that this IS NOT the scale
                             factor that is used for data compression in e.g. netCDF variables (variable attribute scale_factor)
                             but this variable has the purpose to perform a rescaling (e.g. change of units) of the data
                             immediately after it was read (e.g. convert rainfall intensity [mm/h] to [mm/day], the scale_factor would be 1./24.
        @type scale_factor: float

        @param label: the label of the data. This is automatically used e.g. in plotting routines
        @type label: str

        @param unit: specify the unit of the data. Will be used for plotting
        @type unit: str

        @param shift_lon: if this flag is True, then the longitudes are shifted to [-180 ... 180] if they are given
                          in [0 ... 360]
        @type shift_lon: bool

        @param start_time: specifies the start time of the data to be read from file (if None, then beginning of file is used)
        @type start_time: datetime object

        @param stop_time: specifies the stop time of the data to be read from file (if None, then end of file is used)
        @type stop_time: datetime object

        @param mask: mask that is applied to the data when it is read. Needs to have
                     same geometry of data.
        @type mask: array(ny,nx)

        @param time_cycle: specifies size of the periodic cycle of the data. This is needed if
                           deseasonalized anomalies should be calculated. If one works for instance with monthly
                           data, time_cycle needs to be 12. Assumption is that the temporal sampling of
                           the data is regular.
        @type time_cycle: int

        @param squeeze: remove singletone dimensions in the data when it is read
        @type squeeze: bool

        '''

        self.filename     = filename
        self.varname      = varname
        self.scale_factor = scale_factor
        self.lat_name     = lat_name
        self.lon_name     = lon_name
        self.squeeze      = squeeze
        self.squeezed     = False

        self._lon360 = True #assume that coordinates are always in 0 < lon < 360

        self.inmask = mask

        if label == None:
            self.label = self.filename
        else:
            self.label = label

        if unit == None:
            self.unit = None
        else:
            self.unit = unit

        if time_cycle != None:
            self.time_cycle = time_cycle

        if read:
            self.read(shift_lon,start_time=start_time,stop_time=stop_time)


    def set_sample_data(self,a,b,c):
        '''
        fill data matrix with some sample data
        '''
        self.data = plt.rand(a,b,c)

#-----------------------------------------------------------------------

    def _squeeze(self):
        '''
        remove singletone dimensions in data variable
        '''

        print 'SQUEEZING data ... ', self.data.ndim, self.data.shape

        if self.data.ndim > 2:
            self.data = self.data.squeeze()
            self.squeezed = True

        print 'AFTER SQUEEZING data ... ', self.data.ndim, self.data.shape

#-----------------------------------------------------------------------

    def get_zonal_statistics(self,weights,method='mean'):
        '''
        calculate zonal statistics of the data

        returns zonal statistics

        @param weights: grid cell weights for. If weights ARE NOT None, then the aggregation
                        method is automatically set to I{sum}
        @type weights: array(ny,nx)

        @return: returns an array with zonal statistics

        @todo: check weighting of zonal statistics
        '''

        print '... calculating zonal averages'

        if weights != None:
            method = 'sum' #perform weighted sum in case that weights are provided

        dat = self._get_weighted_data(weights)

        if method == 'mean':
            r = dat.mean(axis=self.data.ndim-1)
        elif method == 'sum':
            r = dat.sum(axis=self.data.ndim-1)
        else:
            raise ValueError, 'Invalid option: ' + method

        return r

#-----------------------------------------------------------------------

    def _get_unit(self):
        '''
        get a nice looking string for units

        @return: string with unit like [unit]
        '''
        if self.unit == None:
            u = ''
        else:
            u = '[' + self.unit + ']'

        return u

#-----------------------------------------------------------------------

    def _get_weighted_data(self,weights):
        '''
        calculate area weighted data

        weights are only applied to VALID data
        thus, the weights are renormalized, thus that
        sum(weights_of_valid_data) = 1
        '''

        dat = self.data.copy()

        if weights != None:
            if dat.ndim == 2:
                weights[dat.mask] = 0. #set weights to zero where the data is masked
                sw = weights.sum(); print 'Sum of weights: ', sw
                weights = weights / sw #normalize thus the sum is one
                dat = dat*weights.data
            elif dat.ndim == 3:
                for i in range(len(dat)): #for all timesteps set mask
                    nweights = weights.copy()
                    nweights[dat[i,:,:].mask] = 0. #set weights to zero where the data is masked
                    sw = nweights.sum(); print 'Sum of weights: ', sw
                    nweights = nweights / sw #normalize thus the sum is one
                    print nweights.sum()
                    dat[i,:,:] = dat[i,:,:]*nweights.data


            else:
                raise ValueError, 'Invalid dimensions: not supported yet'

        return dat

#-----------------------------------------------------------------------

    def _shift_lon(self):
        '''
        shift longitude coordinates. Coordinates given as [0...360] are
        converted to [-180...180]

        changes lon field of Data object and sets variable _lon360
        '''
        self.lon[self.lon>=180.] = self.lon[self.lon>=180.]-360.
        self._lon360 = False

#-----------------------------------------------------------------------

    def read(self,shift_lon,start_time=None,stop_time=None):
        '''
        read data from file

        @param shift_lon: if given, longitudes will be shifted
        @type shift_lon: bool

        @param start_time: start time for reading the data
        @type: start_time: datetime object

        @param stop_time: stop time for reading the data
        @type: stop_time: datetime object

        '''
        if not os.path.exists(self.filename):
            sys.exit('Error: file not existing: '+ self.filename)


        #read data
        self.data = self.read_netcdf(self.varname) #o.k.
        #this scaling is related to unit conversion and NOT
        #due to data compression
        print 'scale_factor : ', self.scale_factor
        self.data = self.data * self.scale_factor

        #--- squeeze data to singletone
        if self.squeeze:
            self._squeeze()

        #--- mask data when desired ---
        if self.inmask != None:
            self._apply_mask(self.inmask)

        #read lat/lon
        if self.lat_name != None:
            self.lat = self.read_netcdf(self.lat_name)
        else:
            self.lat = None
        if self.lon_name != None:
            self.lon = self.read_netcdf(self.lon_name)
            #- shift longitudes such that -180 < lon < 180
            if shift_lon:
                self._shift_lon()
        else:
            self.lon = None

        #- read time
        self.time = self.read_netcdf('time') #returns either None or a masked array
        if hasattr(self.time,'mask'):
            self.time = self.time.data
        else:
            self.time = None

        #- determine time
        self.set_time()

        #- lat lon to 2D matrix
        try:
            self._mesh_lat_lon()
        except:
            print 'No lat/lon mesh was generated!'

        #- calculate climatology from ORIGINAL (full dataset)
        if hasattr(self,'time_cycle'):
            self._climatology_raw = self.get_climatology()

        #- perform temporal subsetting
        if self.time != None:
            #- now perform temporal subsetting
            # BEFORE the conversion to the right time is required!
            m1,m2 = self._get_time_indices(start_time,stop_time)
            self._temporal_subsetting(m1,m2)

#-----------------------------------------------------------------------

    def get_yearmean(self,mask=None):
        '''
        This routine calculate the yearly mean of the data field
        A vector with a mask can be provided for further filtering

        e.g. if all the months from JAN-March are masked as TRUE, the
        result will correspnd to the JFM mean for each year

        @param mask: mask [time]
        @type mask : numpy boolean array

        '''

        if mask == None:
            #if not maks is provided, take everything
            mask = np.ones(len(self.time)).astype('bool')
        else:
            if mask.ndim != 1:
                raise ValueError, 'Mask needs to be 1-D of length of time!'
            if len(mask) != len(self.time):
                raise ValueError, 'Mask needs to be 1-D of length of time!'

        #/// get data
        ye = pl.asarray(self._get_years())
        years = pl.unique(ye)
        dat = self.data

        #/// calculate mean
        res = []
        for y in years:
            hlp = (ye == y) & mask
            if self.data.ndim == 1:
                res.append( dat[hlp].mean() )
            else:
                res.append( dat[hlp,:].mean(axis=0) )
        res = pl.asarray(res)

        return years, res

#-----------------------------------------------------------------------

    def get_yearsum(self,mask=None):
        '''
        This routine calculates the yearly sum of the data field
        A vector with a mask can be provided for further filtering

        e.g. if all the months from JAN-March are masked as TRUE, the
        result will correspnd to the JFM sum for each year

        @param mask: mask [time]
        @type mask : numpy boolean array

        '''

        if mask == None:
            #if not maks is provided, take everything
            mask = np.ones(len(self.time)).astype('bool')
        else:
            if mask.ndim != 1:
                raise ValueError, 'Mask needs to be 1-D of length of time!'
            if len(mask) != len(self.time):
                raise ValueError, 'Mask needs to be 1-D of length of time!'

        #/// get data
        ye = pl.asarray(self._get_years())
        years = pl.unique(ye)
        dat = self.data

        #/// calculate mean
        res = []
        for y in years:
            hlp = (ye == y) & mask
            if self.data.ndim == 1:
                res.append( dat[hlp].sum() )
            else:
                res.append( dat[hlp,:].sum(axis=0) )
        res = pl.asarray(res)

        return years, res

#-----------------------------------------------------------------------

    def correlate(self,Y):
        '''
        correlate present data on a grid cell basis
        with another dataset

        @todo: more efficient implementation needed

        @param y: dataset to corrleate the present one with. The
                  data set of self will be used as X in the caluclation
        @type y: C{Data} object

        @return: returns correlation coefficient and its significance
        @rtype: C{Data} objects

        #test routines tcorr1.py, test_corr.py
        '''

        if Y.data.shape != self.data.shape:
            raise ValueError, 'unequal shapes: correlation not possible!'

        #- generate a mask of all samples that are valid in BOTH datasets
        vmask = self.data.mask | Y.data.mask
        vmask = ~vmask

        #- copy original data
        xv = self.data.copy()
        yv = Y.data.copy()
        sdim = self.data.shape

        #- ... and reshape it
        nt = len(self.data)
        xv.shape=(nt,-1); yv.shape=(nt,-1)
        vmask.shape = (nt,-1)

        #- generate new mask for data
        xv.data[xv.mask] = np.nan
        yv.data[yv.mask] = np.nan
        xv[~vmask] = np.nan; yv[~vmask] = np.nan

        xv = np.ma.array(xv,mask=np.isnan(xv))
        yv = np.ma.array(yv,mask=np.isnan(yv))

        #- number of valid data sets where x and y are valid
        nvalid = vmask.sum(axis=0)

        #- calculate correlation only for grid cells with at least 3 valid samples

        mskvalid = nvalid > 2
        r = np.ones(sum(mskvalid)) * np.nan
        p = np.ones(sum(mskvalid)) * np.nan

        xn = xv[:,mskvalid]; yn=yv[:,mskvalid] #copy data
        nv = nvalid[mskvalid]

        #do correlation calculation; currently using np.ma.corrcoef as this
        #supports masked arrays, while stats.linregress doesn't!
        for i in range(sum(mskvalid)):
            #~ slope, intercept, r_value, p_value, std_err = sci.stats.linregress(xn[:,i],yn[:,i])
            r_value = np.ma.corrcoef(xn[:,i],yn[:,i])[0,1]
            r[i]=r_value

            p[i] = get_significance(r_value,nv[i]) #calculate p-value

        #remap to original geometry
        R = np.ones(xv.shape[1]) * np.nan #matrix for results
        P = np.ones(xv.shape[1]) * np.nan #matrix for results
        R[mskvalid] = r; P[mskvalid] = p
        orgshape = (sdim[1],sdim[2])
        R = R.reshape(orgshape) #generate a map again
        P = P.reshape(orgshape)

        R = np.ma.array(R,mask=np.isnan(R))
        P = np.ma.array(P,mask=np.isnan(P))

        RO = self.copy()
        RO.data = R; RO.label = 'correlation' #: ' + self.label + ' ' + y.label
        RO.unit = ''

        PO = self.copy()
        PO.data = P; PO.label = 'p-value' # + self.label + ' ' + y.label
        PO.unit = ''

        return RO,PO





#-----------------------------------------------------------------------

    def get_temporal_mask(self,v,mtype='monthly'):
        '''
        return a temporal mask

        @param v: list of values to be analyzed
        @type v : list of numerical values

        @param mtype: specifies which mask should be applied (valid values: monthly)
        @type mytpe : str

        Example:
        get_temporal_mask([1,2,3],mtype='monthly')
        will return a mask, where the months of Jan-Mar are set to True
        this can be used e.g. further with the routine get_yearmean()

        '''

        valid_types = ['monthly']
        if mtype in valid_types:
            pass
        else:
            raise ValueError, 'Invalid type for mask generation ' + mtype

        #--- get months
        if mtype == 'monthly':
            vals = pl.asarray(self._get_months())
        else:
            raise ValueError, 'Invalid type for mask generation ' + mtype

        #--- generate mask with all months
        mask = pl.zeros(len(self.time)).astype('bool')

        for m in v:
            hlp = vals == m
            mask[hlp] = True

        return pl.asarray(mask)

#-----------------------------------------------------------------------

    def get_climatology(self):
        '''
        calculate climatological mean for a timeincrement
        specified by self.time_cycle
        '''
        if hasattr(self,'time_cycle'):
            pass
        else:
            raise ValueError, 'Climatology can not be calculated without a valid time_cycle'

        if self.data.ndim > 1:
            clim = np.ones(np.shape(self.data[0:self.time_cycle,:])) * np.nan #output grid
        else:
            clim = np.ones(np.shape(self.data[0:self.time_cycle])) * np.nan #output grid

        if clim.ndim == 1:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle] = self.data[i::self.time_cycle].mean(axis=0)
        elif clim.ndim == 2:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle,:] = self.data[i::self.time_cycle,:].mean(axis=0)
        elif clim.ndim ==3:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle,:,:] = self.data[i::self.time_cycle,:,:].mean(axis=0)
        else:
            raise ValueError, 'Invalid dimension when calculating climatology'

        clim = np.ma.array(clim,mask=np.isnan(clim))
        #todo take also into account the mask from the original data set

        return clim

#-----------------------------------------------------------------------

    def get_deseasonalized_anomaly(self,base=None):
        '''
        calculate deseasonalized anomalies

        @param base: specifies the base to be used for the
                     climatology (all: use the WHOLE original dataset
                     as a reference; current: use current data as a reference)
        @type base: str
        '''

        if base == 'current':
            clim = self.get_climatology()
        elif base == 'all':
            clim = self._climatology_raw
        else:
            raise ValueError, 'Anomalies can not be calculated'

        if hasattr(self,'time_cycle'):
            pass
        else:
            raise ValueError, 'Anomalies can not be calculated without a valid time_cycle'

        ret = np.ones(np.shape(self.data)) * np.nan

        if ret.ndim == 1:
            for i in xrange(self.time_cycle):
                ret[i::self.time_cycle]   = self.data[i::self.time_cycle] - clim[i]
        elif ret.ndim == 2:
            for i in xrange(self.time_cycle):
                ret[i::self.time_cycle,:]   = self.data[i::self.time_cycle,:] - clim[i,:]
        elif ret.ndim ==3:
            for i in xrange(self.time_cycle):
                ret[i::self.time_cycle,:,:] = self.data[i::self.time_cycle,:,:] - clim[i,:,:]
        else:
            raise ValueError, 'Invalid dimension when calculating anomalies'

        ret = np.ma.array(ret,mask=np.isnan(ret))

        #return a data object
        res = self.copy(); res.data = ret.copy()
        res.label = self.label + ' anomaly'

        return res

#-----------------------------------------------------------------------

    def set_time(self):
        '''
        This routines sets the timestamp of the data
        to a python type timestamp. Different formats of
        input time are supported
        '''
        if self.time_str == None:
            print 'WARNING: time type can not be determined!'
        elif self.time_str == 'day as %Y%m%d.%f':
            self._convert_time()
        elif 'hours since' in self.time_str:
            basedate = self.time_str.split('since')[1].lstrip()
            self._set_date(basedate,unit='hour')
        elif 'days since' in self.time_str:
            basedate = self.time_str.split('since')[1].lstrip()
            print 'BASEDATE: ', basedate
            self._set_date(basedate,unit='day')
        else:
            print self.filename
            print self.time_str
            sys.exit('ERROR: time type could not be determined!')

#-----------------------------------------------------------------------

    def _temporal_subsetting(self,i1,i2):
        '''
        perform temporal subsetting of the data based on given
        time indices i1,i2

        @param i1: start time index
        @type i1: int

        @param i2: stop time index
        @type i2: int
        '''

        if i2<i1:
            sys.exit('Invalid indices _temporal_subsetting')

        self.time = self.time[i1:i2]
        if self.data.ndim == 3:
            print 'Temporal subsetting for 3D variable ...'
            self.data = self.data[i1:i2,:,:]
        elif self.data.ndim == 2:
            if self.squeezed: #data has already been squeezed and result was 2D (thus without time!)
                print 'Data was already squeezed: not temporal subsetting is performed'
                pass
            else:
                print 'Temporal subsetting for 2D variable ...', self.squeezed
                self.data = self.data[i1:i2,:]
        else:
            sys.exit('Error temporal subsetting')

#-----------------------------------------------------------------------

    def _get_time_indices(self,start,stop):
        '''
        determine time indices start/stop based on data timestamps
        and desired start/stop dates

        @param start: start time
        @type start: datetime object

        @param stop: stop time
        @type stop: datetime object

        @return: returns start/stop indices (int)
        '''

        #- no subsetting
        if (start == None or stop == None):
            return 0, len(self.time)

        if stop < start:
            sys.exit('Error: startdate > stopdate')

        #print start,stop
        s1 = plt.date2num(start); s2=plt.date2num(stop)

        #- check that time is increasing only
        if any(np.diff(self.time)) < 0.:
            sys.exit('Error _get_time_indices: Time is not increasing')

        #- determine indices
        m1 = abs(self.time - s1).argmin()
        m2 = abs(self.time - s2).argmin()

        if m2 < m1:
            sys.exit('Something went wrong _get_time_indices')

        return m1,m2

#-----------------------------------------------------------------------

    def _get_years(self):
        '''
        get years from timestamp

        @return list of years
        '''

        d = plt.num2date(self.time); years = []

        for x in d:
            years.append(x.year)

        return years

#-----------------------------------------------------------------------

    def _get_months(self):
        '''
        get months from timestamp

        @return: returns a list of months
        '''

        d = plt.num2date(self.time); months = []

        for x in d:
            months.append(x.month)

        return months

#-----------------------------------------------------------------------

    def _mesh_lat_lon(self):
        '''
        In case that the Data object lat/lon is given in vectors
        the coordinates are mapped to a 2D field. This routine
        resets the Data object coordinates
        '''
        if (plt.isvector(self.lat)) & (plt.isvector(self.lon)):
            LON,LAT = np.meshgrid(self.lon,self.lat)
            self.lon = LON; self.lat=LAT

#-----------------------------------------------------------------------

    def read_netcdf(self,varname):
        '''
        read data from netCDF file

        @param varname: name of variable to be read
        @type varname: str
        '''
        F=Nio.open_file(self.filename)
        print 'Reading file ', self.filename
        if not varname in F.variables.keys():
            print 'Error: data can not be read. Variable not existing! ', varname
            F.close()
            return None

        var = F.variables[varname]
        data = var.get_value().astype('float').copy()

        self.fill_value = None
        if hasattr(var,'_FillValue'):
            self.fill_value = float(var._FillValue)
            msk = data == self.fill_value
            data[msk] = np.nan #set to nan, as otherwise problems with masked and scaled data
            data = np.ma.array(data,mask=np.isnan(data))
        else:
            data = np.ma.array(data)

        #scale factor
        if hasattr(var,'scale_factor'):
            scal = var.scale_factor
        else:
            scal = 1.

        if hasattr(var,'add_offset'):
            offset = var.add_offset
        else:
            offset = 0.

        data = data * scal + offset

        if 'time' in F.variables.keys():
            tvar = F.variables['time']
            if hasattr(tvar,'units'):
                self.time_str = tvar.units
            else:
                self.time_str = None
        else:
            self.time = None
            self.time_str = None
        F.close()

        return data

#-----------------------------------------------------------------------

    def timmean(self):
        '''
        calculate temporal mean of data field
        '''
        if self.data.ndim == 3:
            return self.data.mean(axis=0)
        if self.data.ndim == 2:
            #no temporal averaging
            return self.data.copy()
        else:
            sys.exit('Temporal mean can not be calculated as dimensions do not match!')

#-----------------------------------------------------------------------

    def timstd(self):
        '''
        calculate temporal standard deviation of data field
        '''
        if self.data.ndim == 3:
            return self.data.std(axis=0)
        if self.data.ndim == 2:
            #no temporal averaging
            return None
        else:
            sys.exit('Temporal standard deviation can not be calculated as dimensions do not match!')

#-----------------------------------------------------------------------

    def timvar(self):
        '''
        calculate temporal variance of data field
        '''
        if self.data.ndim == 3:
            return self.data.var(axis=0)
        if self.data.ndim == 2:
            #no temporal averaging
            return None
        else:
            sys.exit('Temporal variance can not be calculated as dimensions do not match!')

#-----------------------------------------------------------------------

    def timsum(self):
        '''
        calculate temporal sum of data field
        '''
        if self.data.ndim == 3:
            pass
        elif self.data.ndim == 1:
            pass
        else:
            sys.exit('Temporal sum can not be calculated as dimensions do not match!')

        return self.data.sum(axis=0)

#-----------------------------------------------------------------------

    def fldmean(self,return_data = False):
        '''
        calculate mean of the spatial field

        @param return_data: if True, then a C{Data} object is returned
        @type return_data: bool

        @return: vector of spatial mean array[time]
        '''
        if return_data: #return data object
            tmp = np.reshape(self.data,(len(self.data),-1)).mean(axis=1)
            x = np.zeros((len(tmp),1,1))

            x[:,0,0] = tmp
            r = self.copy()
            r.data = np.ma.array(x.copy(),mask=(x-x > 1.) ) #some dummy mask
            return r
        else: #return numpy array
            return np.reshape(self.data,(len(self.data),-1)).mean(axis=1)

#-----------------------------------------------------------------------

    def _get_label(self):
        '''
        return a nice looking label

        @return label (str)
        '''
        u = self._get_unit()
        return self.label + ' ' + u

#-----------------------------------------------------------------------

    def _convert_time(self):
        '''
        convert time that was given as YYYYMMDD.f
        and set time variable of Data object
        '''
        s = map(str,self.time)
        T=[]
        for t in s:
            y = t[0:4]; m=t[4:6]; d=t[6:8]; h=t[8:]
            h=str(int(float(h) * 24.))
            tn = y + '-' + m + '-' + d + ' ' +  h + ':00'
            T.append(tn)
        T=np.asarray(T)
        self.time = plt.datestr2num(T)

#-----------------------------------------------------------------------

    def _set_date(self,basedate,unit='hour'):
        '''
        set Data object time variable

        @param basedate: basis date used for the data;
        @type basedate: str, datestring that can be interpreted by datestr2num()

        @param unit: specify time unit of the time (hour or day)
        @type unit: str
        '''

        if unit == 'hour':
            scal = 24.
        elif unit == 'day':
            scal = 1.
        else:
            raise ValueError, 'Unsupported unit value'

        self.time = (self.time/scal + plt.datestr2num(basedate) ) #+  plt.datestr2num('0001-01-01 00:00:00') #substract basetime of python

#-----------------------------------------------------------------------

    def get_aoi(self,region):
        '''
        region of class Region

        @todo: documentation
        '''

        #copy self
        d = self.copy()

        d.data = region.get_subset(d.data)

        if hasattr(d,'_climatology_raw'):
            d._climatology_raw = region.get_subset(d._climatology_raw)

        if plt.isvector(d.lat):
            d.lat = d.lat[region.y1:region.y2]
        else:
            d.lat  = region.get_subset(d.lat)

        if plt.isvector(d.lon):
            d.lon = d.lon[region.x1:region.x2]
        else:
            d.lon  = region.get_subset(d.lon)

        d.label = d.label + ' (' + region.label + ')'

        return d

#-----------------------------------------------------------------------

    def get_aoi_lat_lon(self,R):
        '''
        get aoi given lat/lon coordinates

        the routine masks all area which
        is NOT in the given area

        coordinates of region are assumed to be in -180 < lon < 180

        given a region R

        @todo: documentation needed
        '''

        #oldmask = self.data.mask[].copy()


        #problem: beim kopieren wird der alte mask value ungültig

        LON = self.lon.copy()
        if self._lon360:
            tmsk = LON>180.
            LON[tmsk] = LON[tmsk] - 360
            del tmsk

        msk_lat = (self.lat >= R.latmin) & (self.lat <= R.latmax)
        msk_lon = (LON      >= R.lonmin) & (LON      <= R.lonmax)
        if R.mask == None: #additional mask in Region object
            msk_region = np.ones(np.shape(msk_lat)).astype('bool')
        else:
            if np.shape(msk_lat) != np.shape(R.mask):
                print np.shape(msk_lat), np.shape(R.mask)
                raise ValueError, 'Invalid geometries for mask'
            else:
                msk_region = R.mask

        msk = msk_lat & msk_lon & msk_region  # valid area

        self._apply_mask(msk)

#-----------------------------------------------------------------------

    def get_valid_mask(self):
        '''
        calculate a mask which is True, when all timestamps
        of the field are valid

        todo this is still not working properly when the data is
        stored in masked arrays, as the mask is not applied in that case!
        '''

        if self.data.ndim == 2:
            return np.ones(self.data.shape).astype('bool')
        elif self.data.ndim == 3:
            n = len(self.data)
            hlp = self.data.copy()
            if hasattr(hlp,'mask'):
                hlp1 = hlp.data.copy()
                hlp1[hlp.mask] = np.nan
            else:
                hlp1 = hlp.data.copy()
            #~ print np.sum(~np.isnan(hlp1),axis=0), n
            msk = np.sum(~np.isnan(hlp1),axis=0) == n
            return msk
        else:
            raise ValueError, 'unsupported dimension!'

#-----------------------------------------------------------------------

    def get_valid_data(self):
        '''
        this routine calculates from the masked array
        only the valid data and returns it together with its
        coordinate as vector
        '''

        n = len(self.time)

        #- vectorize the data
        lon  = self.lon.reshape(-1)
        lat  = self.lat.reshape(-1)
        data = self.data.reshape(n,-1)

        #- extract only valid (not masked data)
        msk = np.sum(~data.mask,axis=0) == n #identify all ngrid cells where all timesteps are valid
        data = data[:,msk]
        lon  = lon[msk]
        lat  = lat[msk]
        del msk

        msk = np.sum(~np.isnan(data),axis=0) == n

        data = data[:,msk]
        lon  = lon[msk]
        lat  = lat[msk]

        return lon,lat,data

#-----------------------------------------------------------------------

    def _apply_mask(self,msk,keep_mask=True):
        '''
        apply a mask to C{Data}. All data where mask==True
        will be masked. Former data and mask will be stored

        @param msk: mask
        @type msk : numpy boolean array

        @param keep_mask: keep old masked
        @type keep_mask : boolean
        '''
        self.__oldmask = self.data.mask.copy()
        self.__olddata = self.data.data.copy()


        print 'Geometry in masking: ', self.data.ndim, self.data.shape

        if self.data.ndim == 2:
            tmp1 = self.data.copy().astype('float') #convert to float to allow for nan support
            tmp1[~msk] = np.nan

            if keep_mask:
                if self.__oldmask.ndim > 0:
                    tmp1[self.__oldmask] = np.nan
            self.data = np.ma.array(tmp1,mask=np.isnan(tmp1))
            del tmp1

        elif self.data.ndim == 3:
            for i in range(len(self.data)):
                tmp              = self.data[i,:,:].copy()
                tmp[~msk]        = np.nan
                self.data[i,:,:] = tmp[:,:]
                del tmp

            if keep_mask:
                if self.__oldmask.ndim > 0:
                    self.data.data[self.__oldmask] = np.nan

            self.data = np.ma.array(self.data.data,mask=np.isnan(self.data.data))
            #self._climatology_raw = np.ma.array(self._climatology_raw,mask=np.isnan(self._climatology_raw))
        else:
            print np.shape(self.data)
            sys.exit('unsupported geometry _apply_mask')

        if hasattr(self,'_climatology_raw'):
            for i in range(len(self._climatology_raw)):
                tmp = self._climatology_raw[i,:,:].copy()
                tmp[~msk] = np.nan
                self._climatology_raw[i,:,:] = tmp[:,:]
                del tmp

        #.... und hier gehts nicht MEHR!



        #~ if self.varname == 'BfCER4e':
            #~
            #~ for i in range(4):
                #~ pl.figure()
                #~ pl.imshow(self.data[i,:,:])
                #~ pl.title('in maskingxxx ' + str(i))
                #~
            #~ print self.__oldmask.shape
            #~ print self.__oldmask.ndim
            #~ print self.__olddata.shape

            #~ for i in range(4):
                #~ pl.figure()
                #~ pl.imshow(self.__oldmask[i,:,:])
                #~
            #~ stop

#-----------------------------------------------------------------------

    def shift_x(self,nx):
        '''
        shift data array in x direction by nx steps

        @param nx: shift by nx steps
        @type nx: int
        '''

        self.data = self.__shift3D(self.data,nx)
        self.lat  = self.__shift2D(self.lat,nx)
        self.lon  = self.__shift2D(self.lon,nx)

#-----------------------------------------------------------------------

    def __shift3D(self,x,n):
        '''
        shift 3D data

        @param x: data to be shifted
        @type x: array(:,:,:)

        @param n: shifting step
        @type n: int
        '''
        tmp = x.copy(); y=x.copy()
        y[:,:,:]=nan
        y[:,:,0:n] = tmp[:,:,-n:]
        y[:,:,n:]  = tmp[:,:,0:-n]

        return y

#-----------------------------------------------------------------------

    def timeshift(self,n):
        '''
        shift data in time by n-steps
        positive numbers mean, that the
        data is shifted leftwards (thus towards earlier times)

        e.g. timeseries 1980,1981,1982, a shift of n=1 would generate
        a dataset with the data ordered as follows
        [1981,1982,1980]

        @param n: lag to shift data (n>=0)
        @type n: int

        @todo: support for n < 0
        '''

        if n == 0:
            return
        if self.data.ndim != 3:
            raise ValueError, 'array of size [time,ny,nx] is needed for temporal shifting!'

        tmp = self.data.copy()
        self.data[:,:,:]=np.nan
        self.data[:-n:,:,:] = tmp[n:,:,:]
        self.data[-n:,:,:]  = tmp[0:n,:,:]



#-----------------------------------------------------------------------

    def __shift2D(self,x,n):
        '''
        shift 2D data

        @param x: data to be shifted
        @type x: array(:,:)

        @param n: shifting step
        @type n: int
        '''
        tmp = x.copy(); y=x.copy()
        y[:,:]=nan
        y[:,0:n] = tmp[:,-n:]
        y[:,n:]  = tmp[:,0:-n]

        return y

#-----------------------------------------------------------------------

    def copy(self):
        '''
        copy complete C{Data} object including all attributes
        @return C{Data} object
        '''
        d = Data(None,None)


        for attr, value in self.__dict__.iteritems():
            try:
                #-copy (needed for arrays)
                cmd = "d." + attr + " = self." + attr + '.copy()'
                exec(cmd)
            except:
                #-copy
                cmd = "d." + attr + " = self." + attr
                exec(cmd)

        return d

#-----------------------------------------------------------------------

    def sub(self,x,copy=True):
        '''
        Substract a C{Data} object from the current object field

        @param x: C{Data} object which will be substracted
        @type  x: Data object

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        '''

        if np.shape(self.data) != np.shape(x.data):
            print 'Inconsistent geometry (sub): can not calculate!'

        if copy:
            d = self.copy()
        else:
            d = self

        d.data = d.data - x.data
        d.label = self.label + ' - ' + x.label

        return d

#-----------------------------------------------------------------------

    def subc(self,x,copy=True):
        '''
        Substract a constant value from the current object field

        @param x: constant
        @type  x: float

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        '''



        if copy:
            d = self.copy()
        else:
            d = self
        d.data -= x
        return d

#-----------------------------------------------------------------------

    def div(self,x,copy=True):
        '''
        Divide current object field by field of a C{Data} object

        @param x: C{Data} object in the denominator
        @type  x: C{Data} object

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        '''

        if np.shape(self.data) != np.shape(x.data):
            print 'Inconsistent geometry (sub): can not calculate!'

        if copy:
            d = self.copy()
        else:
            d = self

        d.data = d.data / x.data
        d.label = self.label + ' / ' + x.label

        return d

#-----------------------------------------------------------------------

    def _sub_sample(self,step):
        '''
        subsample data of current C{Data} object

        @param step: stepsize for subsampling
        @type step: int
        '''
        self.data = self.data[:,::step,::step]
        self.lat  = self.lat [::step,::step]
        self.lon  = self.lon [::step,::step]

#-----------------------------------------------------------------------

    def corr_single(self,x,pthres=1.01,mask=None):
        '''
        The routine correlates a data vector with
        all data of the current object.

        @param x: the data vector correlations should be calculated with
        @type  x: numpy array [time]

        @param pthres: significance threshold. All values below this
                       threshold will be returned as valid
        @type pthres:  float

        @param mask: mask to flag invalid data
        @type mask: array(:,:)

        For efficiency reasons, the calculations are
        performed row-wise for all grid cells using
        corrcoef()

        @return: correlation coefficient for each grid point

        @todo significance of correlation (is the calculation correct?)
        '''

        if self.data.ndim != 3:
            raise ValueError, 'Invalid geometry!'

        nt,ny,nx = sz = np.shape(self.data)

        if nt != len(x):
            raise ValueError, 'Inconsistent geometries'

        R=np.ones((ny,nx))*np.nan #output matrix for correlation
        S=np.ones((ny,nx))*np.nan #output matrix for slope
        I=np.ones((ny,nx))*np.nan #output matrix for intercept
        CO=np.ones((ny,nx))*np.nan #output matrix for covariance

        print 'Calculating correlation ...'
        for i in range(ny): #todo how to further increase efficiency?
            if i % 25 == 0:
                print i, '/', ny
            c = np.vstack((self.data[:,i,:].T,x))
            r=np.corrcoef(c)
            CO[i,:] = np.cov(c)[0:-1,-1]

            poly = np.polyfit(x,self.data[:,i,:],1  ) #self.data is dependent variable

            R[i,:] = r[0:-1,-1]
            S[i,:] = poly[0,:]
            I[i,:] = poly[1,:]


        #--- calculate significance (assuming no gaps in the timeseries)
        p_value = get_significance(R,len(x))


        #--- prepare output data objects
        Rout = self.copy() #copy obkect to get coordinates
        Rout.label = 'correlation'
        msk = (p_value > pthres) | np.isnan(R)
        Rout.data = np.ma.array(R,mask=msk).copy()

        Sout = self.copy() #copy object to get coordinates
        Sout.label = 'slope'
        Sout.data = np.ma.array(S,mask=msk).copy()

        Iout = self.copy() #copy object to get coordinates
        Iout.label = 'intercept'
        Iout.data = np.ma.array(I,mask=msk).copy()

        Pout = self.copy() #copy object to get coordinates
        Pout.label = 'p-value'
        Pout.data = np.ma.array(p_value).copy()

        Cout = self.copy() #copy object to get coordinates
        Cout.label = 'covariance'
        Cout.data = np.ma.array(CO,mask=msk).copy()

        if mask != None:
            #apply a mask
            Rout._apply_mask(mask)
            Sout._apply_mask(mask)
            Iout._apply_mask(mask)
            Pout._apply_mask(mask)
            Cout._apply_mask(mask)

        return Rout,Sout,Iout,Pout, Cout




