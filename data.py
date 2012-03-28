#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

import os,sys

import Nio

import numpy as np

from matplotlib import pylab as plt

from statistic import get_significance



class Data():
    '''
    generic data handling class for pyCMBS
    '''
    def __init__(self,filename,varname,lat_name=None,lon_name=None,read=False,scale_factor = 1.,label=None,unit=None,shift_lon=False,start_time=None,stop_time=None,mask=None,time_cycle=None):
        self.filename     = filename
        self.varname      = varname
        self.scale_factor = scale_factor
        self.lat_name     = lat_name
        self.lon_name     = lon_name
        
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
        
            
    def _shift_lon(self):
        #~ for i in range(len(self.data)):
            #~ d,l = shiftgrid(180.,self.data[i,:,:],self.lon,start=False)
            #~ #self.data[i,:,:] = d[:,:]
        self.lon[self.lon>=180.] = self.lon[self.lon>=180.]-360.
        self._lon360 = False
        
    def read(self,shift_lon,start_time=None,stop_time=None):
        '''
        read data from file
        '''
        if not os.path.exists(self.filename):
            sys.exit('Error: file not existing: '+ self.filename)
            
            
        #read data
        self.data = self.read_netcdf(self.varname)
        #this scaling is related to unit conversion and NOT
        #due to data compression
        self.data = self.data * self.scale_factor 
        
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
            #~ print 'found indices: ', m1,m2
            self._temporal_subsetting(m1,m2)
            #~ print 'After temporal subsetting', np.shape(self.data), m1,m2
        
        
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

        return clim
        
        
    def get_deseasonalized_anomaly(self,base=None):
        '''
        calculate deseasonalized anomalies
        
        base:  all --> use WHOLE original dataset as a reference
            current --> use current dataset as a reference
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
        res = self.copy()
        res.data = ret.copy()
        res.label = self.label + ' anomaly'

        return res
        

        
        
    def set_time(self):
        '''
        set time
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
            self._set_date(basedate,unit='day')
        else:
            print self.filename
            print self.time_str
            sys.exit('ERROR: time type could not be determined!')
        
        
    def _temporal_subsetting(self,i1,i2):
        '''
        perform temporal subsetting
        '''
        
        if i2<i1:
            sys.exit('Invalid indices _temporal_subsetting')
        
        self.time = self.time[i1:i2]
        if self.data.ndim == 3:
            self.data = self.data[i1:i2,:,:]
        elif self.data.ndim == 2:
            self.data = self.data[i1:i2,:]
        else:
            sys.exit('Error temporal subsetting')
            
        
    def _get_time_indices(self,start,stop):
        '''
        determine time indices start/stop
        
        start/stop needs to be a datetime object
        
        returns start/stop index
        '''
        
        #- no subsetting
        if (start == None or stop == None):
            return 0, len(self.time)
        
        
        if stop < start:
            sys.exit('Error: startdate > stopdate')
        
        print start,stop
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
        
        
    def _get_years(self):
        '''
        get years from timestamp
        '''
        
        d = plt.num2date(self.time)
        years = []
        
        for x in d:
            years.append(x.year)
            
        return years
        
        
        
    def _mesh_lat_lon(self):
        if (plt.isvector(self.lat)) & (plt.isvector(self.lon)):
            LON,LAT = np.meshgrid(self.lon,self.lat)
            self.lon = LON; self.lat=LAT
        
        
    def read_netcdf(self,varname):
        '''
        read data from netCDF file
        '''
        F=Nio.open_file(self.filename)
        #print F
        print 'Reading file ', self.filename
        if not varname in F.variables.keys():
            print 'Error: data can not be read. Variable not existing! ', varname
            print F
            F.close()
            return None
        #print F
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
        print '    FillValue: ', self.fill_value
        
        
        #scale factor
        if hasattr(var,'scale_factor'):
            scal = var.scale_factor
        else:
            scal = 1.
            #todo add_offset
        
        print '    scale_factor: ', scal
        data = data * scal
        
        
        print F
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
        
    def timmean(self):
        if self.data.ndim == 3:
            return self.data.mean(axis=0)
        if self.data.ndim == 2:
            #no temporal averaging
            return self.data.copy()
        else:
            sys.exit('Temporal mean can not be calculated as dimensions do not match!')
    
        
        
    def timsum(self):
        if self.data.ndim != 3:
            sys.exit('Temporal sum can not be calculated as dimensions do not match!')
    
        return self.data.sum(axis=0)        
        
    def fldmean(self):
        return np.reshape(self.data,(len(self.data),-1)).mean(axis=1)
        
        
    def _get_label(self):
        if self.unit == None:
            u = ''
        else:
            u = '[' + self.unit + ']'
        return self.label + ' ' + u
        
    def _convert_time(self):
        '''
        convert time that was given as YYYYMMDD.f
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
        
    def _set_date(self,basedate,unit='hour'):
        '''
        basedate = datestring that can be interpreted by datestr2num()
        unit = 'hours', specifies in which units the original data us
        '''
        
        if unit == 'hour':
            scal = 24.
        elif unit == 'day':
            scal = 1.
        else:
            raise ValueError, 'Unsupported unit value'
        
        self.time = plt.datestr2num(basedate) + self.time/scal
        
    def get_aoi(self,region):
        '''
        region of class Region
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
        
    def get_aoi_lat_lon(self,R):
        '''
        get aoi given lat/lon coordinates
        
        the routine masks all area which
        is NOT in the given area
        
        coordinates of region are assumed to be in -180 < lon < 180
        
        given a region R
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
        
        if self.data.ndim == 2:
            tmp = self.data.copy()
            tmp[~msk] = np.nan
            if keep_mask:
                tmp[self.__oldmask] = np.nan
            self.data = np.ma.array(tmp,mask=np.isnan(tmp))
                
        elif self.data.ndim == 3:
            for i in range(len(self.data)):
                tmp = self.data[i,:,:].copy()
                tmp[~msk] = np.nan
                self.data[i,:,:] = tmp[:,:]
                del tmp
                
            if keep_mask:
                self.data.data[self.__oldmask] = np.nan
                
            self.data = np.ma.array(self.data.data,mask=np.isnan(self.data.data))
                
            #self._climatology_raw = np.ma.array(self._climatology_raw,mask=np.isnan(self._climatology_raw))
        else:
            sys.exit('unsupported geometry _apply_mask')
            
            
        if hasattr(self,'_climatology_raw'):
            for i in range(len(self._climatology_raw)):
                tmp = self._climatology_raw[i,:,:].copy()
                tmp[~msk] = np.nan
                self._climatology_raw[i,:,:] = tmp[:,:]
                del tmp
            

        

    def shift_x(self,nx):
        '''
        shift data array in x direction by nx steps
        '''
        
        self.data = self.__shift3D(self.data,nx)
        self.lat  = self.__shift2D(self.lat,nx)
        self.lon  = self.__shift2D(self.lon,nx)

        
    def __shift3D(self,x,n):
        tmp = x.copy()
        y=x.copy()
        y[:,:,:]=nan
        y[:,:,0:n] = tmp[:,:,-n:]
        y[:,:,n:]  = tmp[:,:,0:-n]
        
        return y
        
    def __shift2D(self,x,n):
        tmp = x.copy()
        y=x.copy()
        y[:,:]=nan
        y[:,0:n] = tmp[:,-n:]
        y[:,n:]  = tmp[:,0:-n]
        
        return y
        
        
        
        
        
    def copy(self):
        '''
        copy all attributes
        of self to a new Data object
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
        
        
    def sub(self,x,copy=True):
        '''
        
        @param x: A C{Data} object which will be substracted
        @type  x: Data object
        
        Substract variable x from current field
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
        
    def _sub_sample(self,step):
        '''
        subsample data
        '''
        self.data = self.data[:,::step,::step]
        self.lat  = self.lat [::step,::step]
        self.lon  = self.lon [::step,::step]
        
        
        
    def corr_single(self,x,pthres=1.01,mask=None):
        '''
        The routine correlates a data vector with
        all data of the current object.
        
        @param x: the data vector correlations should be calculated with
        @type  x: numpy array [time]
        
        For efficiency reasons, the calculations are
        performed rowwise for all grid cells using
        corrcoef()
        
        Output: correlation coefficient for each grid point
        
        @todo significance of correlation
        
        '''

        if self.data.ndim != 3:
            raise ValueError, 'Invalid geometry!'
        
        nt,ny,nx = sz = np.shape(self.data)
        
        if nt != len(x):
            raise ValueError, 'Inconsistent geometries'

        R=np.ones((ny,nx))*np.nan #output matrix for correlation
        S=np.ones((ny,nx))*np.nan #output matrix for slope
        I=np.ones((ny,nx))*np.nan #output matrix for intercept
        
        print 'Calculating correlation ...'
        for i in range(ny): #todo how to further increase efficiency?
            if i % 25 == 0:
                print i, '/', ny
            c = np.vstack((self.data[:,i,:].T,x))
            r=np.corrcoef(c)
            
            poly = np.polyfit(x,self.data[:,i,:],1  ) #self.data is dependent variable
            
            R[i,:] = r[0:-1,-1]
            S[i,:] = poly[0,:]
            I[i,:] = poly[1,:]
        
        
        #--- calculate significance (assuming no gaps in the timeseries)
        p_value = get_significance(R,len(x))
        
        
        #--- prepare output data objects
        Rout = self.copy() #copy obkect to get coordinates
        Rout.label = 'correlation'
        Rout.data = np.ma.array(R,mask=p_value>pthres).copy()
        
        Sout = self.copy() #copy object to get coordinates
        Sout.label = 'slope'
        Sout.data = np.ma.array(S,mask=p_value>pthres).copy()
        
        Iout = self.copy() #copy object to get coordinates
        Iout.label = 'intercept'
        Iout.data = np.ma.array(I,mask=p_value>pthres).copy()
        
        Pout = self.copy() #copy object to get coordinates
        Pout.label = 'p-value'
        Pout.data = np.ma.array(p_value).copy()
        
        if mask != None:
            #apply a mask
            Rout._apply_mask(mask)
            Sout._apply_mask(mask)
            Iout._apply_mask(mask)
            Pout._apply_mask(mask)

            
            
        return Rout,Sout,Iout,Pout
        
        
        
        
