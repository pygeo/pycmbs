#!/usr/bin/env python
# -*- coding: utf-8 -*-
from jockey.ui import bool

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
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

import os

import sys

import Nio
import numpy as np
from matplotlib import pylab as plt
from statistic import get_significance, ttest_ind
import matplotlib.pylab as pl
from scipy import stats
from pyCDO import *

class Data():
    """
    Data class: main class
    """
    def __init__(self,filename,varname,lat_name=None,lon_name=None,read=False,scale_factor = 1.,label=None,unit=None,shift_lon=False,start_time=None,stop_time=None,mask=None,time_cycle=None,squeeze=False,level=None,verbose=False,cell_area=None,time_var='time'):
        """
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

        @param level: specify level to work with (needed for 4D data)
        @type level: int

        @param verbose: verbose mode for printing
        @type verbose: bool

        @param cell_area: area size [m**2] of each grid cell. These weights are uses as an INITIAL weight and are renormalized in accordance to the valid data only.
        @type cell_area: numpy array

        """

        self.filename     = filename; self.varname      = varname
        self.scale_factor = scale_factor; self.lat_name     = lat_name
        self.lon_name     = lon_name; self.squeeze      = squeeze
        self.squeezed     = False; self.detrended    = False

        self.verbose = verbose
        self._lon360 = True #assume that coordinates are always in 0 < lon < 360

        self.inmask = mask; self.level = level

        if label is None:
            self.label = self.filename
        else:
            self.label = label

        if unit is None:
            self.unit = None
        else:
            self.unit = unit

        if time_cycle is not None:
            self.time_cycle = time_cycle

        self.cell_area = cell_area #[m**2]

        self.lat = None
        self.lon = None

        if read:
            self.read(shift_lon,start_time=start_time,stop_time=stop_time,time_var=time_var)

#-----------------------------------------------------------------------

    def __set_sample_data(self,a,b,c):
        """
        fill data matrix with some sample data
        """
        self.data = plt.rand(a,b,c)

#-----------------------------------------------------------------------

    def save(self,filename,varname=None,format='nc',delete=False):
        """
        saves the data object to a file

        @param filename: filename of output file
        @param varname: name of output variable
        @param delete: delete file if existing without asking
        @param format: output format ['nc','txt']
        """

        print 'Saving object in file ' + filename

        if format == 'nc':
            self._save_netcdf(filename,varname=varname,delete=delete)
        elif format == 'ascii':
            self._save_ascii(filename,varname=varname,delete=delete)
        else:
            raise ValueError, 'This output format is not defined yet!'


#-----------------------------------------------------------------------

    def _save_ascii(self,filename,varname=None,delete=False):
        """
        saves the data object to an ASCII file
        (unittest)

        @param filename: filename of output file
        @param varname: name of output variable; this explicitely overwrites self.varname, which is tried to be used as a first order
        @param delete: delete file if existing without asking
        """

        raise ValueError, 'This is not implemented yet!!!'

        #/// check if output file already there
        if os.path.exists(filename):
            if delete:
                os.remove(filename)
            else:
                raise ValueError, 'File already existing. Please delete manually or use DELETE option: ' + filename

        #/// variable name
        if varname == None:
            if self.varname is None:
                varname = 'var1'
            else:
                varname = self.varname


        F = open(filename,'w')

        F.write(str(len(self.time)) + '\n'   ) #number of timesteps
        for i in xrange(len(self.time)):
            F.write(str(pl.num2date(self.time[i])) + ' , ' + str(self.data[i,:,:]) )

        F.close()


        #write 2D field for each timestep






#-----------------------------------------------------------------------

    def _save_netcdf(self,filename,varname=None,delete=False):
        """
        saves the data object to a netCDF file
        (unittest)

        @param filename: filename of output file
        @param varname: name of output variable; this explicitely overwrites self.varname, which is tried to be used as a first order
        @param delete: delete file if existing without asking
        """

        #/// check if output file already there
        if os.path.exists(filename):
            if delete:
                os.remove(filename)
            else:
                raise ValueError, 'File already existing. Please delete manually or use DELETE option: ' + filename

        #/// variable name
        if varname == None:
            if self.varname is None:
                varname = 'var1'
            else:
                varname = self.varname

        #/// create new file
        F = Nio.open_file(filename,mode='w')

        #/// create dimensions
        if self.data.ndim == 3:
            if self.time == None:
                raise ValueError, 'No time variable existing! Can not write 3D data!'
            nt,ny,nx = self.data.shape
            F.create_dimension('time',nt)
        elif self.data.ndim == 2:
            ny,nx = self.data.shape

        F.create_dimension('ny',ny)
        F.create_dimension('nx',nx)

        #/// create variable
        if self.time != None:
            F.create_variable('time','d',('time',))
            F.variables['time'].units = 'days since 0001-01-01 00:00:00 UTC'

        if self.data.ndim == 3:
            F.create_variable(varname,'d',('time','ny','nx'))
        elif self.data.ndim == 2:
            F.create_variable(varname,'d',('ny','nx'))

        if self.lat != None:
            F.create_variable('lat','d',('ny','nx'))
            F.variables['lat'].units = 'degrees_north'
            F.variables['lat'].axis  = "Y"
            F.variables['lat'].long_name  = "latitude"

        if self.lon != None:
            F.create_variable('lon','d',('ny','nx'))
            F.variables['lon'].units = 'degrees_east'
            F.variables['lon'].axis  = "X"
            F.variables['lon'].long_name  = "longitude"


        #/// write data
        F.variables['time'] .assign_value(self.time)
        F.variables[varname].assign_value(self.data)
        if self.lat != None:
            F.variables['lat'].assign_value(self.lat)
        if self.lon != None:
            F.variables['lon'].assign_value(self.lon)

        F.variables[varname].long_name    = self.long_name
        F.variables[varname].units        = self.unit
        F.variables[varname].scale_factor = 1.
        F.variables[varname].add_offset   = 0.

        #/// global attributes
        F.source = self.filename #store original filename

        #/// close file
        F.close()

#-----------------------------------------------------------------------

    def _get_bounding_box(self):
        """
        estimates bounding box of valid data. It returns the indices
        of the bounding box which frames all valid data

        @return: returns indices for bounding box
        """

        msk = self.get_valid_mask() #gives a 2D mask

        #estimate boundary box indices
        xb = msk.sum(axis=0)
        yb = msk.sum(axis=1)

        j1=None; j2=None
        for i in xrange(len(xb)):
            if (j1 == None) and (xb[i]> 0.):
                j1 = i
            if (j2 == None) and (xb[i]> 0.) and (j1 != i):
                j2 = i

        raise ValueError, 'Not finished yet!!!' #todo

        return yb,xb,j1,j2


#-----------------------------------------------------------------------

    def _squeeze(self):
        '''
        remove singletone dimensions in data variable
        '''

        if self.verbose:
            print 'SQUEEZING data ... ', self.data.ndim, self.data.shape

        if self.data.ndim > 2:
            self.data = self.data.squeeze()
            self.squeezed = True

        if self.verbose:
            print 'AFTER SQUEEZING data ... ', self.data.ndim, self.data.shape

#-----------------------------------------------------------------------

    def _set_cell_area(self):
        """
        set cell area size. If a cell area was already given (either by user or from file)
        nothing will happen. Otherwise it will be tried to calculate cell_area from
        coordinates using the CDO's

        @todo: implement the calculation of cell_area in a pythonic way
        """

        if not self.cell_area == None:
            return

        if (self.lat == None) or (self.lon == None):
            print '        WARNING: cell area can not be calculated (missing coordinates)!'
            return

        #--- calculate cell area from coordinates ---
        cell_file = self.filename[:-3]+'_cell_area.nc'

        if not os.path.exists(cell_file): #calculate grid area using CDO's
            cdo = pyCDO(self.filename,None,None)
            cell_file = cdo.gridarea()

        #--- read cell_area file ---
        F=Nio.open_file(cell_file,'r')
        self.cell_area = F.variables['cell_area'].get_value().astype('float').copy()


#-----------------------------------------------------------------------

    def get_zonal_mean(self,return_object=False):
        """
        calculate zonal mean statistics of the data for each timestep
        returns zonal statistics [time,ny]

        uses area weighting of data

        gives exact same results as function 'zonmean' in cdo's

        @return: returns an array with zonal statistics
        @rtype numpy array

        @param return_object: return Data object
        @type return_object: bool

        @todo: implement check if latitudes in y-axis direction are all the same! Otherwise the routine does not make sense
        """

        if self.cell_area is None:
            print 'WARNING: no cell area given, zonal means are based on equal weighting!'
            w = np.ones(self.data.shape)
        else:
            w = self._get_weighting_matrix()

        #/// weight data
        dat = self.data * w

        #/// calculate zonal mean
        if dat.ndim == 2:
            r = dat.sum(axis=1) / w.sum() #zonal mean
        elif dat.ndim == 3:
            nt,ny,nx = dat.shape
            r = np.ones((nt,ny))*np.nan
            W = np.ones((nt,ny))*np.nan

            for i in xrange(nt):
                r[i] = dat[i,:,:].sum(axis=1) / w[i,:,:].sum(axis=1)
                W[i] = w[i,:,:].sum(axis=1)

            r = np.ma.array(r,mask=W==0.)

        else:
            print dat.shape
            raise ValueError, 'Unsupported geometry'

        if return_object:
            res = self.copy()
            res.label = self.label + ' zonal mean'
            res.data  =  r.T #[lat,time]
            res.lat   = self.lat[:,0] #latitudes as a vector
        else:
            res = r


        return res

#-----------------------------------------------------------------------

    def get_percentile(self,p,return_object = True):

        """
        calculate percentile

        uses:
        scipy.stats.mstats.scoreatpercentile(data, per, limit=(), alphap=0.4, betap=0.4
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mstats.scoreatpercentile.html#scipy.stats.mstats.scoreatpercentile

        @param p: percentile value to obtain, e.g. 0.05 corresponds to 5% percentil
        @type p: float

        @param return_object: specifies of a C{Data} object shall be returned [True] or a numpy array [False]
        @type return_object: bool

        @return returns the percentiles as either C{Data} object or as numpy array
        @rtype C{Data} object or numpy array

        @todo: get mask of pixls with at least a few valid samples, so performance is better!
        """

        if self.data.ndim != 3:
            raise ValueError, 'Percentile calculation only supported for 3D data!'

        nt = len(self.data)
        x = self.data.copy(); x.shape = (nt,-1)

        #--- calculate percentile ---
        res = stats.mstats.scoreatpercentile(x,p*100.)

        #--- reshape data array ---
        res.shape = np.shape(self.data[0,:,:])
        res = np.ma.array(res,mask=np.isnan(res))

        #--- return
        if return_object:
            r = self.copy()
            r.label = self.label + '\n percentile: ' + str(round(p,2))
            r.data = res
            return r
        else:
            return res

#-----------------------------------------------------------------------

    def _get_unit(self):
        """
        get a nice looking string for units
        @return: string with unit like [unit]
        """
        if self.unit == None:
            u = ''
        else:
            u = '[' + self.unit + ']'

        return u


#-----------------------------------------------------------------------

    def _shift_lon(self):
        """
        shift longitude coordinates. Coordinates given as [0...360] are
        converted to [-180...180]

        changes lon field of Data object and sets variable _lon360
        """
        self.lon[self.lon>=180.] = self.lon[self.lon>=180.]-360.
        self._lon360 = False

#-----------------------------------------------------------------------

    def read(self,shift_lon,start_time=None,stop_time=None,time_var='time'):
        """
        read data from file
        @param shift_lon: if given, longitudes will be shifted
        @type shift_lon: bool

        @param start_time: start time for reading the data
        @type: start_time: datetime object

        @param stop_time: stop time for reading the data
        @type: stop_time: datetime object

        @param time_var: name of time variable field
        @type time_var: str

        """
        if not os.path.exists(self.filename):
            sys.exit('Error: file not existing: '+ self.filename)
        else:
            print ''
            print 'FILE: ', self.filename


        #read data
        self.data = self.read_netcdf(self.varname) #o.k.
        #this scaling is related to unit conversion and NOT
        #due to data compression
        if self.verbose:
            print 'scale_factor : ', self.scale_factor

        self.data = self.data * self.scale_factor

        #--- squeeze data to singletone
        if self.squeeze:
            self._squeeze()

        #--- mask data when desired ---
        if self.inmask != None:
            self._apply_mask(self.inmask)

        #read lat/lon
        if self.lat_name == None: #try reading lat/lon using default names
            self.lat_name = 'lat'
        if self.lon_name == None:
            self.lon_name = 'lon'
        if self.lat_name != None:
            self.lat = self.read_netcdf(self.lat_name)
        else:
            self.lat = None
        if not self.lon_name is None:
            self.lon = self.read_netcdf(self.lon_name)
            #- shift longitudes such that -180 < lon < 180
            if shift_lon:
                self._shift_lon()
        else:
            self.lon = None

        #- read time
        if time_var != None:
            self.time = self.read_netcdf(time_var) #returns either None or a masked array
            if hasattr(self.time,'mask'):
                self.time = self.time.data
            else:
                self.time = None
        else:
            self.time = None

        #- determine time
        if self.time != None:
            self.set_time()

        #- lat lon to 2D matrix
        try:
            self._mesh_lat_lon()
        except:
            if self.verbose:
                print '        WARNING: No lat/lon mesh was generated!'

        #- cell_area
        #  check if cell_area is already existing. if not, try to calculate from coordinates
        self._set_cell_area()

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

    def get_yearmean(self,mask=None,return_data=False):
        """
        This routine calculate the yearly mean of the data field
        A vector with a mask can be provided for further filtering

        e.g. if all the months from JAN-March are masked as TRUE, the
        result will correspnd to the JFM mean for each year
        (unittest)

        @param mask: temporal mask [time]
        @type mask : numpy boolean array

        @param return_data: specifies if results should be returned as C{Data} object
        @type return_data: bool
        """

        if mask == None:
            #if not mask is provided, take everything
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
        if self.data.ndim == 1:
            res = np.zeros(len(years))*np.nan
            su  = np.zeros(len(years))*np.nan
        elif self.data.ndim == 3:
            nt,ny,nx = self.data.shape
            res = np.zeros((len(years),ny,nx))
            su  = np.zeros((len(years),ny,nx))
        else:
            raise ValueError, 'Unsupported dimension!'

        for i in xrange(len(years)):
            y = years[i]
            hlp = (ye == y) & mask
            if self.data.ndim == 1:
                res[i] = dat[hlp].mean()
                su [i] = dat[hlp].sum()  #calculate sum also (needed for masking in the end)
            else:
                res[i,:,:] = dat[hlp,:].mean(axis=0)
                su [i,:,:] = dat[hlp].sum(axis=0)  #calculate sum also (needed for masking in the end)

        res = np.ma.array(res,mask= (su == 0.)) #this is still not the best solution, but works

        if return_data:
            #generate data object
            r = self.copy()
            r.data = res
            r.time = pl.datestr2num(np.asarray([str(years[i])+'-01-01' for i in range(len(years))]))
            r.time_cycle = 1
            return r
        else:
            return years, res

#-----------------------------------------------------------------------

    def get_yearsum(self,mask=None,return_data=False):
        """
        This routine calculates the yearly sum of the data field
        A vector with a mask can be provided for further filtering

        e.g. if all the months from JAN-March are masked as TRUE, the
        result will correspnd to the JFM sum for each year

        @param mask: mask [time]
        @type mask : numpy boolean array

        """


        if mask is None:
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
        msk = dat.count(0) == 0

        for i in xrange(len(res)):
            res[i,msk] = np.nan

        res = np.ma.array(res,mask=np.isnan(res)) #mask all data that contained no single valid value!


        #res = pl.asarray(res)

        #return years, res xxxxxxxxxxxxxx



        if return_data:
            #generate data object
            r = self.copy()
            r.data = res
            r.time = pl.datestr2num(np.asarray([str(years[i])+'-01-01' for i in range(len(years))]))
            r.time_cycle = 1
            return r
        else:
            return years, res





#-----------------------------------------------------------------------

    def correlate(self,Y,pthres=1.01):
        """
        correlate present data on a grid cell basis
        with another dataset

        @todo: more efficient implementation needed

        @param Y: dataset to corrleate the present one with. The
                  data set of self will be used as X in the caluclation
        @type Y: C{Data} object

        @param pthres: threshold for masking insignificant pixels
        @type pthres: float

        @return: returns correlation coefficient and its significance
        @rtype: C{Data} objects

        #test routines tcorr1.py, test_corr.py

        @todo: implement faster correlation calculation
        @todo: slope calculation as well ???
        @todo: significance correct ???

        """

        if not Y.data.shape == self.data.shape:
            print Y.data.shape, self.data.shape
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
            if p[i] > pthres:
                r[i] = np.nan

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
        """
        return a temporal mask

        @param v: list of values to be analyzed
        @type v : list of numerical values

        @param mtype: specifies which mask should be applied (valid values: monthly)
        @type mytpe : str

        Example:
        get_temporal_mask([1,2,3],mtype='monthly')
        will return a mask, where the months of Jan-Mar are set to True
        this can be used e.g. further with the routine get_yearmean()

        """

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

    def get_climatology(self,return_object=False,nmin=1):
        """
        calculate climatological mean for a time increment
        specified by self.time_cycle

        Note: one can not assume that the climatology starts from January if you use a time_cycle = 12
        Instead, the climatology simply starts with the value which corresponds to the first value of the data

        @param return_object: specifies if a C{Data} object shall be returned
        @type return_object: bool

        @param nmin: specifies the minimum number of datasets used for climatology; else the result is masked
        @type nmin: bool
        """
        if hasattr(self,'time_cycle'):
            pass
        else:
            raise ValueError, 'Climatology can not be calculated without a valid time_cycle'

        #generate output fields
        if self.data.ndim > 1:
            clim = np.ones(np.shape(self.data[0:self.time_cycle,:])) * np.nan #output grid
            slim = np.ones(np.shape(self.data[0:self.time_cycle,:])) * np.nan #output grid
        else:
            clim = np.ones(np.shape(self.data[0:self.time_cycle])) * np.nan #output grid
            slim = np.ones(np.shape(self.data[0:self.time_cycle])) * np.nan #output grid

        if clim.ndim == 1:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle] = self.data[i::self.time_cycle].mean(axis=0)
                slim[i::self.time_cycle] = self.data[i::self.time_cycle].sum(axis=0)
        elif clim.ndim == 2:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle,:] = self.data[i::self.time_cycle,:].mean(axis=0)
                slim[i::self.time_cycle,:] = self.data[i::self.time_cycle,:].sum(axis=0)
        elif clim.ndim ==3:
            for i in xrange(self.time_cycle):
                clim[i::self.time_cycle,:,:] = self.data[i::self.time_cycle,:,:].mean(axis=0)
                slim[i::self.time_cycle,:,:] = self.data[i::self.time_cycle,:,:].sum(axis=0)
        else:
            raise ValueError, 'Invalid dimension when calculating climatology'

        n = slim / clim; del slim #number of data taken into account for climatology
        clim = np.ma.array(clim,mask=( np.isnan(clim) | (n < nmin) | np.isnan(n)) ); del n

        if return_object:
            r = self.copy()
            r.label = r.label + ' - climatology'
            r.data = clim
            r.time = []
            for i in xrange(self.time_cycle):
                r.time.append(self.time[i])
            r.time = np.asarray(r.time)

            if len(r.time) != len(r.data):
                print len(r.time)
                print len(r.data)
                raise ValueError, 'Data and time are inconsistent in get_climatology()'

            return r
        else:
            return clim




#-----------------------------------------------------------------------

    def get_deseasonalized_anomaly(self,base=None):
        """
        calculate deseasonalized anomalies

        @param base: specifies the base to be used for the
                     climatology (all: use the WHOLE original dataset
                     as a reference; current: use current data as a reference)
        @type base: str

        @return: returns a C{Data} object with the anomalies
        @rtype: C{Data}
        """

        if base == 'current':
            clim = self.get_climatology()
        elif base == 'all':
            clim = self._climatology_raw
        else:
            raise ValueError, 'Anomalies can not be calculated, invalid BASE'

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

        ret = np.ma.array(ret,mask=(np.isnan(ret) | self.data.mask) )

        #--- return a data object
        res = self.copy(); res.data = ret.copy()
        res.label = self.label + ' anomaly'

        return res

#-----------------------------------------------------------------------

    def set_time(self):
        """
        This routines sets the timestamp of the data
        to a python type timestamp. Different formats of
        input time are supported
        """
        if self.time_str is None:
            print '        WARNING: time type can not be determined!'
        elif self.time_str == 'day as %Y%m%d.%f':
            self._convert_time()
        elif 'hours since' in self.time_str:
            basedate = self.time_str.split('since')[1].lstrip()
            self._set_date(basedate,unit='hour')
        elif 'days since' in self.time_str:
            basedate = self.time_str.split('since')[1].lstrip()
            if self.verbose:
                print 'BASEDATE: ', basedate
            self._set_date(basedate,unit='day')
        elif 'months since' in self.time_str:
            basedate = self.time_str.split('since')[1].lstrip()
            print 'BASEDATA: ', basedate
            if self.verbose:
                print 'BASEDATE: ', basedate
            self._set_date(basedate,unit='month')
        else:
            print self.filename
            print self.time_str
            sys.exit("ERROR: time type could not be determined!")

#-----------------------------------------------------------------------

    def _temporal_subsetting(self,i1,i2):
        """
        perform temporal subsetting of the data based on given
        time indices i1,i2

        @param i1: start time index
        @type i1: int

        @param i2: stop time index
        @type i2: int
        """

        if i2<i1:
            sys.exit('Invalid indices _temporal_subsetting')

        self.time = self.time[i1:i2]
        if self.data.ndim == 3:
            if self.verbose:
                print 'Temporal subsetting for 3D variable ...'
            self.data = self.data[i1:i2,:,:]
        elif self.data.ndim == 2:
            if self.squeezed: #data has already been squeezed and result was 2D (thus without time!)
                print 'Data was already squeezed: not temporal subsetting is performed!'
                pass
            else:
                self.data = self.data[i1:i2,:]
        elif self.data.ndim == 1: #single temporal vector assumed
            self.data = self.data[i1:i2]
        else:
            sys.exit('Error temporal subsetting: invalid dimension!')

#-----------------------------------------------------------------------

    def interp_time(self,t,method='linear'):
        """
        interpolate data matrix in time. The existing data is interpolated to a new temporal spaceing that
        is specified by the time vector argument 't'
        The interpolation is done, by constructing a weighting matrix which basically performs a linear
        interpolation as y = w*x(1) + (1-w)*x(2)

        @param t: vector of time where the data should be interpolated to. The vector is expected to correspond
                  to numbers which correspond to the python standard for time (see datestr2num documentation). The array
                  needs to be in ascending order, otherwise an error occurs.
        @type t: numpy array

        @param method: option to specify interpolation method. At the moment, only linear interpolation is supported!
        @type method: str

        @return: returns a new C{Data} object that contains the interpolated values
        @rtype: Data

        @todo: still some boundary effects for last timestep
        """

        if method != 'linear':
            raise ValueError, 'Only linear interpolation supported at the moment so far!'

        #/// checks
        if self.data.ndim !=3:
            raise ValueError, 'Interpolation currently only supported for 3D arrays!'
        if not np.all(np.diff(t) > 0):
            raise ValueError, 'Input time array is not in ascending order! This must not happen! Please ensure ascending order'
        if not np.all(np.diff(self.time) > 0):
            raise ValueError, 'Time array of data is not in ascending order! This must not happen! Please ensure ascending order'

        nt0,ny,nx = self.data.shape #original dimensions
        nt = len(t) #target length of time

        #/// copy data
        X = self.data.copy(); X.shape = (nt0,-1) #[time,npix]
        nt0,npix = X.shape

        #/// preliminary checks
        f_err = False

        #A) all data is BEFORE desired period
        if self.time.max() < t.min():
            print 'WARNING: specified time period is BEFORE any data availability. NO INTERPOLATION CAN BE DONE!'
            f_err = True

        #B) all data is AFTER desired period
        if self.time.min() > t.max():
            print 'WARNING: specified time period is AFTER any data availability. NO INTERPOLATION CAN BE DONE!'
            f_err = True

        if f_err:
            r = self.copy(); r.data = np.zeros(nt,ny,nx)
            return r

        #/// construct weighting matrix
        W = np.zeros((nt,nt0))
        i1 = 0; i2 = 1 #indices in original data
        f_init=True
        for i in xrange(nt-1):
            #1) find start of interpolation period
            if f_init:
                while self.time[i2] <= t[0]: #do nothing while data coverage not reached yet
                    i2+=1
                    continue
            f_init=False
            i1 = i2 - 1
            if i1 < 0:
                raise ValueError, 'Invalid index i1:'


            #/// increment
            if i2 < nt0:
                if self.time[i2] < t[i]:
                    if i2 <= nt0-1:
                        i2 += 1
            if i1 < nt0-1:
                if self.time[i1+1]< t[i]:
                    i1 += 1
            else:
                continue


            #2) check consistency
            if self.time[i1] > t[i]:
                #the first timeperiod with valid data has not been reached yet
                # ... loop
                continue
            if i2 >=nt0:
                continue




            if self.time[i2] < t[i]:
                print self.time[i1], t[i], self.time[i2]
                print i1,i,i2
                raise ValueError, 'interp_time: this should not happen!'


            if i2 > nt0-1:
                break

            #... here we have valid data
            #print i,i1,i2, nt0
            t1   = self.time[i1]; t2=self.time[i2]
            W[i,i1] = (t2 - t[i]) / (t2-t1)
            W[i,i2] = 1.-W[i,i1]

            #... now increment if needed
            if i < (nt0-1):
                if t2 < t[i+1]:
                    i1 +=1; i2+=1

        #/// generate interpolation Matrix and perform interpolation
        N = np.ma.dot(W,X) #could become a problem for really large matrices!
        N[nt-1,:] = np.nan #avoid boundary problem (todo: where is the problem coming from ??)
        #mask all data that is outside of valid time period
        msk = (t < self.time.min()) | (t > self.time.max())
        N[msk,:] = np.nan
        N.shape = (nt,ny,nx)

        res = self.copy()
        res.time = t
        res.data = np.ma.array(N,mask=np.isnan(N))
        del N

        return res

    #-----------------------------------------------------------------------

    def _get_time_indices(self,start,stop):
        """
        determine time indices start/stop based on data timestamps
        and desired start/stop dates

        @param start: start time
        @type start: datetime object

        @param stop: stop time
        @type stop: datetime object

        @return: returns start/stop indices (int)
        @rtype: int
        """

        #- no subsetting
        if (start == None or stop == None):
            return 0, len(self.time)
        if stop < start:
            sys.exit('Error: startdate > stopdate')

        s1 = plt.date2num(start); s2=plt.date2num(stop)

        #- check that time is increasing only
        if any(np.diff(self.time)) < 0.:
            sys.exit('Error _get_time_indices: Time is not increasing')

        #- determine indices
        m1 = abs(self.time - s1).argmin(); m2 = abs(self.time - s2).argmin()

        if self.time[m1] < s1:
            m1 += 1
        if self.time[m2] > s2:
            m2 -= 1
        if m2 < m1:
            sys.exit('Something went wrong _get_time_indices')

        return m1,m2

#-----------------------------------------------------------------------

    def _get_years(self):
        """
        get years from timestamp

        @return list of years
        """
        d = plt.num2date(self.time); years = []
        for x in d:
            years.append(x.year)
        return years

#-----------------------------------------------------------------------

    def _get_months(self):
        """
        get months from timestamp
        @return: returns a list of months
        """
        d = plt.num2date(self.time); months = []
        for x in d:
            months.append(x.month)
        return months

#-----------------------------------------------------------------------

    def _mesh_lat_lon(self):
        """
        In case that the Data object lat/lon is given in vectors
        the coordinates are mapped to a 2D field. This routine
        resets the Data object coordinates
        """
        if (plt.isvector(self.lat)) & (plt.isvector(self.lon)):
            LON,LAT = np.meshgrid(self.lon,self.lat)
            self.lon = LON; self.lat=LAT
        else:
            pass

#-----------------------------------------------------------------------

    def read_netcdf(self,varname):
        """
        read data from netCDF file

        @param varname: name of variable to be read
        @type varname: str
        """
        F=Nio.open_file(self.filename,'r')
        if self.verbose:
            print 'Reading file ', self.filename
        if not varname in F.variables.keys():
            if self.verbose:
                print '        WARNING: data can not be read. Variable not existing! ', varname
            F.close()
            return None

        var = F.variables[varname]
        data = var.get_value().astype('float').copy()

        if data.ndim > 3:
            if self.level == None:
                print data.shape
                raise ValueError, '4-dimensional variables not supported yet! Either remove a dimension or specify a level!'
            else:
                data = data[:,self.level,:,:] #[time,level,ny,nx ] --> [time,ny,nx]

        self.fill_value = None
        if hasattr(var,'_FillValue'):
            self.fill_value = float(var._FillValue)
            msk = data == self.fill_value
            data[msk] = np.nan #set to nan, as otherwise problems with masked and scaled data
            data = np.ma.array(data,mask=np.isnan(data))
        else:
            data = np.ma.array(data,mask=np.zeros(data.shape).astype('bool'))

        #scale factor
        if hasattr(var,'scale_factor'):
            if plt.is_string_like(var.scale_factor):
                scal = float(var.scale_factor.replace('.f','.'))
            else:
                scal = var.scale_factor
        else:
            scal = 1.

        if hasattr(var,'add_offset'):
            offset = float(var.add_offset)
        else:
            offset = 0.

        data = data * scal + offset


        if hasattr(var,'long_name'):
            self.long_name = var.long_name
        else:
            self.long_name = '-'


        #check if file has cell_area attribute and only use it if it has not been set by the user
        if 'cell_area' in F.variables.keys() and self.cell_area == None:
            self.cell_area = F.variables['cell_area'].get_value().astype('float').copy() #unit should be in m**2

    #set units if possible; if given by user, this is taken
        #otherwise unit information from file is used if available
        if self.unit == None and hasattr(var, 'units'):
            self.unit = var.units

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

    def temporal_trend(self,return_object=False, pthres=1.01):
        """
        calculate temporal trend of the data over time
        the slope of the temporal trend
        (unittest)

        @param return_object; specifies if a C{Data} object shall be returned [True]
                              or if a numpy array shall be returned [False]

        @param pthres: specifies significance threshold; all values above this threshold will be masked
        @return: returns either C{Data} object or a numpy array. The following variables are returned: correlation, slope, intercept, p-value
                 the slope which is returned has unit [dataunit/day]
        """
        x = self.time
        R,S,I,P,C = self.corr_single(x,pthres=pthres)

        if return_object:
            S.unit += ' / day'
            return R,S,I,P
        else:
            return R.data,S.data,I.data,P.data

#-----------------------------------------------------------------------

    def timmean(self,return_object=False):
        """
        calculate temporal mean of data field

        @param return_object: specifies if a C{Data} object shall be returned [True]; else a numpy array is returned
        @type return_object: bool
        """
        if self.data.ndim == 3:
            res = self.data.mean(axis=0)
        elif self.data.ndim == 2:
            #no temporal averaging
            res = self.data.copy()
        else:
            print self.data.ndim
            sys.exit('Temporal mean can not be calculated as dimensions do not match!')

        if return_object:
            tmp = self.copy(); tmp.data = res
            return tmp
        else:
            return res

#-----------------------------------------------------------------------

    def timcv(self,return_object=True):
        """
        calculate temporal coefficient of variation

        @param return_object: specifies if a C{Data} object shall be returned [True]; else a numpy array is returned
        @type return_object: bool
        """
        res = self.timstd (return_object=False) / self.timmean(return_object=False)
        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy(); tmp.data = res; tmp.label=self.label + ' (CV)'; tmp.unit='-'
                return tmp
        else:
            return res


#-----------------------------------------------------------------------

    def timstd(self,return_object=False):
        """
        calculate temporal standard deviation of data field

        @param return_object: specifies if a C{Data} object shall be returned [True]; else a numpy array is returned
        @type return_object: bool
        """
        if self.data.ndim == 3:
            res = self.data.std(axis=0)
        elif self.data.ndim == 2:
            #no temporal averaging
            res = None
        else:
            sys.exit('Temporal standard deviation can not be calculated as dimensions do not match!')

        if return_object:
            if res is None:
                return res
            else:
                tmp = self.copy(); tmp.data = res
                return tmp
        else:
            return res

#-----------------------------------------------------------------------

    def timvar(self):
        """
        calculate temporal variance of data field
        """
        if self.data.ndim == 3:
            return self.data.var(axis=0)
        if self.data.ndim == 2:
            #no temporal averaging
            return None
        else:
            sys.exit('Temporal variance can not be calculated as dimensions do not match!')

#-----------------------------------------------------------------------

    def timsum(self):
        """
        calculate temporal sum of data field
        """
        if self.data.ndim == 3:
            pass
        elif self.data.ndim == 1:
            pass
        else:
            sys.exit('Temporal sum can not be calculated as dimensions do not match!')

        return self.data.sum(axis=0)

#-----------------------------------------------------------------------

    def timn(self):
        """
        calculate number of samples
        done via timmean and timsum to take
        into account the valid values only
        """
        return self.timsum() / self.timmean()

#-----------------------------------------------------------------------

    def _get_weighting_matrix(self):
        """
        (unittest)

        get matrix for area weigthing of grid cells. For each timestep
        the weights are calculated as a function  of the number of valid
        grid cells.

        The returned array contains weights for each timestep. The sum
        of these weights is equal to one for each timestep.

        @return weighting matrix in same geometry as original data
        @rtype numpy array
        """

        w = np.zeros(self.data.shape)

        if self.data.ndim == 2:
            m = ~self.data.mask
            w[m] = self.cell_area[m] / self.cell_area[m].sum()
            return w
        elif self.data.ndim == 3:
            nt = len(self.data)

            #1) repeat cell area nt-times
            cell_area = self.cell_area.copy()
            s = np.shape(cell_area)
            cell_area.shape = (-1)
            if len(s) == 2:
                w = cell_area.repeat(nt).reshape((s[0]*s[1],nt)).T
            elif len(s) == 1:
                w = cell_area.repeat(nt).reshape((1,nt)).T
            else:
                raise ValueError, 'Invalid geometry!'
            print self.label
            print w.shape
            print self.data.shape
            w.shape = self.data.shape #geometry is the same now as data

            #2) mask areas that do not contain valid data
            w = np.ma.array(w,mask=self.data.mask)

            #3) calculate for each time the sum of all VALID grid cells --> normalization factor
            no = w.reshape(nt,-1).sum(axis=1) #... has size nt

            #4) itterate over all timesteps and calculate weighting matrix
            for i in xrange(nt):
                w[i,:,:] /= no[i]
            return w
        else:
            raise ValueError, 'weighting matrix not supported for this data shape'


#-----------------------------------------------------------------------

    def mean(self,apply_weights=True):
        """
        calculate mean of the spatial field using weighted averaging

        @param apply_weights: apply weights when calculating area weights
        @type apply_weights: bool
        """
        if apply_weights:
            #area weighting
            w = self._get_weighting_matrix() #get weighting matrix for each timestep (taking care of invalid data)
            w *= self.data #multiply the data with the weighting matrix in memory efficient way
            return w.sum() #... gives weighted sum = mean
        else:
            #no area weighting
            return self.data.mean()

#-----------------------------------------------------------------------

    def fldmean(self,return_data = False,apply_weights=True):
        """
        calculate mean of the spatial field for each time using weighted averaging
        results are exactly the same as one would obtain with the similar
        cdo function

        @param return_data: if True, then a C{Data} object is returned
        @type return_data: bool

        @param apply_weights: apply weights when calculating area weights
        @type apply_weights: bool

        @return: vector of spatial mean array[time]
        @rtype: C{Data} object or numpy array
        """

        if apply_weights:
            #area weighting
            w = self._get_weighting_matrix() #get weighting matrix for each timestep (taking care of invalid data)
            w *= self.data #multiply the data with the weighting matrix in memory efficient way
            w.shape = (len(self.data),-1)
            tmp = w.sum(axis=1) #... gives weighted sum
        else:
            #no area weighting
            tmp = np.reshape(self.data,(len(self.data),-1)).mean(axis=1)

        #////
        if return_data: #return data object
            x = np.zeros((len(tmp),1,1))
            x[:,0,0] = tmp
            assert(isinstance(tmp,np.ma.masked_array))
            r = self.copy()
            r.data = np.ma.array(x.copy(),mask=tmp.mask ) #use mask of array tmp (important if all values are invalid!)

            #return cell area array with same size of data
            r.cell_area = np.array([1.])

            return r
        else: #return numpy array
            return tmp


#-----------------------------------------------------------------------

    def fldstd(self,return_data = False,apply_weights=True):
        """
        calculate stdv of the spatial field using area weighting

        returns exactly same results as the same CDO function

        @param return_data: if True, then a C{Data} object is returned
        @type return_data: bool

        @return: vector of spatial std array[time]
        """

        if apply_weights:
            #calculate weighted standard deviation.
            #http://en.wikipedia.org/wiki/Mean_square_weighted_deviation
            #(adapted from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy)

            #calculate weighting matrix
            w = self._get_weighting_matrix() #get weighting matrix for each timestep (taking care of invalid data)

            #some temp. variables
            h2 = self.data*w  #wx
            h1 = self.data*h2 #w*x**2

            #do calculation
            nt,ny,nx = self.data.shape

            s = np.ones(nt)*np.nan #generate output array
            for i in xrange(nt):
                s[i] = h1[i,:,:].sum() * w[i,:,:].sum() - h2[i,:,:].sum()**2
                s[i] /= (w[i,:,:].sum()**2  - (w[i,:,:]*w[i,:,:]).sum()  )
            tmp = np.sqrt(s)

        else:
            #no area weighting
            tmp = np.reshape(self.data,(len(self.data),-1)).std(axis=1)


        #--- return either a Data object or a numpy array
        if return_data: #return data object
            tmp = np.reshape(self.data,(len(self.data),-1)).std(axis=1)
            x = np.zeros((len(tmp),1,1))

            x[:,0,0] = tmp
            r = self.copy()
            r.data = np.ma.array(x.copy(),mask=tmp.mask )
            return r
        else: #return numpy array
            return tmp


#-----------------------------------------------------------------------

    def _get_label(self):
        """
        return a nice looking label

        @return: label
        @rtype: str
        """

        if hasattr(self,'label'):
            pass
        else:
            self.label = ''

        u = self._get_unit()
        return self.label + ' ' + u

#-----------------------------------------------------------------------

    def _convert_time(self):
        """
        convert time that was given as YYYYMMDD.f
        and set time variable of Data object
        """
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
    def adjust_time(self,day=None,month=None,year=None):
        """
        correct all timestamps and assign
        same day and/or month

        (unittest)

        @param day: day to apply to all timestamps
        @type day: int

        @param month: year to apply to all timestamps
        @type month: int
        """

        o = []
        for t in self.time:
            d = plt.num2date(t)
            s = str(d) #convert to a string
            if day is not None:
                s = s[0:8] + str(day).zfill(2) + s[10:] #replace day
            if month is not None:
                s = s[0:5] + str(month).zfill(2) + s[7:] #replace day
            if year is not None:
                s = str(year).zfill(4) + s[4:]
            o.append(plt.datestr2num(s))

        o = np.asarray(o)
        self.time = o.copy()

#-----------------------------------------------------------------------

    def timsort(self,return_object = False):
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

        (unittest)

        @param return_object: return a C{Data} object as result
        @type return_object: bool

        @return: either a C{Data} object is returned or the current data object is modified
        """

        #- checks
        if self.time == None:
            raise ValueError, 'Time array needed for timsort()'
        if self.data.ndim != 3:
            raise ValueError, '3D array needed for timsort()'

        #- specify object ot work on
        if return_object:
            x = self.copy()
        else:
            x = self

        #- do the sorting
        s = np.argsort(x.time)
        x.data = x.data[s,:,:]
        x.time = x.time[s]
        if hasattr(x,'std'): #standard deviation
            x.std = x.std[s,:,:]
        if hasattr(x,'n'):   #number of datasets
            x.n = x.n[s,:,:]

        #- result
        if return_object:
            return x





#-----------------------------------------------------------------------

    def _set_date(self,basedate,unit='hour'):
        """
        set C{Data} object time variable

        @param basedate: basis date used for the data;
        @type basedate: str, datestring that can be interpreted by datestr2num()

        @param unit: specify time unit of the time (hour or day)
        @type unit: str
        """

        #explicit conversion of datestr, instead of using datestr2num(), as datestr2num can NOT
        #handle appropriate basedates like 0001-01-01 00:00:00, as its interpretation is that this
        #corresponds to a date of 01/01/2001 !!!
        from datetime import datetime

        try: #check if the date is in the format YYYYDDMM HH:MM:SS, this is needed to avoid that string like YMMDD cause an error
            b = basedate.split(' ')
            c = b[0].split('-')
            d = b[1].split(':')
            basedate = c[0].zfill(4)+'-'+c[1].zfill(2)+'-'+c[2].zfill(2)+' '+d[0].zfill(2)+':'+d[1].zfill(2)+':'+d[2].zfill(2)
            bdate = datetime.strptime(basedate,'%Y-%m-%d %H:%M:%S')
        except:
            raise ValueError, 'basedate is formatted in an unexpected way: ' + basedate

        if unit == 'hour':
            scal = 24.
            self.time = (self.time/scal + plt.date2num(bdate) )
        elif unit == 'day':
            scal = 1.
            self.time = (self.time/scal + plt.date2num(bdate) )
            if self.verbose:
                print 'print set_time: time', self.time
                print 'Basedate: ', basedate, bdate
        elif unit == 'month':
            #months since basedate
            from dateutil.rrule import rrule, MONTHLY
            from datetime import datetime
            bdate = plt.num2date(plt.datestr2num(basedate))
            sdate = [d for d in rrule(MONTHLY,dtstart=datetime(bdate.year,bdate.month,bdate.day),count=self.time[0]+1)] #calculate date of first dataset
            sdate=sdate[-1] #last date as start date

            interval = np.diff(self.time)
            msk = interval == interval[0]
            interval = map(int,interval)

            if ~all(msk):
                #--- the months are not equally spaced, therefore generate a list manually
                from datetime import date
                from dateutil.relativedelta import relativedelta
                x=[]
                bdate0 = datetime(bdate.year,bdate.month,bdate.day)
                print bdate0
                for t in self.time:
                    x.append(bdate0 + relativedelta( months = int(t) ))
                self.time = plt.date2num(x)
            else:
                #--- monthly timeseries for equidistant data
                x=[d for d in rrule(MONTHLY,dtstart=datetime(sdate.year,sdate.month,sdate.day),count=len(self.time),interval=interval[0] )] #calculate series of data
                self.time = plt.date2num(x)
        else:
            raise ValueError, 'Unsupported unit value'



#-----------------------------------------------------------------------

    def get_aoi(self,region):
        """
        region of class Region

        @todo: documentation
        """

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

    def get_aoi_lat_lon(self,R,apply_mask=True):
        """
        get aoi given lat/lon coordinates

        the routine masks all area which
        is NOT in the given area

        coordinates of region are assumed to be in -180 < lon < 180

        given a region R

        @todo: documentation needed
        """

        #oldmask = self.data.mask[].copy()


        #problem: beim kopieren wird der alte mask value ungültig

        LON = self.lon.copy()
        if self._lon360:
            tmsk = LON>180.
            LON[tmsk] -= 360.
            del tmsk

        msk_lat = (self.lat >= R.latmin) & (self.lat <= R.latmax)
        msk_lon = (LON      >= R.lonmin) & (LON      <= R.lonmax)
        if (R.mask == None) | (apply_mask==False): #additional mask in Region object
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
        """
        calculate a mask which is True, when all timestamps
        of the field are valid

        todo this is still not working properly when the data is
        stored in masked arrays, as the mask is not applied in that case!
        """

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

    def get_valid_data(self,return_mask=False,mode='all'):
        """
        this routine calculates from the masked array
        only the valid data and returns it together with its
        coordinate as vector

        valid means that ALL timestamps need to be valid!

        @param return_mask: specifies if the mask applied to the original data should be returned as well
        @type return_mask: bool

        @param mode: analysis mode: 'all'=all timestamps need to be valid, 'one'=at least a single dataset needs to be valid
        """

        n = len(self.time)

        #- vectorize the data
        if hasattr(self,'lon'):
            if self.lon != None:
                lon  = self.lon.reshape(-1)
            else:
                lon = None
        else:
            lon = None
        if hasattr(self,'lat'):
            if self.lat != None:
                lat  = self.lat.reshape(-1)
            else:
                lat = None
        else:
            lat = None

        data = self.data.reshape(n,-1)
        data.mask[np.isnan(data.data)] = True

        #- extract only valid (not masked data)
        if mode == 'all':
            msk = np.sum(~data.mask,axis=0) == n #identify all ngrid cells where all timesteps are valid
        elif mode == 'one':
            msk = np.sum(~data.mask,axis=0) > 0  #identify ONE grid cell where all timesteps are valid
        else:
            print mode
            raise ValueError, 'Invalid option in get_valid_data()'

        data = data[:,msk]
        if lon is not None:
            lon  = lon[msk]
        if lat is not None:
            lat  = lat[msk]

        if return_mask:
            return lon,lat,data,msk
        else:
            return lon,lat,data

#-----------------------------------------------------------------------

    def _apply_mask(self,msk1,keep_mask=True):
        """
        apply a mask to C{Data}. All data where mask==True
        will be masked. Former data and mask will be stored

        @param msk: mask
        @type msk : numpy boolean array or Data

        @param keep_mask: keep old masked
        @type keep_mask : boolean
        """

        if isinstance(msk1,Data):
            msk = msk1.data
        else:
            msk = msk1

        self.__oldmask = self.data.mask.copy()
        self.__olddata = self.data.data.copy()

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


#-----------------------------------------------------------------------

    def shift_x(self,nx):
        """
        shift data array in x direction by nx steps

        @param nx: shift by nx steps
        @type nx: int
        """

        self.data = self.__shift3D(self.data,nx)
        self.lat  = self.__shift2D(self.lat,nx)
        self.lon  = self.__shift2D(self.lon,nx)

#-----------------------------------------------------------------------

    def __shift3D(self,x,n):
        """
        shift 3D data

        @param x: data to be shifted
        @type x: array(:,:,:)

        @param n: shifting step
        @type n: int
        """
        tmp = x.copy(); y=x.copy()
        y[:,:,:]=np.nan
        y[:,:,0:n] = tmp[:,:,-n:]
        y[:,:,n:]  = tmp[:,:,0:-n]

        return y

#-----------------------------------------------------------------------

    def timeshift(self,n,return_data = False):
        """
        shift data in time by n-steps
        positive numbers mean, that the
        data is shifted leftwards (thus towards earlier times)

        e.g. timeseries 1980,1981,1982, a shift of n=1 would generate
        a dataset with the data ordered as follows
        [1981,1982,1980]

        @param n: lag to shift data (n>=0)
        @type n: int

        @param return_data: specifies if a NEW C{Data} object shall be returned
        @type return_data: bool

        @todo: support for n < 0
        """

        if n == 0:
            return
        if self.data.ndim != 3:
            raise ValueError, 'array of size [time,ny,nx] is needed for temporal shifting!'

        #--- copy data
        tmp = self.data.copy()

        #--- generate output
        if return_data: #... a new data object is returned
            res           = self.copy()
            res.data[:,:,:]    = np.nan
            res.data[:-n:,:,:] = tmp[n:,:,:]
            res.data[-n:,:,:]  = tmp[0:n,:,:]
            return res
        else: #... the data object is changed
            self.data[:,:,:]    = np.nan
            self.data[:-n:,:,:] = tmp[n:,:,:]
            self.data[-n:,:,:]  = tmp[0:n,:,:]
            return None



#-----------------------------------------------------------------------

    def _set_valid_range(self,vmin,vmax):
        """
        sets the valid range of the data

        only data with vmin <= data <= vmax will be kept as valid

        @param vmin: minimum valid value
        @type vmin: float
        @param vmax: maximum valid value
        @type vmax: float
        """
        self.data = np.ma.array(self.data,mask = ((self.data < vmin) | (self.data > vmax))   )


#-----------------------------------------------------------------------

    def __shift2D(self,x,n):
        """
        shift 2D data

        @param x: data to be shifted
        @type x: array(:,:)

        @param n: shifting step
        @type n: int
        """
        tmp = x.copy(); y=x.copy()
        y[:,:]=np.nan
        y[:,0:n] = tmp[:,-n:]
        y[:,n:]  = tmp[:,0:-n]

        return y

#-----------------------------------------------------------------------

    def copy(self):
        """
        copy complete C{Data} object including all attributes
        @return C{Data} object
        """
        d = Data(None,None)


        for attr, value in self.__dict__.iteritems():
            try:
                #-copy (needed for arrays)
                cmd = "d." + attr + " = self." + attr + '.copy()'
                exec cmd
            except:
                #-copy
                cmd = "d." + attr + " = self." + attr
                exec cmd

        return d

#-----------------------------------------------------------------------

    def add(self,x,copy=True):
        """
        Add a C{Data} object to the current object field

        @param x: C{Data} object which will be added
        @type  x: Data object

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if np.shape(self.data) != np.shape(x.data):
            raise ValueError, 'Inconsistent geometry (add): can not calculate!'

        if copy:
            d = self.copy()
        else:
            d = self

        d.data = d.data + x.data
        d.label = self.label + ' + ' + x.label

        return d

#-----------------------------------------------------------------------

    def sub(self,x,copy=True):
        """
        Substract a C{Data} object from the current object field

        @param x: C{Data} object which will be substracted
        @type  x: Data object

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if np.shape(self.data) != np.shape(x.data):
            raise ValueError, 'Inconsistent geometry (sub): can not calculate!'

        if copy:
            d = self.copy()
        else:
            d = self

        d.data = d.data - x.data
        d.label = self.label + ' - ' + x.label

        return d

#-----------------------------------------------------------------------------------------------------------------------

    def diff(self,x,axis=0,equal_var=True,mask_data = False,pthres=0.05):
        """
        Difference between two C{Data} objects

        This routine calculates the difference between the data of two
        datasets. It calculates the significance and returns
        the mean differences and their corresponding significances

        The significance is calculated using a two-tailored t-test or a welch test
        in case of different variances. Independent samples are assumed!
        (unittest test_diff)

        @param x: C{Data} object which will be substracted from self
        @type  x: Data object

        @param axis: axis along which the data will be aggregated (typically axis=0 corresponds to time)
        @type axis: int

        @param equal_var: specifies if the two input datasets (self,x) are expected to have same variance.
                          dependent on this parameter. If the variance is equal, then a t-test is applied,
                          if not, then a welch test is used.
        @type equal_var: bool

        @param mask_data: specifies if the data field which is returned should be masked for all areas that do *not* show significant changes!
        @type mask_data: bool

        @param pthres: threshold for significant p-value; a value of e.g. 0.05 corresponds to the 95% significance level.
                       this threshold will be used for the generation of a mask that might be used e.g. as an overlay in map_plot()
        @type pthres: float

        @return: returns a C{Data} object that includes a) the difference map, b) the p-value, c) a mask that can be used e.g. as an overlay for map_plot()
        @rtype: C{Data}

        @todo: implementation of welch test. This should be actually already be implemented in stats.ttest_ind, but is not available in my python installation!
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html#scipy.stats.ttest_ind
        but the code would be available here: https://github.com/scipy/scipy/blob/v0.11.0/scipy/stats/stats.py#L2957


        """

        #/// check consistency
        if np.shape(self.data) != np.shape(x.data):
            raise ValueError, 'Inconsistent geometry (sub): can not calculate!'
        if axis > len(self.data.shape)-1:
            raise ValueError, 'Invalid axis parameter: ' + str(axis)
        if axis < 0:
            raise ValueError, 'Invalid axis parameter: ' + str(axis)

        #/// create new data object
        d = self.copy(); d.label = self.label + ' - ' + x.label

        #/// calculate statistical significance of the difference
        if isinstance(d.data,np.ma.masked_array):
            sta = stats.mstats
            t,p = ttest_ind(d.data, x.data ,axis=axis) #use routine in pyCMBS.statistic.py
        else:
            t,p = stas.ttest_ind(d.data, x.data ,axis=axis) #todo equal var for welch test not part of my psthon installation!

        p   = 1.- p #invert p-value, as a p-value of 1. would correspond to the same data

        #/// mean difference masked if p-value too low
        mask      = p <= pthres
        if mask_data:
            d.data    = np.ma.array(self.timmean() - x.timmean(),mask=~mask) #mean difference as masked array
        else:
            d.data    = np.ma.array(self.timmean() - x.timmean(),mask=np.zeros(self.timmean().shape).astype('bool') ) #mean difference as masked array
        d.p_value = p
        d.p_mask  = mask #masks the grid cells that show significant changes (todo check this again!) needs additional validation
        d.t_value = t

        return d


#-----------------------------------------------------------------------------------------------------------------------

    def subc(self,x,copy=True):
        """
        Substract a constant value from the current object field

        @param x: constant (can be either a scalar or a field that has
                  the same geometry as the second and third dimension of self.data)
        @type  x: float

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if copy:
            d = self.copy()
        else:
            d = self

        if np.isscalar(x):
            d.data -= x
        elif x.ndim == 2: #x is an array
            for i in xrange(len(self.time)):
                d.data[i,:,:] -= x[:,:]
        else:
            raise ValueError, 'Invalid geometry in detrend()'
        return d

#-----------------------------------------------------------------------

    def addc(self,x,copy=True):
        """
        Add a constant value to the current object field
        (unittest)

        @param x: constant
        @type  x: float

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if copy:
            d = self.copy()
        else:
            d = self
        d.data += x
        return d

#-----------------------------------------------------------------------

    def divc(self,x,copy=True):
        """
        Divide current data by a constant

        (unittest)

        @param x: constant
        @type  x: float

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if copy:
            d = self.copy()
        else:
            d = self
        d.data /= x
        return d

#-----------------------------------------------------------------------

    def div(self,x,copy=True):
        """
        Divide current object field by field of a C{Data} object

        (unittest)

        @param x: C{Data} object in the denominator
        @type  x: C{Data} object (data needs to have either same geometry
                  as self.data or second and third dimension need to match)

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if np.shape(self.data) != np.shape(x.data):
            if self.data.ndim == 3:
                if x.data.ndim == 2:
                    if np.shape(self.data[0,:,:]) == np.shape(x.data):
                        #second and third dimension match
                        pass
                    else:
                        print np.shape(self.data)
                        print np.shape(x.data)
                        raise ValueError, 'Inconsistent geometry (div): can not calculate!'
                else:
                        print np.shape(self.data)
                        print np.shape(x.data)
                        raise ValueError, 'Inconsistent geometry (div): can not calculate!'
            else:
                print np.shape(self.data)
                print np.shape(x.data)
                raise ValueError, 'Inconsistent geometry (div): can not calculate!'


        if copy:
            d = self.copy()
        else:
            d = self
        if np.shape(d.data) == np.shape(x.data):
            d.data /=  x.data
        elif np.shape(d.data[0,:,:]) == np.shape(x.data):
            for i in xrange(len(self.time)):
                d.data[i,:,:] /=  x.data
        else:
            raise ValueError, 'Can not handle this geometry in div()'

        d.label = self.label + ' / ' + x.label

        return d

#-----------------------------------------------------------------------

    def mul(self,x,copy=True):
        """
        Multiply current object field by field by a C{Data} object

        @param x: C{Data} object in the denominator
        @type  x: C{Data} object (data needs to have either same geometry
                  as self.data or second and third dimension need to match)

        @param copy: if True, then a new data object is returned
                     else, the data of the present object is changed
        @type copy: bool
        """

        if np.shape(self.data) != np.shape(x.data):
            if self.data.ndim == 3:
                if x.data.ndim == 2:
                    if np.shape(self.data[0,:,:]) == np.shape(x.data):
                        #second and third dimension match
                        pass
                    else:
                        print np.shape(self.data)
                        print np.shape(x.data)
                        raise ValueError, 'Inconsistent geometry (div): can not calculate!'
                else:
                        print np.shape(self.data)
                        print np.shape(x.data)
                        raise ValueError, 'Inconsistent geometry (div): can not calculate!'
            else:
                print np.shape(self.data)
                print np.shape(x.data)
                raise ValueError, 'Inconsistent geometry (div): can not calculate!'


        if copy:
            d = self.copy()
        else:
            d = self
        if np.shape(d.data) == np.shape(x.data):
            d.data = d.data * x.data
        elif np.shape(d.data[0,:,:]) == np.shape(x.data):
            for i in xrange(len(self.time)):
                d.data[i,:,:] = d.data[i,:,:] * x.data
        else:
            raise ValueError, 'Can not handle this geometry in div()'

        d.label = self.label + ' * ' + x.label

        return d

#-----------------------------------------------------------------------

    def _sub_sample(self,step):
        """
        subsample data of current C{Data} object

        @param step: stepsize for subsampling
        @type step: int
        """
        if self.data.ndim == 3:
            self.data = self.data[:,::step,::step]
        elif self.data.ndim == 2:
            self.data = self.data[::step,::step]
        else:
            raise ValueError, 'Data Dimension not supported!'
        self.lat  = self.lat [::step,::step]
        self.lon  = self.lon [::step,::step]

#-----------------------------------------------------------------------

    def corr_single(self,x,pthres=1.01,mask=None):
        """
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

        @return: list of C{Data} objects for correlation, slope, intercept, p-value, covariance
        @rtype: list

        @todo significance of correlation (is the calculation correct? currently it assumes that all data is valid!)
        @todo: it is still not ensured that masked data is handled properly!!!
        """

        if self.data.ndim != 3:
            raise ValueError, 'Invalid geometry!'

        nt,ny,nx = sz = np.shape(self.data)

        if nt != len(x):
            raise ValueError, 'Inconsistent geometries'

        #--- get data with at least one valid value
        lo,la,dat,msk = self.get_valid_data(return_mask=True,mode='one')
        xx,n = dat.shape
        #~ print msk.shape, sum(msk)
        if self.verbose:
            print '   Number of grid points: ', n
        #~ print dat.shape

        R=np.ones((ny,nx))*np.nan #output matrix for correlation
        P=np.ones((ny,nx))*np.nan #output matrix for p-value
        S=np.ones((ny,nx))*np.nan #output matrix for slope
        I=np.ones((ny,nx))*np.nan #output matrix for intercept
        CO=np.ones((ny,nx))*np.nan #output matrix for covariance

        R.shape = (-1); S.shape = (-1)
        P.shape = (-1); I.shape = (-1)
        CO.shape = (-1)

        print 'Calculating correlation ...'
        res = [stats.mstats.linregress(x,dat[:,i]) for i in range(n)] #@todo: still rather inefficient for masked arrays
        res = np.asarray(res)

        slope = res[:,0]; intercept = res[:,1]
        r_value = res[:,2]; p_value = res[:,3]
        std_err = res[:,4]

        R[msk] = r_value; P[msk] = p_value; I[msk] = intercept; S[msk] = slope
        R.shape = (ny,nx); P.shape = (ny,nx); I.shape = (ny,nx); S.shape = (ny,nx)

        #--- prepare output data objects
        Rout = self.copy() #copy obkect to get coordinates
        Rout.label = 'correlation'
        msk = (P > pthres) | np.isnan(R)
        Rout.data = np.ma.array(R,mask=msk).copy()

        Sout = self.copy() #copy object to get coordinates
        Sout.label = 'slope'
        Sout.data = np.ma.array(S,mask=msk).copy()

        Iout = self.copy() #copy object to get coordinates
        Iout.label = 'intercept'
        Iout.data = np.ma.array(I,mask=msk).copy()

        Pout = self.copy() #copy object to get coordinates
        Pout.label = 'p-value'
        Pout.data = np.ma.array(P).copy()

        Cout = self.copy() #copy object to get coordinates
        Cout.label = 'covariance'
        Cout.data = np.ma.array(np.ones(P.shape)*np.nan,mask=msk).copy() #currently not supported: covariance!

        if mask is not None:
            #apply a mask
            Rout._apply_mask(mask); Sout._apply_mask(mask)
            Iout._apply_mask(mask); Pout._apply_mask(mask)
            Cout._apply_mask(mask)

            Rout.unit = None; Sout.unit = None; Iout.unit = None
            Pout.unit = None; Cout.unit = None

        return Rout,Sout,Iout,Pout, Cout

#-----------------------------------------------------------------------

    def detrend(self,return_object=True):
        """
        detrend data timeseries by removing linear trend over time.
        It is assumed that the timesamples have equidistant spacing.
        This assumption is important to consider, as regression is calculated
        only using the length of the data vector!

        @param return_object: specifies if C{Data} object will be returned (default0True)
        @type return_object: bool
        """

        print 'Detrending data ...'

        if self.data.ndim != 3:
            raise ValueError, 'Can not detrend data other than 3D!'

        #generate dummy vector for linear correlation (assumes equally spaced data!!!!!) todo: generate unittest for this
        x = np.arange(len(self.time)) #@todo: replace this by using actual timestamp for regression calcuclation

        #correlate and get slope and intercept
        Rout,Sout,Iout,Pout, Cout = self.corr_single(x)

        #calculate regression field
        reg = Data(None,None)
        reg.data = np.zeros(self.data.shape) * np.nan
        reg.label = 'trend line'
        nt = self.data.shape[0]
        for i in range(nt):
            reg.data[i,:,:] = Sout.data * i + Iout.data

        #substract regression line
        res = self.sub(reg)
        res.label = self.label + '(det.)'

        if return_object:
            res.detrended = True
            return res
        else:
            self.data = res.data
            self.detrended = True
            return None

#-----------------------------------------------------------------------


