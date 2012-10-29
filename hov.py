#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


'''
HOVMOELLER CLASS
class to generate hovmoeller diagrams
'''

from pylab import *
import pylab as pl
import matplotlib.dates as mdates
import sys
import matplotlib.pyplot as pyplot

def agg_hourly(d,v,timestamp='mid',mode='mean'):
    '''
    calculate hourly mean values in a very efficient way

    INPUT
       d - array of datetime objects
       v - array of values

       optional arguments
       timestamp - specifies if timestamp is at beginning / mid or end of hour
       --> first / mid / last

    OUTPUT
       returns pandas timeseries object

    example from
    http://stackoverflow.com/questions/6467832/how-to-get-the-correlation-between-two-timeseries-using-pandas
    '''

    import pandas as pa
    s = pa.Series(v,index=d)
    #mean hourly values using pandas group function
    if mode == 'mean':
        #r = s.groupby(lambda date: date.replace(second=0,microsecond=0,minute=0)).mean() #works only for most recent version of pandas
        r = s.groupby(lambda date: date.replace(second=0,microsecond=0,minute=0)).agg(mean)
    else:
        sys.exit('agg_hourly - invalid mode: ' + mode)

    if timestamp == 'first':
        pass
    elif timestamp == 'mid':
        r=r.rename(lambda date: date.replace(minute=30))
    elif timestamp == 'last':
        r=r.rename(lambda date: date.replace(minute=59,second=59))
    else:
        sys.exit('agg_hourly: invalid timestamp option' )
    return r


def align(x,y):
    '''
    aligns two pandas timeseries and returns
    values two pandas timeseries which are aligned

    this is a hack as I dont know how to do it better at the moment
    '''

    t=x+y #just add the two series and then substract the individuals again

    return t-y,t-x #corresponds to x,y


def generate_monthly_timeseries(t,sday='01'):
    '''
    generate a vector monthly timeseries

    t: time (numeric)
    '''


    tn=[]; d=num2date(t)
    for i in range(len(t)):
        y = str(d[i].year).zfill(4)
        m = str(d[i].month).zfill(2)
        s = y + '-' + m + '-' + sday

        tn.append(datestr2num(s))

    return tn





class hovmoeller:
    def __init__(self,time,value,var_unc=None,rescalex=1,rescaley=1,lat=None,lon=None,transpose=False):
        '''
        construction of hovmoeller class

        time: datetime structure
        value: timeseries to be used

        if the argument lat is provided it is assumed that lat/lon are 2D matrices
        In this case the value is expected to be a 3D variables as
        value(time,ny,nx)


        EXAMPLES:

        1) a minimum example for a time-latitude plot
        #get data from e.g. file; here random data
        t1=datestr2num('2011-05-01'); t2=datestr2num('2011-08-31')

        #create a hovmoeller object
        myhov = hovmoeller(d,dat,lat=LATS,rescaley=10,rescalex=10)
        ## hov2.time_to_lat(dlat=200.,yticksampling=1)
        ## hov2.plot(title='HI alex',ylabel='lat',xlabel='days',origin='lower',xtickrotation=30)
        ## show()


        var: additional array with variance information (can be used to show uncertainties)

        '''

        #/// check consistency ///
        if len(time) == len(value):
            pass
        else:
            sys.exit('Inconsistent sizes of time and value (hovmoeller)')

        if lat != None:
            if shape(lat) != shape(value[0,:,:]):
                print shape(lat), shape(value[0,:,:])
                sys.exit('Inconsistent latitudes and data (hovmoeller)')


        #/// set values of class ///
        self.time = time
        self.transpose=transpose
        ntim = len(self.time)

        t=date2num(time)
        self.t_min=num2date(floor(t.min()))
        self.t_max=num2date(ceil(t.max()))

        self.value = value.copy()
        self.value = ma.array(self.value,mask=isnan(self.value))
        self.value.shape = (ntim,-1) #now reshape everything to [time,ngridcells]

        if var_unc == None:
            self.var_unc = None
        else:
            self.var_unc = var_unc.copy()
            self.var_unc.shape = (ntim,-1)

        self.hov_var = None

        self.rescalex = rescalex
        self.rescaley = rescaley

        if lat != None:
            self.lat = lat.copy()
            self.lat.shape = (-1)
        else:
            self.lat = None


        if lon != None:
            self.lon = lon.copy()
            self.lon.shape = (-1)
        else:
            self.lon = None



    def plot(self,xticks=None,xlabel=None,ylabel=None,title='',grid=True,climits=None,figsize=None,origin=None,xtickrotation=0,cmap='jet',showcolorbar=True,ax=None,show_uncertainties=False,norm_uncertainties=False) :

        '''
        plot result
        clim: tuple

        norm_uncertainties: divide value by variance ==> contourplot of self.hov / self.hov_var is generated
        '''


        if climits == None:
            sys.exit('Hovmoeller, please specify climits')

        if xlabel == None:
            self.xlabel='x-label'
        else:
            self.xlabel=xlabel

        if figsize == None:
            if self.transpose:
                figsize=(6,11)
            else:
                figsize=(12,4)

        if ylabel == None:
            self.ylabel='y-label'
        else:
            self.ylabel=ylabel

        self.title=title


        if ax == None:
            self.fig = figure(figsize=figsize)
            self.ax  = self.fig.add_subplot(111)
        else:
            self.ax = ax
            self.fig = self.ax.figure

        if self.transpose:
            arr = self.hov.repeat(self.rescalex,axis=0).repeat(self.rescaley,axis=1)
            #~ print 'HOVMOELLER SHAPE ARRAY1: ', arr.shape
            self.im=self.ax.imshow(arr.T,interpolation='Nearest',origin=origin,vmin=climits[0],vmax=climits[1],cmap=pyplot.get_cmap(cmap) )
        else:
            #rescale array for visualization
            arr = self.hov.repeat(self.rescaley,axis=0).repeat(self.rescalex,axis=1)
            #~ print 'HOVMOELLER SHAPE ARRAY2: ', arr.shape
            self.im=self.ax.imshow(arr,interpolation='Nearest',origin=origin,cmap=pyplot.get_cmap(cmap),vmin=climits[0],vmax=climits[1])

            if (show_uncertainties) & (self.hov_var != None):
                arr1 = self.hov_var.repeat(self.rescaley,axis=0).repeat(self.rescalex,axis=1)

                if norm_uncertainties:
                    #normalized by variance
                    arr1 = arr / arr1
                self.ax.contour(arr1,linestyles='-',colors='black')


        if xlabel != None: self.ax.set_xlabel(self.xlabel)
        if ylabel != None: self.ax.set_ylabel(self.ylabel)
        #self.fig.colorbar()
        self.ax.set_title(self.title)

        #/// ticks
        if self.transpose:
            self.ax.yaxis.set_major_locator(self.x_major_locator)
            self.ax.yaxis.set_major_formatter(self.x_major_formatter)

            self.ax.xaxis.set_major_locator(self.y_major_locator)
            self.ax.xaxis.set_major_formatter(self.y_major_formatter)
        else:
            self.ax.yaxis.set_major_locator(self.y_major_locator)
            self.ax.yaxis.set_major_formatter(self.y_major_formatter)

            self.ax.xaxis.set_major_locator(self.x_major_locator)
            self.ax.xaxis.set_major_formatter(self.x_major_formatter)

        pl.xticks(rotation=xtickrotation)

        if xticks != None:
            nx = 2
            set_label(nx)

        if grid: #show grid
            self.ax.grid(color='white',linewidth=1,linestyle='-')

        #~ if climits == None:
            #~ pass
        #~ else:
            #~ clim(climits)

        #~ self.ax.set_ylim(top=5*self.rescaley, bottom=19*self.rescaley) #this CAN NOT WORK!



        #/// colorbar
        if showcolorbar:
            if self.transpose:
                self.fig.colorbar(self.im,orientation='vertical',shrink=0.5)
            else:
                self.fig.colorbar(self.im,orientation='horizontal',shrink=0.5,aspect=30)


        #        cax = fig.add_axes([0.2, 0.05, 0.5, 0.03])
        #fig.colorbar(im[0], cax, orientation='horizontal')

    def show(self):
        self.fig.show()


    def set_label(self,nx):
        '''
        return labels
        '''
        xticks = self.ax.xaxis.get_xticks()
        nticks = len(xticks)

        #self.ax.xaxis.set_xticks([])


    def time_to_lat(self,dlat=1.,mode='average',monthsamp=1,yticksampling=1,monthly=False,yearonly=False  ):
        '''
        convert timeseries and latitude to a hovmoeller matrix

        dlat: sampling of latitudes [degree]
        mode: method how to aggregate the data
        monthsamp: method how to sample the data

        monthly: aggregated data to monthly (always mid of month)

        todo: implement more properly handling with masked arrays!

        todo: implement area weighting!
        '''

        self.yearonly = yearonly

        if self.lat == None:
            sys.exit('Error time_to_lat in hovmoeller: no latitude specified')

        #/// preprocessing: extract only VALID data
        value = self.value.copy()
        lat   = self.lat.copy()

        #remove all lats where all data is invalid
        npix = len(lat)
        if shape(value)[1] != npix:
            print 'Invalid geometry time_to_lat', npix,shape(value)
            sys.exit()

        pixmsk=ones(npix)
        for i in range(npix):
            if all(value[:,i].mask): #check if whole timeseries is masked
                pixmsk[i]=0.

        #~ print len(pixmsk),sum(pixmsk)

        #print pixmsk

        #/// generate latitudes by rounding
        lats = ((lat - 0.5*dlat) / dlat).round() * dlat
        ulats = unique(lats[pixmsk==1.]); nlat = len(ulats)
        #~ print 'Number of lats: ', nlat
        print ulats

        d = self.time; tnum = date2num(d)
        t_min = date2num(self.t_min)
        t_max = date2num(self.t_max)

        data_days = floor(date2num(d))
        days = floor(date2num(d))

        if monthly:
            #/// monthly ///
            data_days = generate_monthly_timeseries(tnum) #convert times to monthly
            all_days  = unique(data_days)
        else:
            #/// per days ///
            all_days = linspace(t_min,t_max,t_max-t_min+1) #every day

        #/// init output arrays ///
        outsum = zeros((nlat,len(all_days)))
        outn   = zeros((nlat,len(all_days)))



        #/// loop over all days
        tidx = arange(len(data_days))
        for i in range(len(all_days)):
            m = data_days == all_days[i] #multiple days possible if sub-day sampling!
            #~ print sum(m)

            actidx = list(tidx[m]) #[0]
            #~ print 'INDEX: ', actidx, array(value[actidx,:]).shape
            #~ tmp = value[actidx,:].copy()
            #~ print mean(tmp), tmp.mean()
            v = mean(value[actidx,:].copy(),axis=0) #average all timestamps of same day
            thelats = lats.copy()
            #print i,actidx, unique(v), sum(m),tidx[m]
            #print shape(mean(v,axis=0))


            #/// loop over all latitudes
            for j in range(len(ulats)):
                ml = thelats == ulats[j] #determine all pixels for a given lat ...
                #print 'shape ml: ', shape(ml), shape(v)
                if sum(ml) > 0:
                    v1 = v[ml]
                    v1 = v1[~isnan(v1)] #only valid data

                    if hasattr(v1,'mask'):
                        #~ print 'has mask v1'
                        v1 = v1.data[~v1.mask] #subset data only for valid data; important  to calculate approp. the mean value!

                    #~ print ''
                    #~ print num2date(all_days[i])
                    #~ print 'latmask', sum(ml)
                    #~ print v1

                    #print unique(v1), ulats[j]
                    if len(v1)>0.:
                        outsum[j,i] = outsum[j,i] + sum(v1) #... and average them
                        outn  [j,i] = outn  [j,i] + len(v1)
                    #~ print sum (v1), outsum[j,i], outn[j,i]
                else:
                    print 'Da sollten wir beim Testen nicht landen todo'


        if mode == 'average':
            out = outsum / outn
        else:
            sys.exit('Unknown average mode in time_to_lat')

        #/// assign data matrix ///
        self.hov = out.copy()

        #/// labels or y-axis uses matplotlib.ticker
        lattick = arange(ulats.min(),ulats.max()+dlat,dlat)
        #~ print lattick
        #~ print ulats.min(),ulats.max(),dlat

        yticklabels = map(str,lattick)

        if self.transpose:
            scal = self.rescalex
        else:
            scal = self.rescaley

        yticks = ((lattick - lattick.min()) / dlat )     * scal

        #/// sumsampling of ticks
        yticks = yticks[::yticksampling]
        yticklabels = yticklabels[::yticksampling]

        self.y_major_locator=FixedLocator(yticks)
        self.y_major_formatter=FixedFormatter(yticklabels)

        #/// x-ticks and labels ///
        self.generate_xticks(all_days,monthsamp=monthsamp)


    def time_to_day_hour2(self,mode='average',monthsamp=1,yticksampling=1):
        d=self.time
        v=self.value


    def time_to_day_hour_fast(self,yticksampling=1,monthsamp=1,yearonly=False):
        '''
        routine to make a fast calculation of hovmoeller data using pandas
        '''

        print 'time_to_day_hour_fast not checked for new geometry yet!'
        stop

        import pandas as pa1
        arr = None

        self.yearonly = yearonly

        #1) generate hourly timeseries from original data and store it in a dataframe
        r = agg_hourly(self.time,self.value)
        df_data = pa1.DataFrame(r.values,index=r.index,columns=['data'])

        dmin=df_data.index.min(); dmax=df_data.index.max()
        t1 = str(dmin.month).zfill(2) + '/' + str(dmin.day).zfill(2) + '/' + str(dmin.year).zfill(4)
        t2 = str(dmax.month).zfill(2) + '/' + str(dmax.day).zfill(2) + '/' + str(dmax.year).zfill(4)

        for i in range(24):
            #2) generate a vector of days for a particular hour and store it in a dataframe
            hr=str(i).zfill(2)
            d_hours   = pa1.DateRange(t1 + ' ' + hr + ':30:00+00:00', t2 + ' ' + hr + ':30:00+00:00', offset=pa1.datetools.Hour(n=24))
            df_hours  = pa1.DataFrame(rand(len(d_hours))*nan,index=d_hours,columns=['nix'])  #reference data frame

            #3) merge the two dataframes: result is a tiemseries of the same length as the data, but only for the single hour
            df_j = df_hours.join(df_data)

            #4) store results for actual hour in output array
            if arr==None:
                ndays  = len(df_j.values)
                nhours = 24
                arr=zeros((nhours,ndays))*nan #create output array
            else:
                pass

            arr[i-1,:]=df_j['data'].values

        self.hov = arr

        #/// labels or y-axis uses matplotlib.ticker
        all_days = date2num(df_hours['nix'].index)

        if self.transpose:
            scal = self.rescalex
        else:
            scal = self.rescaley

        yticks = array([0,6,12,18,24]) * scal #array for y-ticks
        yticklabels = ['00','06','12','18','24']

        #/// sumsampling of ticks
        yticks = yticks[::yticksampling]
        yticklabels = yticklabels[::yticksampling]

        self.y_major_locator=FixedLocator(yticks)
        self.y_major_formatter=FixedFormatter(yticklabels)

        #/// x-ticks and labels ///
        self.generate_xticks(all_days,monthsamp=monthsamp)




    def time_to_day_hour(self,mode='average',monthsamp=1,yticksampling=1,yearonly=False):
        '''
        convert timeseries to hovemoeller matrix
        '''

        print 'time_to_day_hour not checked for new geometry yet!'
        stop

        self.yearonly = yearonly

        d=self.time

        tnum=date2num(d)
        all_days = unique(floor(date2num(d)))
        t_min = date2num(self.t_min)
        t_max = date2num(self.t_max)
        all_days = linspace(t_min,t_max,t_max-t_min+1)
        days = floor(date2num(d))

        out = ones((24.,len(all_days)))*nan
        outsum = zeros((24.,len(all_days)))

        if (self.var_unc != None):
            if (shape(self.var_unc) == shape(self.value)):
                f_var = True
            else:
                print 'HOV: inconsistent values for data array and variance array!', shape(self.var_unc), shape(self.value)
                sys.exit()
        else:

            f_var = False

        if f_var: #process also variance information
            outvar = out.copy()

        for i in range(len(all_days)):
            m = days == all_days[i]
            t1 = tnum[m] #time of a particular day
            v = self.value[m]
            if f_var:
                vars = self.var_unc[m]

            for j in range(len(t1)):
                t = num2date(t1[j])
                if isnan(out[t.hour,i]):
                    out[t.hour,i] = v[j]
                    outsum[t.hour,i] = 1.
                    if f_var:
                        outvar[t.hour,i] = vars[j]
                else:
                    out[t.hour,i] = out[t.hour,i] + v[j]
                    if f_var:
                        outvar[t.hour,i] = outvar[t.hour,i] + vars[j] #variance of the sum is sum of variances (assumed uncorrelated data)

                    outsum[t.hour,i] += 1.

        if mode == 'average':
            out = out / outsum #average value per hour
        else:
            print 'Unknown mode! (hov)'
            stop

        self.hov = out
        if f_var:
            self.hov_var = outvar

        #/// labels or y-axis uses matplotlib.ticker
        if self.transpose:
            scal = self.rescalex
        else:
            scal = self.rescaley

        yticks = array([0,6,12,18,24]) * scal #array for y-ticks
        yticklabels = ['00','06','12','18','24']

        #/// sumsampling of ticks
        yticks = yticks[::yticksampling]
        yticklabels = yticklabels[::yticksampling]

        self.y_major_locator=FixedLocator(yticks)
        self.y_major_formatter=FixedFormatter(yticklabels)

        #/// x-ticks and labels ///
        self.generate_xticks(all_days,monthsamp=monthsamp)


    def generate_xticks(self,all_days,monthsamp=1):
        '''
        generate xticks

        all_days: array of times (days since 01-0-0)
        '''

        #~ print all_days

        dd=num2date(all_days)
        xticks=[]
        xticklabels=[]
        last_month = dd[0].month
        last_year  = dd[0].year

        xticks.append(all_days[0])
        strmon = str(dd[0].month).zfill(2)
        #if dd[0].month < 10:
        #    strmon = '0' + strmon
        if self.yearonly:
            xticklabels.append(str(dd[0].year))
        else:
            xticklabels.append(strmon + '/' + str(dd[0].year))
        for d in dd:
            #~ print d, last_month
            if (d.month != last_month) | (d.year != last_year):
                #~ print d
                xticks.append(date2num(d)) #save tick location
                mstr = str(d.month).zfill(2)
                #if d.month < 10:
                #    mstr = '0'+mstr
                if self.yearonly:
                    xticklabels.append(str(d.year))
                else:
                    xticklabels.append(mstr + '/' + str(d.year))
                last_month = d.month
                last_year  = d.year

        #~ print 'XTICKS: ', xticks

        if self.transpose:
            scal = self.rescaley
        else:
            scal = self.rescalex
        xticks=asarray(xticks)

        #/// remap to image coordinates
        ny,nx = shape(self.hov)
        w = (xticks - date2num(self.t_min)) / (date2num(self.t_max) - date2num(self.t_min)) #weights for positions on x-axis
        xticks = w * nx * scal; xticks=list(xticks)

        #here we have monthly ticks which we can now subsample
        xticks = xticks[::monthsamp]
        xticklabels=xticklabels[::monthsamp]


        #~ print 'xticklabels: ', xticklabels

        self.x_major_locator   = FixedLocator(xticks)
        self.x_major_formatter = FixedFormatter(xticklabels)



        #xlabels could be made with matplotlib.dates



## '''
## This is a sample program hov to use the hovmoeller class
## '''
#~
#~ close('all')
#~
#~ #generate some time series with gaps
#~ t=linspace(0,365*2.,365*2*24)
#~ t[2000:5000]=nan
#~ m=t-floor(t)>0.75
#~ t[m]=nan
#~ m=~isnan(t)
#~ t=t[m]
#~ d=datestr2num('2001-04-01')+t
#~ d=num2date(d)
#~ y=sin(t*12*pi/365.)
#~
#~
#~ v=y*0.1 * rand(len(y))
#~
#~
#~ #/// adn plot is as a hovmoeller diagram ///
#~ hov = hovmoeller(d,y,rescaley=10,transpose=False,var=v)
#~ #hov = hovmoeller(d,y,rescalex=10,transpose=True)
#~
#~
#~
#~
#~
#~ #hov.time_to_day_hour_fast()
#~ #hov.plot(title='HI alex',ylabel='hour',xlabel='days',climits=(-1.,1.))
#~
#~ figure()
#~
#~
#~
#~
#~ hov.time_to_day_hour()
#~ hov.plot(title='HI alex',ylabel='hour',xlabel='days',climits=(-1.,1.),show_uncertainties=True,norm_uncertainties=True )
#~
#~ stop
#~
#~
#~ figure()
#~ subplot(211)
#~ imshow(hov.hov); colorbar()
#~ subplot(212)
#~ imshow(hov.arr2); colorbar()
#~
#~ stop
#~
#~ hov.plot(title='HI alex',ylabel='hour',xlabel='days',climits=(-1.,1.))
#~


## #/// test for latitude hovmoeller diagram
## #result should be a latitude plot
## lat=rand(len(y))*90.-30.
## y1 = lat.copy()


## hov1 = hovmoeller(d,y1,rescaley=10,lat=lat)
## hov1.time_to_lat(dlat=5.)
## hov1.plot(title='HI alex',ylabel='lat',xlabel='days',origin='lower',xtickrotation=30,cmap='RdBu_r')

## hov2 = hovmoeller(d,y1,lat=lat,rescalex=20,transpose=True)
## hov2.time_to_lat(dlat=5.,yticksampling=5)
## hov2.plot(title='HI alex',ylabel='lat',xlabel='days',origin='lower',xtickrotation=30,cmap='RdBu_r')




## hov.show()
## hov1.show()
## hov2.show()


#test for 2D resampling
## import Nio
## close('all')

## filename='/media/Data/Data/ISCCP/sub.nc'
## F=Nio.open_file(filename,'r')

## dat = F.variables['ci'][:,:,:].copy() #read a 3D data set
## dummy = F.variables['ci']._FillValue
## dat[dat==dummy]=nan
## t=F.variables['time'][:]
## #dat[dat==dummy]=nan

## nt,ny,nx=shape(dat)
## #nt = 100

## #dummy lat/lon
## lats=arange(ny)*10.
## lons=arange(nx)*100.
## LONS,LATS=meshgrid(lons,lats)

## #dummy data for testing
## #dat=[]; thetime=zeros(nt)  #thetime=zeros((nt,ny,nx))
## #for i in range(nt):
## #    dat.append(LATS)
## #    thetime[i]=t[i]

## dat=dat[0:nt]
## thetime=t[0:nt]

## dat=asarray(dat) #data for different time steps
## thetime=asarray(thetime)

## #thelons=LONS.flatten(); thelats=LATS.flatten()
## #thetime=thetime.flatten()

## #thelons=thelons.repeat(nt)
## #thelats=thelats.repeat(nt)

## #dat=dat.flatten()
## d=num2date(thetime)
## LONS,LATS = meshgrid(lons,lats)


## hov2 = hovmoeller(d,dat,lat=LATS,rescaley=10,rescalex=10)
## hov2.time_to_lat(dlat=200.,yticksampling=1)
## hov2.plot(title='HI alex',ylabel='lat',xlabel='days',origin='lower',xtickrotation=30)
## show()

## #map 2D field to vector for each time step
