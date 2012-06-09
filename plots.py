#!/usr/bin/pythonl
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"


from matplotlib import pylab as plt

from matplotlib.patches import Polygon

from mpl_toolkits.basemap import Basemap,shiftgrid

from scipy import stats

import numpy as np

from matplotlib.patches import Circle

import sys

from diagnostic import *

from scipy.spatial import cKDTree as KDTree #import the C version of KDTree (faster)

#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes 


#todo:
# - implement Glecker plots
# - implement writing of statistics to an ASCII file as export


class CorrelationAnalysis():
        '''
        perform correlation analysis
        and plot results
        '''
        
        def __init__(self,X,Y,mask=None,ax=None):
            self.x = X
            self.y = Y
            self.mask = mask
            
            if ax == None:
                f = plt.figure()
                self.ax = f.add_subplot(111)
            else:
                self.ax = ax
            
        def do_analysis(self):
            '''
            perform correlation analysis
            '''
            
            #todo: implement area weighting
            
            
            #--- calculate diagnostics
            D = Diagnostic(self.x,y=self.y)
            D._mat2vec(mask = self.mask) #here is the point fo rregional statistics
            rmse = D.get_rmse_value()
            r    = D.get_correlation_value()
            n    = D. get_n()
            
            print 'RMSE: ', rmse
            print 'R   : ', r
            print 'N   : ', n
        




class ReichlerPlot():
    def __init__(self,ax=None):
        '''
        plotting Reichler index
        '''
        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax
            
        self.e2 = [] #list to store RMS error results
        self.labels = []
        self.colors=[]
            
    def add(self,e2,label,color=None):
        self.e2.append(e2)
        self.labels.append(label)
        self.colors.append(color)
        
    def bar(self,vmin=None,vmax=None,title='',**kwargs):
        '''
        generate barplot
        '''
        print 'Doing Reichler plot as barplot ...'
        self._normalize()
        x = np.arange(len(self.e2_norm))
        self.ax.bar(x,self.e2_norm*100.,**kwargs)
        self.ax.set_xticks(x+0.5)
        self.ax.set_xticklabels(self.labels)
        if (vmin !=None) & (vmax != None):
            self.ax.set_ylim(vmin,vmax)
        self.ax.set_ylabel('relative error \n to multimodel \n mean [%]')
        self.ax.grid()
        self.ax.set_title(title)
        
        
    def simple_plot(self):
        for i in np.arange(len(self.e2)):
            self.ax.plot(self.e2[i],'o',label=self.labels[i])
        
    def circle_plot(self):
        print 'Doing Reichler plot as circle plot ...'
        self._normalize()
        
        dx=0.
        tsize=10.
        for i in np.arange(len(self.e2)): #over all timestamps
            print i, self.labels[i], self.e2_norm[i]*100.
            #~ print self.e2_norm
            circle = Circle( (self.e2_norm[i],0.), 0.1)

            circle.set_color(self.colors[i])
            circle.set_alpha(0.4)
            circle.set_label(self.labels[i])
            circle.set_edgecolor('k')
            self.ax.add_artist(circle)
            
            self.ax.text(0.1+dx, 0.5, self.labels[i], size=tsize, rotation=0.,
             ha="center", va="center",
             bbox = dict(boxstyle="round",
                         ec=self.colors[i],
                         fc=self.colors[i],
                         alpha = 0.4,
                         ))
            dx += 0.15
            
            
            
            
        self.ax.set_ylim(-1.,1.)
        self.ax.set_xlim(-1.,1.)
        self.ax.set_xlabel('relative error to multimodel mean')
        self.ax.legend()
        

        
        
    def _normalize(self):
        '''
        normalize results from different models
        Glecker et al, eq. 2
        '''
        
        n  = len(self.e2[0])
        E2 = []
        
        for e in self.e2:
            if len(e) != n:
                print 'WARNING: non consistent length in error statistics!!!'
            E2.append(np.nansum(e)) #temporal aggregation
        
        E2 = np.asarray(E2);  EM = E2.mean()
        self.e2_norm =  (E2 - EM) / EM #see Glecker et al, eq.2
        
        
########################################################################        

class ScatterPlot():
    
    def __init__(self,x,ax=None):
        '''
        x data object. This will be used for the x-variable
        '''
        
        
        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax
        
        self.figure = self.ax.figure
        self.x = x
        self.lines = []
        self.labels = []
        

        
    def plot(self,y,regress=True,**kwargs):
        label=y.label
        
        xdat = self.x.fldmean(); ydat = y.fldmean()
        
        if regress: #calculate linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(xdat,ydat)
            label = label + ' (r=' + str(round(r_value,2)) + ', p=' + str(round(p_value,2)) + ')'

        l = self.ax.plot(xdat,ydat,'.',label=label,**kwargs)[0]
        if regress:
            self.ax.plot(xdat,xdat*slope+intercept,'--',color=l.get_color()) 
        self.lines.append(l); self.labels.append(label)
        
        self.ax.set_xlabel(self.x._get_label() )
        self.ax.set_ylabel(y._get_unit())
        
    def legend(self):
        '''
        plot legend
        '''
        self.ax.legend(self.lines,self.labels,prop={'size':8})


class LinePlot():
    '''
    class for a pyCMBS Line Plot
    '''
    def __init__(self,ax=None,regress=False,title=None,show_xlabel=True,show_ylabel=True):
        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax
        
        self.figure = self.ax.figure
        self.regress = regress
        self.title = title
        self.show_xlabel = show_xlabel
        self.show_ylabel = show_ylabel
        
        self.lines = []
        self.labels = []
        
    def legend(self):
        '''
        plot legend
        '''
        self.ax.legend(self.lines,self.labels,prop={'size':8})
        
        
        
    def plot(self,x,ax=None,vmin=None,vmax=None,label = None, **kwargs):
        '''
        x: object of Data class
        '''
        
        if len(x.time) > 0:
            
            if ax == None:
                ax = self.ax
            else:
                ax = ax
            
            y = x.fldmean()
            if label == None:
                label = x.label
            
            if self.regress: #calculate linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(x.time,y)
                label = label + ' (r=' + str(round(r_value,2)) + ', p=' + str(round(p_value,2)) + ')'
                
            self.labels.append(label)

            p = ax.plot(plt.num2date(x.time), y , label=label, **kwargs)[0]
            self.lines.append(p)
            if self.regress:
                ax.plot(x.time,x.time*slope+intercept,'--',color=p.get_color()) #plot regression line

            if self.show_ylabel:
                ax.set_ylabel(x._get_unit())
            if self.show_xlabel:
                ax.set_xlabel('time')
            
            if self.title != None:
                ax.set_title(self.title)
                
            if vmin != None:
                if vmax != None:
                    ax.set_ylim(vmin,vmax)
            

            
            

class ZonalPlot():
    def __init__(self,ax=None,dir='y'):
        '''
        constructor for zonal plot
        
        dir - specifies direction for aggregation: y = zonal, x = meridional aggregation
        
        CAUTION: the class simply aggregates x/y. Thus the user needs to ensure, that the data is projected
        in a way that all lon/lat are on the same row/col
        '''
        
        #--- directionalities
        if dir == 'y': #zonal plot
            self.dir = 'y'
        elif dir == 'x':
            self.dir = 'x'
        else:
            raise ValueError, 'Invalid value for agregation direction (ZonalPlot): ', dir
        
        #--- set axis
        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

    def plot(self,x,areaweights,xlim=None):
        '''
        plot zonal plot
        x : data object
        '''
        
        #check if all latitudes are the same
        lu = x.lat.mean(axis=1)
        if any( abs(lu - x.lat[:,0]) > 1.E-5):
            print 'WARNING: latitudes are not unique!!!'
            print lu.shape,x.lat.shape
            print lu
            
            print x.lat[:,0]
            print x.lat[:,0] - lu
            
            stop
        
        if self.dir == 'y':
            dat = x.get_zonal_statistics(areaweights) #no area weighting performed
        else:
            raise ValueError, 'Invalid option'
            
        if dat.shape[x.data.ndim-2] != x.lat.shape[0]:
            print 'Inconsistent shapes!'
            print dat.shape
            print x.lat.shape
            sys.exit()
        
        
        #--- plot zonal statistics
        if dat.ndim == 1:
            self.ax.plot(dat,x.lat[:,0])
        elif dat.ndim == 2:
            for i in range(len(dat)):
                print 'Time in zonal: ', i
                print dat[i,:]
                #~ self.ax.plot(dat[i,:],label='time='+str(i))
                self.ax.plot(dat[i,:],x.lat[:,0],label='time='+str(i))
        
        self.ax.set_ylabel('latitude [deg]')
        self.ax.set_ylim(-90.,90.)
        
        if xlim != None:
            self.ax.set_xlim(xlim)
        
        self.ax.grid()
        


#=======================================================================

def __basemap_ancillary(m):
    latvalues=np.arange(-90.,120.,30.)
    lonvalues= np.arange(-180.,180.,90.)
    m.drawcountries(); m.drawcoastlines()
    m.drawlsmask(lakes=True)
    m.drawmapboundary() # draw a line around the map region
    m.drawparallels(latvalues,labels=[1, 0, 0, 0])
    m.drawmeridians(lonvalues,labels=[0, 0, 0, 1]) # draw meridians

#=======================================================================

def map_plot(x,use_basemap=False,ax=None,cticks=None,region=None,nclasses=10,cmap_data='jet', title=None,regions_to_plot = None,logplot=False,logoffset=None, **kwargs):
    '''
    produce a plot with
    values for each dataset
    and the difference between the two
    
    
    regions_to_plot (list): list of Region objects. If this argument is given, then the boundaries of each region is plotted in the map together with its label
    '''
    
    #--- create new figure
    if ax == None:
        fig = plt.figure(); ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    
    #--- create colormap
    cmap = plt.cm.get_cmap(cmap_data, nclasses)   
    
    #--- temporal mean fields
    xm = x.timmean()
    
    #--- logscale plot ?
    if logplot:
        if logoffset == None:
            if xm.min() < 0.:
                logoffset = abs(xm.min())*1.01
            else:
                logoffset = 0.
        else:
            logoffset = logoffset
        
        print 'logoffset: ', logoffset
        
        xm = np.log10( xm + logoffset )
    
    #--- set projection parameters
    proj='robin'; lon_0=0.; lat_0=0.
    
    #--- plot using basemap
    if use_basemap:
        llcrnrlon=None; llcrnrlat=None; urcrnrlon=None; urcrnrlat=None
        if region !=None:
            if not hasattr(region,'lonmin'):
                print 'WARNING map boundaries can not be set, as region ' + region.label.upper() + ' has not lat/lon information'
            else:
                dlat = (region.latmax-region.latmin)*0.25; dlon = (region.lonmax-region.lonmin)*0.25
                #~ di = max(dlat,dlon)
                di = 0. #with 0 it works; for other values problems may occur for negative lon!
                llcrnrlon=region.lonmin - di; llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di; urcrnrlat=region.latmax + di
                
                #~ print llcrnrlon, urcrnrlon, llcrnrlat,urcrnrlat
                #~ stop
                
                
                proj='tmerc'
        #generate map
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        
        #use KDTRee nearest neighbor resampling to avoid stripes in plotting
        lons = np.unique(x.lon); lats = np.unique(x.lat)
        lons.sort(); lats.sort()
        TLON,TLAT = np.meshgrid(lons,lats)  #generate target coordinates
        XT,YT = m1(TLON,TLAT)
        X=XT.copy(); Y=YT.copy()
        shape0 = np.shape(XT) 
        XT.shape = (-1); YT.shape = (-1) #... vectorize them for inertpolation
        tree = KDTree(zip(XT,YT)) #generate tree from TARGET coordinates
        
        #prepare data and interpolate
        xmap,ymap = m1(x.lon,x.lat)
        xmap.shape = (-1); ymap.shape = (-1)
        pts  = zip(xmap,ymap) #generate points to interpolate from source data
        dist,idx = tree.query(pts,k=1)     #perform nearest neighbor interpolation (returns distance and indices)
                
        #- map data to output matrix for plotting
        Z = np.ones(shape0)*np.nan; Z.shape = (-1) #generate target vector
        omask = np.ones(shape0).astype('bool'); omask.shape = (-1) 
        
        msk1 = xm.mask.copy(); msk1.shape = (-1); omask[idx] = msk1

        xm1 = xm.copy(); xm1.shape = (-1)
        Z[idx]   = xm1 #assign data and reshape it and set generate masked array

        #~ print Z.shape, omask.shape
        Z[omask] = np.nan
        Z = np.reshape(Z,shape0); Z = np.ma.array(Z,mask=np.isnan(Z))
        
        #~ def _unique_lon(lon):
            #~ ''' check if all longitudes are the same '''
            #~ if lon.ndim == 1:
                #~ return lon.copy()
            #~ elif lon.ndim == 2:
                #~ if (np.any(abs(lon - lon[0,:]) > 0.) )   :
                    #~ raise ValueError, 'No unique longitudes given!'
                #~ else:
                    #~ return lon[0,:].copy()
            #~ else:
                #~ raise ValueError, 'Can not determine unique longitudes!'
            

        #~ if shift: #shift grid
            #~ thelon = x.lon.copy()
            #~ msk = thelon < 0.
            #~ thelon[msk] = thelon[msk] + 360. #0...360
            #~ 
            #~ thelon = _unique_lon(thelon)
            #~ 
            #~ #if masked array then keep mask
            #~ if hasattr(xm,'mask'):
                #~ orgmask = xm.mask.copy()
            #~ else:
                #~ orgmask = np.ones(xm.shape).astype('bool')
            #~ xm,thelon1  = shiftgrid(180.,xm,thelon)
            #~ orgmask,xxx = shiftgrid(180.,orgmask,thelon)
            #~ 
            #~ tmp = x.lat[:,0].copy()
            #~ THELON,TMP = np.meshgrid(thelon1,tmp) #generate 2D field of lon again
            #~ 
            #~ xm = np.ma.array(xm,mask=orgmask)
#~ 
            #~ 
        #~ else:
            #~ THELON = x.lon
        
        #here is still a problem in the plotting over land; masking does not work properly,
        #while the data as such is o.k.!
        #~ im1=m1.pcolormesh(xmap,ymap,xm,cmap=cmap,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        im1=m1.pcolormesh(X,Y,Z,cmap=cmap,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m1)

    else: #use_basemap = False
        #- normal plots
        im1=ax.imshow(xm,cmap=cmap,**kwargs)
        ax.set_xticks([]); ax.set_yticks([])

        
    #set legend aligned with plot (nice looking)
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal("5%", pad=0.05, axes_class=maxes.Axes)
    ax.figure.add_axes(cax) 
    norm = mpl.colors.Normalize(vmin=im1.get_clim()[0], vmax=im1.get_clim()[1])
    cb   = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,ticks=cticks)

    def _add_region(m,r,color='red'):
        '''
        plot region r on top of basemap map m
        '''
        corners = r.get_corners() #get list of corner coordinates
        corners = np.asarray(corners)
        lons = corners[:,0]; lats=corners[:,1]
        x,y = m(lons,lats)
        xy = list(zip(x,y))
        mapboundary = Polygon(xy,edgecolor=color,linewidth=1,fill=False,linestyle='dashed')
        m.ax.add_patch(mapboundary)

    #--- plot regions in the map ---
    if regions_to_plot != None:
        if use_basemap:
            for region in regions_to_plot:
                if region.type=='latlon':
                    _add_region(m1,region)

    
    
    #--- set title
    if title == None:
        ax.set_title(x._get_label(),size=8)
    else:
        ax.set_title(title,size=8)
    
    return fig


#=======================================================================
            
def hov_difference(x,y,climits=None,dlimits=None,**kwargs):
    '''
    
    class to plot hovmoeller diagrams of two datasets
    and their difference
    
    x,y two Data structures
    '''
    
    if climits == None:
        sys.exit('Please specify climits for hovmoeller')
    if dlimits == None:
        sys.exit('Please specify dlimits for hovmoeller')
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    hov1 = hovmoeller(num2date(x.time),x.data,rescaley=6,lat=x.lat)
    hov2 = hovmoeller(num2date(y.time),y.data,rescaley=6,lat=y.lat)
    
    hov1.time_to_lat(**kwargs)
    hov2.time_to_lat(**kwargs)
    
    hov1.plot(title=x._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap='jet',ax=ax1,showcolorbar=False,climits=climits)
    hov2.plot(title=y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap='jet',ax=ax2,showcolorbar=False,climits=climits)
    
    plt.colorbar(hov1.im,ax=ax1,shrink = 0.5,orientation='vertical')
    plt.colorbar(hov2.im,ax=ax2,shrink = 0.5,orientation='vertical')
    
    if x.data.shape == y.data.shape:
        hov3 = hovmoeller(num2date(y.time),x.data - y.data,rescaley=6,lat=y.lat)
        hov3.time_to_lat(**kwargs)
        hov3.plot(title=x._get_label() + ' - ' + y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap='RdBu',ax=ax3,showcolorbar=False,climits=dlimits)
        plt.colorbar(hov3.im,ax=ax3,shrink = 0.5,orientation='vertical')
    else:
        msg = 'Difference plot not possible as data has different shape'
        ax3.text(0.5, 0.5,msg,
             horizontalalignment='center',
             verticalalignment='center') #,
             #transform = ax.transAxes)
        ax3.set_xticks([]); ax3.set_yticks([])
        
    return fig
    
    
    
    
    
def map_difference(x,y,dmin=None,dmax=None,use_basemap=False,ax=None,shift=False,title=None,cticks=None,region=None,nclasses=10,cmap_data='jet',cmap_difference = 'RdBu_r',rmin=-1.,rmax=1., **kwargs):
    '''
    produce a plot with
    values for each dataset
    and the difference between the two
    '''
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    
    #--- get colormap
    cmap = plt.cm.get_cmap(cmap_data, nclasses)   
    
    #- temporal mean fields
    xm = x.timmean(); ym = y.timmean()
    
    proj='robin'; lon_0=0.; lat_0=0.
    
    #- plot first dataset
    map_plot(x,use_basemap=use_basemap,ax=ax1,cticks=cticks,region=region,nclasses=nclasses,cmap_data=cmap_data, title=title, shift=shift, **kwargs)
    
    #- plot second dataset
    map_plot(y,use_basemap=use_basemap,ax=ax2,cticks=cticks,region=region,nclasses=nclasses,cmap_data=cmap_data, title=title, shift=shift, **kwargs)
    
    #-first minus second dataset
    map_plot(x.sub(y),use_basemap=use_basemap,ax=ax3,vmin=dmin,vmax=dmax,cticks=None,region=region,nclasses=nclasses,cmap_data=cmap_difference, title='absolute difference [' + x.unit + ']', shift=shift)
    
    #- relative error
    map_plot(y.div(x).subc(1.),use_basemap=use_basemap,ax=ax4,vmin=rmin,vmax=rmax,title='relative difference',cticks=None,region=region ,nclasses=nclasses,cmap_data=cmap_difference, shift=shift)
    

    return fig
