#!/usr/bin/pythonl
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

#TODO
#- analyze only dry/wet years




from matplotlib import pylab as plt

from mpl_toolkits.basemap import Basemap,shiftgrid

from scipy import stats

import numpy as np

from matplotlib.patches import Circle

import sys

from diagnostic import *


#todo:
# - implement Glecker plots
#- implement writing of statistics to an ASCII file as export


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
        
    def bar(self,vmin=None,vmax=None,**kwargs):
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
        self.ax.set_ylabel('relative error to multimodel mean [%]')
        self.ax.grid()
        
        
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
        self.regress = regress
        self.title = title
        self.show_xlabel = show_xlabel
        self.show_ylabel = show_ylabel
        
        
    def plot(self,x,ax=None,vmin=None,vmax=None,**kwargs):
        '''
        x: object of Data class
        '''
        
        if len(x.time) > 0:
            
            if ax == None:
                ax = self.ax
            else:
                ax = ax
            
            y = x.fldmean()
            label = x.label
            
            if self.regress: #calculate linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(x.time,y)
                label = label + ' (r=' + str(round(r_value,2)) + ', p=' + str(round(p_value,2)) + ')'

            p = ax.plot(plt.num2date(x.time), y , label=label,**kwargs)[0]
            if self.regress:
                ax.plot(x.time,x.time*slope+intercept,'--',color=p.get_color()) #plot regression line

            if self.show_ylabel:
                ax.set_ylabel(x.unit)
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
            
        #~ print dat
        #~ print dat.shape
        #~ print x.lat.shape
        #~ stop
        
        #~ plt.figure()
        #~ plt.plot(dat[0,:])
        #~ plt.plot(dat[2,:])
        
        if dat.shape[x.data.ndim-2] != x.lat.shape[0]:
            print 'Inconsistent shapes!'
            print dat.shape
            print x.lat.shape
            sys.exit()
        
        #~ print dat
        
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
        
        #~ fig = self.ax.figure
        #~ fig.legend()
        
        
        
     
        



#=======================================================================

def __basemap_ancillary(m):
    latvalues=np.arange(-90.,120.,30.)
    lonvalues= np.arange(-200.,200.,30.)
    m.drawcountries(); m.drawcoastlines()
    m.drawlsmask(lakes=True)
    m.drawmapboundary() # draw a line around the map region
    m.drawparallels(latvalues,labels=[1, 0, 0, 0])
    m.drawmeridians(lonvalues,labels=[0, 0, 0, 1]) # draw meridians

#=======================================================================

def map_plot(x,use_basemap=False,ax=None,cticks=None,region=None,nclasses=10,cmap_data='jet', title=None, **kwargs):
    '''
    produce a plot with
    values for each dataset
    and the difference between the two
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
                di = max(dlat,dlon)
                di = 0. #with 0 it works; for other values problems may occur for negative lon!
                llcrnrlon=region.lonmin - di; llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di; urcrnrlat=region.latmax + di
                
                #~ print llcrnrlon, urcrnrlon, llcrnrlat,urcrnrlat
                #~ stop
                
                
                proj='tmerc'
        #generate map
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        xmap, ymap = m1(x.lon,x.lat)
        
        #~ plt.figure(); plt.imshow(xmap);
        #~ plt.figure(); plt.imshow(ymap);
        #~ stop
        
        #~ plt.figure()
        #~ plt.imshow(xm)
        #~ stop
        
        
        im1=m1.pcolormesh(xmap,ymap,xm,cmap=cmap,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m1)
        #~ plt.colorbar(im1,ax=ax,ticks=cticks)

    else: #use_basemap = False
        #- normal plots
        im1=ax.imshow(xm,cmap=cmap,**kwargs); plt.colorbar(im1,ax=ax,ticks=cticks)
        ax.set_xticks([]); ax.set_yticks([])

    #--- set title
    if title == None:
        ax.set_title(x._get_label())
    else:
        ax.set_title(title)
    
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
    
    
    
    
    
def map_difference(x,y,dmin=None,dmax=None,use_basemap=False,ax=None,cticks=None,region=None,nclasses=10,cmap_data='jet',cmap_difference = 'RdBu_r', **kwargs):
    '''
    produce a plot with
    values for each dataset
    and the difference between the two
    '''
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    #--- get colormap
    cmap = plt.cm.get_cmap(cmap_data, nclasses)   
    
    #- temporal mean fields
    xm = x.timmean(); ym = y.timmean()
    
    proj='robin'; lon_0=0.; lat_0=0.
    
    if use_basemap:
        #draw_map(x.lat,x.lon,xm,'none.png',vmin=-2.,vmax=2.,show_grids=False,show_colorbar=True,tit=x.label,resolution='l',proj=proj,save=False,outdpi=300.,ax=ax1,ccmap=plt.cm.jet)
        
        #set map boundaries by region
        llcrnrlon=None; llcrnrlat=None; urcrnrlon=None; urcrnrlat=None
        if region !=None:
            if not hasattr(region,'lonmin'):
                print 'WARNING map boundaries can not be set, as region ' + region.label.upper() + ' has not lat/lon information'
            else:
                dlat = (region.latmax-region.latmin)*0.25
                dlon = (region.lonmax-region.lonmin)*0.25
                di = max(dlat,dlon)
                llcrnrlon=region.lonmin - di
                llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di
                urcrnrlat=region.latmax + di
                
                proj='tmerc'

            
        
        

            
            
        
        #plot 1
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax1,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        xmap, ymap = m1(x.lon,x.lat)
        im1=m1.pcolormesh(xmap,ymap,xm,cmap=cmap**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m1)
        plt.colorbar(im1,ax=ax1,ticks=cticks)
        
        #plot 2
        m2=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax2,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        xmap, ymap = m2(y.lon,y.lat)
        im2=m2.pcolormesh(xmap,ymap,ym,cmap=cmap**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m2)
        plt.colorbar(im2,ax=ax2,ticks=cticks)
        
        if xm.shape == ym.shape:
            #plot 2
            m3=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax3,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            xmap, ymap = m3(y.lon,y.lat)
            cmap_diff = plt.cm.get_cmap(cmap_difference, nclasses)    # discrete colors
            #im3=m3.pcolormesh(xmap,ymap,xm-ym,vmin=dmin,vmax=dmax,cmap=plt.cm.RdBu) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
            im3=m3.pcolormesh(xmap,ymap,xm-ym,vmin=dmin,vmax=dmax,cmap=cmap_diff) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
            __basemap_ancillary(m3)
            plt.colorbar(im3,ax=ax3)
            
            
        else:
            msg = 'Difference plot not possible as data has different shape'
            ax3.text(0.5, 0.5,msg,
                 horizontalalignment='center',
                 verticalalignment='center') #,

        
    else: #use_basemap
        #- normal plots
        cmap_diff = plt.cm.get_cmap(cmap_difference, nclasses)    # discrete colors
        im1=ax1.imshow(xm,cmap=cmap,**kwargs); plt.colorbar(im1,ax=ax1,ticks=cticks)
        im2=ax2.imshow(ym,cmap=cmap,**kwargs); plt.colorbar(im2,ax=ax2,ticks=cticks)
        
        if xm.shape == ym.shape:
            if (dmin != None) & (dmax != None):
                im3=ax3.imshow(xm - ym,cmap=cmap_diff,vmin=dmin,vmax=dmax); plt.colorbar(im3,ax=ax3)
            else:
                im3=ax3.imshow(xm - ym,cmap=cmap_diff); plt.colorbar(im3,ax=ax3)
            ax3.set_title(x._get_label() + ' - ' + y._get_label() )
        else:
            msg = 'Difference plot not possible as data has different shape'
            ax3.text(0.5, 0.5,msg,
                 horizontalalignment='center',
                 verticalalignment='center') #,
                 #transform = ax.transAxes)
        
        ax1.set_xticks([]); ax1.set_yticks([])
        ax2.set_xticks([]); ax2.set_yticks([])
        ax3.set_xticks([]); ax3.set_yticks([])
        
    ax1.set_title(x._get_label())
    ax2.set_title(y._get_label())
    
    return fig
