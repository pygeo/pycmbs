#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"


from matplotlib import pylab as plt

from mpl_toolkits.basemap import Basemap,shiftgrid

from scipy import stats

import numpy as np

from matplotlib.patches import Circle

#~ from python.hov import *


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
        
        
    def circle_plot(self):
        print 'Doing Reichler plot as circle plot ...'
        self._normalize()
        
        dx=0.
        tsize=10.
        for i in np.arange(len(self.e2)):
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
    def __init__(self,ax):
        self.ax = ax
        
        
    def plot(self,x,**kwargs):
        '''
        x: object of Data class
        '''
        
        if len(x.time) > 0:
            self.ax.plot(plt.num2date(x.time), x.fldmean(), label=x.label,**kwargs)
            
            
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
    
    
    
    
    
def map_difference(x,y,dmin=None,dmax=None,use_basemap=False,ax=None,cticks=None,region=None,nclasses=10,**kwargs):
    '''
    produce a plot with
    values for each dataset
    and the difference between the two
    '''
    
    fig = plt.figure()
    
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
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
                di = max(dlat,dmin)
                llcrnrlon=region.lonmin - di
                llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di
                urcrnrlat=region.latmax + di
                
                proj='tmerc'

            
        
        
        def __basemap_ancillary(m):
            latvalues=np.arange(-90.,120.,10.)
            lonvalues= np.arange(-200.,200.,10.)
            m.drawcountries(); m.drawcoastlines()
            m.drawlsmask(lakes=True)
            m.drawmapboundary() # draw a line around the map region
            m.drawparallels(latvalues,labels=[1, 0, 0, 0])
            m.drawmeridians(lonvalues,labels=[0, 0, 0, 1]) # draw meridians
            
        #plot 1
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax1,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        xmap, ymap = m1(x.lon,x.lat)
        im1=m1.pcolormesh(xmap,ymap,xm,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m1)
        plt.colorbar(im1,ax=ax1,ticks=cticks)
        
        #plot 2
        m2=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax2,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        xmap, ymap = m2(y.lon,y.lat)
        im2=m2.pcolormesh(xmap,ymap,ym,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m2)
        plt.colorbar(im2,ax=ax2,ticks=cticks)
        
        if xm.shape == ym.shape:
            #plot 2
            m3=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax3,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
            xmap, ymap = m3(y.lon,y.lat)
            cmap = plt.cm.get_cmap('RdBu_r', nclasses)    # discrete colors
            #im3=m3.pcolormesh(xmap,ymap,xm-ym,vmin=dmin,vmax=dmax,cmap=plt.cm.RdBu) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
            im3=m3.pcolormesh(xmap,ymap,xm-ym,vmin=dmin,vmax=dmax,cmap=cmap) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
            __basemap_ancillary(m3)
            plt.colorbar(im3,ax=ax3)
            
            
        else:
            msg = 'Difference plot not possible as data has different shape'
            ax3.text(0.5, 0.5,msg,
                 horizontalalignment='center',
                 verticalalignment='center') #,

        
    else: #use_basemap
        #- normal plots
        cmap = plt.cm.get_cmap('RdBu_r', nclasses)    # discrete colors
        im1=ax1.imshow(xm,**kwargs); plt.colorbar(im1,ax=ax1,ticks=cticks)
        im2=ax2.imshow(ym,**kwargs); plt.colorbar(im2,ax=ax2,ticks=cticks)
        
        if xm.shape == ym.shape:
            if (dmin != None) & (dmax != None):
                im3=ax3.imshow(xm - ym,cmap=cmap,vmin=dmin,vmax=dmax); plt.colorbar(im3,ax=ax3)
            else:
                im3=ax3.imshow(xm - ym,cmap=cmap); plt.colorbar(im3,ax=ax3)
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
