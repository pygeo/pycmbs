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
Module that contains relevant classes for diagnostic plots

@todo: implement writing of statistics to an ASCII file as export
@todo: implement taylor plots
@todo: faster implementation of Basemap plots. For large number of grid cells, the current KTree implementation is by far too slow!

'''

from data import *

from hov import *

from matplotlib import pylab as plt

from matplotlib.patches import Polygon

from mpl_toolkits.basemap import Basemap,shiftgrid

from scipy import stats

import numpy as np

from matplotlib.patches import Circle

import sys


from scipy.spatial import cKDTree as KDTree #import the C version of KDTree (faster)

from matplotlib.ticker import MaxNLocator

import matplotlib.gridspec as gridspec

#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes

import matplotlib as mpl

#-----------------------------------------------------------------------

def thin_xticks(ax,n):
    """
    thin xticks of axis

    If there are too many xticks in a plot or the labels
    are overlapping, it makes sense to thin the mńumber of labels

    @param ax: axis that will be treated
    @type ax: matplotlib axis

    @param n: number of ticks to plot
    @type n: int
    """
    ax.xaxis.set_major_locator(MaxNLocator(n+1))

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class CorrelationAnalysis():
        """
        perform correlation analysis between two datasets
        and plot results
        """

        def __init__(self,X,Y,mask=None,ax=None):
            """
            constructor of class

            @param X: x dataset (either [time,sample] or [time,sample,sample]
            @type X: numpy array

            @param Y: y dataset (either [time,sample] or [time,sample,sample]
            @type Y: numpy array

            @param mask: mask to be applied to the data
            @type mask: numpy array(:,:) or (:)

            @param ax: axis to plot results to; new figure will be generated if ax==None
            @type ax: matplotlib axis
            """

            self.x = X; self.y = Y
            self.mask = mask

            if ax == None:
                f = plt.figure()
                self.ax = f.add_subplot(111)
            else:
                self.ax = ax

#-----------------------------------------------------------------------

        def do_analysis(self):
            """
            perform correlation analysis

            @todo: implement area weighting
            @todo: implement regional (condition) statisitcs based on a mask
            @todo: return value
            """

            #--- calculate diagnostics
            D = Diagnostic(self.x,y=self.y)
            D._mat2vec(mask = self.mask) #here is the point fo regional statistics
            rmse = D.get_rmse_value()
            r    = D.get_correlation_value()
            n    = D. get_n()

            print 'RMSE: ', rmse
            print 'R   : ', r
            print 'N   : ', n

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class HovmoellerPlot():
    def __init__(self,D,rescaley=10,rescalex=10,yticksampling=1,monthly=False,ax=None):
        """
        D : C{Data} object

        if the argument lat is provided it is assumed that lat/lon are 2D matrices
        In this case the value is expected to be a 3D variables as
        value(time,ny,nx)
        """
        if ax == None:
            self.figure = pl.figure()
            self.ax = self.figure.add_subplot(111)
        else:
            self.figure = self.ax.figure

        self.hov = hovmoeller(pl.num2date(D.time),None,rescaley=rescaley,rescalex=rescalex)
        #self.hov.time_to_lat(dlat=dlat,yticksampling=yticksampling,monthly=monthly)
        self.x = D

    def plot(self,title=None,climits=None,showxticks=True,showcolorbar=True,cmap='jet',xtickrotation=90,ylim=None):
        if climits is None:
            raise ValueError, 'CLIMITS needs to be specified!'
        self.hov.plot(input=self.x,ax=self.ax,title=title,ylabel='lat',xlabel='days',origin='lower',xtickrotation=xtickrotation,climits=climits,showxticks=showxticks,showcolorbar=showcolorbar,cmap=cmap)
        if ylim != None:
            self.ax.set_ylim(ylim)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class ReichlerPlot():
    """
    class for Reichler plot generation

    @todo: add example how to use Reichler plotting
    @todo: provide references Glecher and Reichler + Kim
    """
    def __init__(self,ax=None):
        """
        constructor for Reichler plot

        @param ax: axis to plot data to; if None, new figure will be generated
        @type ax: matplotlib axis
        """
        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

        self.e2 = [] #list to store RMS error results
        self.labels = []; self.colors=[]

#-----------------------------------------------------------------------

    def add(self,e2,label,color=None):
        """
        register data to be plotted

        @param e2: reichler index that was already calculated
        @type e2: list

        @param label: label to be used for plotting
        @type label: str

        @param color: color to be used for plotting
        @type color: str
        """
        self.e2.append(e2)
        self.labels.append(label)
        self.colors.append(color)

#-----------------------------------------------------------------------

    def bar(self,vmin=None,vmax=None,title='',**kwargs):
        """
        generate barplot which shows results from all diagnostic
        values (e.g. different model results

        it calculates the mean error of all model and then plots
        the relative error of each particular model
        compared to the multimodel mean

        @param vmin: minimum value for plotting
        @type vmin: int

        @param vmax: maximum value for plotting
        @type vmax: int

        @param title: title for the plot
        @type title: str
        """

        if len(self.e2) == 0: #no valid data
            return self.ax.figure

        #- normalize results (relative modle performance)
        self._normalize()
        x = np.arange(len(self.e2_norm))
        y1 = self.e2_norm*100.; y2 = self.e2_norm*100.
        y1[y1 < 0.] = np.nan  #posistive values only
        y2[y2 >= 0.] = np.nan #negative values only
        self.ax.bar(x,y1,color='red' ,edgecolor='None',**kwargs)
        self.ax.bar(x,y2,color='blue',edgecolor='None',**kwargs)

        self.ax.set_xticks(x+0.5)
        self.ax.set_xticklabels(self.labels,rotation=90.)
        if (vmin !=None) & (vmax != None):
            self.ax.set_ylim(vmin,vmax)
        else:
            vmin,vmax = self.ax.get_ylim() #equal axes
            if abs(vmax) > abs(vmin):
                vmin = -vmax
                vmax = vmax
            else:
                vmin = vmin
                vmax = -vmin
            self.ax.set_ylim(vmin,vmax)
        self.ax.set_ylabel('$\\epsilon / \\bar{\\epsilon}$ [%]')
        self.ax.grid(); self.ax.set_title(title)

        return self.ax.figure


#-----------------------------------------------------------------------

    def simple_plot(self):
        '''
        do a very simple plot of diagnostics
        '''
        for i in np.arange(len(self.e2)):
            self.ax.plot(self.e2[i],'o',label=self.labels[i])

#-----------------------------------------------------------------------

    def circle_plot(self):
        '''
        nice looking plot of Reichler index
        '''
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


        self.ax.set_ylim(-1.,1.); self.ax.set_xlim(-1.,1.)
        self.ax.set_xlabel('$\\epsilon / \\bar{\\epsilon}$ [%]')
        self.ax.legend()

#-----------------------------------------------------------------------

    def _normalize(self):
        """
        normalize results from different models
        Glecker et al, eq. 2
        """

        n  = len(self.e2[0])
        E2 = []

        for e in self.e2:
            if len(e) != n:
                print 'WARNING: non consistent length in error statistics!!!'
            E2.append(np.nansum(e)) #temporal aggregation

        E2 = np.asarray(E2);  EM = E2.mean()
        self.e2_norm =  (E2 - EM) / EM #see Glecker et al, eq.2

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class ScatterPlot():
    '''
    Class for generation of scatterplots
    '''
    def __init__(self,x,ax=None,ticksize=10,normalize_data=False,show_xlabel=True):
        '''
        constructor of class C{ScatterPlot}

        @param x: Variable that will be used as the x-variable
        @type x: C{Data} object

        @param normalize_data: if True, then the dataseries is normalizued internally so that it has zero mean and a std of 1
        @type normalize_data: bool

        @param show_xlabel: show xlabel in plot
        @type show_xlabel: bool
        '''

        self.show_xlabel = show_xlabel


        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

        self.figure = self.ax.figure
        self.x = x
        self.lines = []; self.labels = []
        self.ticksize=ticksize

        self.normalize = normalize_data


#-----------------------------------------------------------------------
    def __normalize_data(self,x):
        '''
        normmalize timeseries
        '''
        return (x-x.mean()) / x.std()


#-----------------------------------------------------------------------

    def plot(self,y,regress=True,**kwargs):
        '''
        add a dataset to the scatterplot and plot
        it. It also allows to perform a regression analysis

        @param y: data to be plotted on y-axis
        @type y: C{Data} object

        @param regress: Perform linear regression analysis
        @type regress: bool
        '''
        label=y.label
        xdat = self.x.fldmean(); ydat = y.fldmean()

        if self.normalize:
            xdat = self.__normalize_data(xdat)
            ydat = self.__normalize_data(ydat)

        #- calculate linear regression
        #~ print xdat.shape
        if regress:
            slope, intercept, r_value, p_value, std_err = stats.linregress(xdat,ydat)
            if p_value < 0.01:
                spvalue = 'p < 0.01'
            else:
                spvalue = 'p=' + str(round(p_value,2))
            label = label + ' (r=' + str(round(r_value,2)) + ', ' + spvalue + ')'
            rms_error = np.mean(((xdat-ydat)**2))
            std_error = np.std(xdat-ydat)

        l = self.ax.plot(xdat,ydat,'.',label=label,**kwargs)[0]
        if regress:
            self.ax.plot(xdat,xdat*slope+intercept,'--',color=l.get_color())
        self.lines.append(l); self.labels.append(label)

        if self.show_xlabel:
            self.ax.set_xlabel(self.x._get_label(),size=self.ticksize )
        self.ax.set_ylabel(y._get_unit(),size=self.ticksize)

        self._change_ticklabels()

        if regress:
            return r_value,p_value,rms_error, std_error
        else:
            return None

    def _change_ticklabels(self):
        for tick in self.ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(self.ticksize)
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)

#-----------------------------------------------------------------------

    def legend(self):
        """
        plot legend
        """
        self.ax.legend(self.lines,self.labels,prop={'size':8})

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class LinePlot():
    """
    class for a pyCMBS Line Plot

    This class is usefull for plotting timeseries
    """
    def __init__(self,ax=None,regress=False,title=None,show_xlabel=True,show_ylabel=True,ticksize=10,normx=1.,show_equation=True,xtickrotation=90):
        """
        constructor of LinePlot

        @param ax: axis to plot data to. If I{None} then a new figure is generated
        @type ax: matplotlib axis

        @param regress: perform a linear regression on the data to be plotted
        @type regress: bool

        @param title: title of the plot
        @type title: str

        @param show_xlabel: show x-label for the plot
        @type show_xlabel: bool

        @param show_ylabel: show y-label for the plot
        @type show_ylabel: bool

        @param normx: normalization constant for x-variable (needed e.g. if you want to normalize a timevector for regression analysis)
        @type normx: float

        @param xtickrotation: rotation for xtick labels
        @type xtickrotation: float

        """

        if ax == None:
            f = plt.figure(figsize=(8,7))
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

        self.figure = self.ax.figure
        self.regress = regress
        self.title = title
        self.show_xlabel = show_xlabel
        self.show_ylabel = show_ylabel

        self.lines = []; self.labels = []

        self.ticksize = ticksize

        self.normx=normx
        self.show_equation = show_equation

        self.xtickrotation = xtickrotation

        #~ self.showxticks = showxticks

#-----------------------------------------------------------------------

    def legend(self,prop={'size':8},**kwargs):
        """
        plot legend
        """
        self.ax.legend(self.lines,self.labels,prop=prop,**kwargs)

#-----------------------------------------------------------------------

    def _change_ticklabels(self,ax=None):
        if ax == None:
            ax = self.ax

        #yaxis
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)
        #xaxis
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)
            tick.label.set_rotation(self.xtickrotation)

#-----------------------------------------------------------------------

    def plot(self,x,ax=None,vmin=None,vmax=None,label = None, norm_std = False, set_ytickcolor=True, **kwargs):
        """
        plot LinePlot data. If a spatial field is provided, this is aggregated
        using the fldmean() function of C{Data}

        @param x: data to be plotted
        @type x: C{Data}

        @param ax: axis to plot to. If None, then a new figure is generated
        @type ax: matplotlib axis

        @param vmin: minimum value for y-axis
        @type vmin: float

        @param vmax: maximum value for y-axis
        @type vmax: float

        @param label: label to be used for current plot. If None, then
                      the label of the provided C{Data} object is used
        @type label: str

        @param norm_std: normalize timeseries with its stdv. This is a useful option when comparing trends of variables with different amplitudes
        @type norm_std: bool
        """

        if len(x.time) > 0:

            if ax == None:
                ax = self.ax
                set_axiscolor=False
            else:
                ax = ax
                set_axiscolor=True

            y = x.fldmean() #gives timeseries

            if norm_std:
                y /= y.std()


            if label == None:
                label = x.label

            if self.regress: #calculate linear regression
                slope_print, intercept_print, r_value, p_value, std_err = stats.linregress(x.time/self.normx,y)
                slope, intercept, r_value, p_value, std_err = stats.linregress(x.time,y)
                self.tmp_slope = slope
                self.tmp_corr = r_value

                if p_value < 0.01:
                    spvalue = 'p < 0.01'
                else:
                    spvalue = 'p=' + str(round(p_value,2))

                if self.show_equation:
                    label = label + ' (y=' + "%.1e" % slope_print + 'x+' + "%.1e" % intercept_print  + ', r=' + str(round(r_value,2)) + ', ' + spvalue + ')'
                else:
                    label = label + ' (r=' + str(round(r_value,2)) + ', ' + spvalue + ')'

            self.labels.append(label)

            p = ax.plot(plt.num2date(x.time), y , label=label, **kwargs)[0]
            self.lines.append(p)
            if self.regress:
                ax.plot(x.time,x.time*slope+intercept,'--',color=p.get_color()) #plot regression line

            if self.show_ylabel:
                ax.set_ylabel(x._get_unit(),size=self.ticksize)
            if self.show_xlabel:
                ax.set_xlabel('time',size=self.ticksize)

            if self.title != None:
                ax.set_title(self.title,size=10)

            if vmin != None and vmax != None:
                ax.set_ylim(vmin,vmax)

            if set_ytickcolor:
                for tl in ax.get_yticklabels():
                    tl.set_color(p.get_color())


            self._change_ticklabels(ax)

            #~ if not self.showxticks: #showxticklabels?
                #~ for tick in ax.xaxis.get_major_ticks():
                    #~ print 'tick'
                    #~ tick.label.set_label('')

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class GlobalMeanPlot():
    """
    plots timeseries of global mean field
    """

    def __init__(self,ax=None,climatology=True,ax1=None):
        """
        @param ax: specifies axis to plot the data to
        @type ax: axis

        @param climatology: specifies if a second plot for a climatological mean value shall be generated
        @type climatology: bool
        """
        if climatology:
            nplots = 2
        else:
            nplots = 1
        self.climatology = climatology

        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(nplots,1,1)
            if self.climatology:
                if ax1 == None:
                    raise ValueError, 'For climatology plot, we need a ax1 as an argument!'
                else:
                    self.ax1 = ax1
        else:
            if self.climatology:
                self.ax  = ax.figure.add_subplot(nplots,1,1)
                self.ax1 = ax.figure.add_subplot(nplots,1,2)
            else:
                self.ax = ax

        self.labels=[]; self.plots=[]

    def plot(self,D1,color=None,linewidth=1,show_std=False,label=None,linestyle='-',mask=None):
        """
        generate global mean plot. The plot includes the temporal evolution
        of the global mean field and also (as an option) its stdv

        assumes data structure of [time,ny,nx]

        labels need to be unique and are derived either from the optional
        argument or, if a C{Data} object is given they are derived from
        the Data label. In case of a duplication of the data labels,
        no plot will be done!

        @param: D1 data field to be plotted
        @type: C{Data} or (time,data) tuple

        @param color: color of the line
        @type color: str

        @param linewidth: width of the line
        @type linewidth: float

        @param show_std: shows standard deviation
        @type show_std: bool

        @param mask: mask to be applied to the data prior to final analyis
        @type mask: either numpy bool array or C{Data} object
        """

        if ((label==None) and (D1.label in self.labels)):
            #print 'Label already existing: ', D.label, ' skipping analysis'
            return
        elif ((label != None) and (label in self.labels)):
            #print 'Label already existing: ', label, ' skipping analysis'
            return

        #ensure to work with a data object
        if 'tuple' in str(type(D1)): #only a vector is provided as data together with time (time,data)
            D = D1[2] #(time,data,orgdata)
        else:
            D = D1

        if D.data.ndim != 3:
            raise ValueError, 'Global plot only supported for 3D data'

        if mask != None:
            D._apply_mask(mask)

        #mean field
        m = D.fldmean(return_data=True) #mean
        mdata = m.data.flatten()

        #time
        t = plt.num2date(m.time)

        #--- plot generation ---
        if color == None:
            p = self.ax.plot(t,mdata,linewidth=linewidth,linestyle=linestyle)
        else:
            p = self.ax.plot(t,mdata,color=color,linewidth=linewidth,linestyle=linestyle)

        if show_std:
            s = D.fldstd (return_data=True) #std
            sdata = s.data.flatten()
            self.ax.fill_between(t,mdata-sdata,y2=mdata+sdata,color=p[0].get_color(),alpha=0.5)

        #- plot climatology if desired
        if self.climatology:
            tmp = D.get_climatology(return_object=True)
            tmp.adjust_time(year=1700,day=15)
            tmp.timsort()
            m = tmp.fldmean(return_data=False).flatten()

            self.ax1.plot(np.arange(1,13),m,linestyle=linestyle)
            self.ax1.set_xlim(0.,13.)
            self.ax1.set_ylabel(tmp._get_unit())
            del tmp

            self.ax1.set_xlabel('months')

        #- store information for legend
        self.plots.append(p[0])
        if label==None:
            self.labels.append(D.label)
        else:
            self.labels.append(label)

        #- labels
        self.ax.set_ylabel(D._get_unit())
        self.ax.set_xlabel('time')

        #- legend
        self.ax.legend(self.plots,self.labels,loc='lower center',ncol=2,fancybox=True)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class HistogrammPlot():
    """
    class to plot histograms based on C{Data} objects
    """
    def __init__(self,ax=None,bins=10):
        """
        @param ax: axis to plot data to. If not specified, then a new figure is created
        @type: ax: axis

        @param bins: bins for histogram calculation, either int or a list
        @type bins: int or list or array
        """

        #- Figure init
        if ax == None:
            self.figure = pl.figure()
            self.ax = self.figure.add_subplot(111)
        else:
            self.ax = ax
            self.figure = self.ax.figure

        self.bins = bins

    def plot(self,X,color='black',linestyle='-',linewidth=1.,label='',shown=False,show_legend=False,**kwargs):
        """
        plot data to histogram

        @param X: data to be plotted as histogram
        @type X: Data or np.array

        @param color: color for line plot
        @type color: str

        @param linestyle: style of line to plot
        @type linestyle: str

        @param linewidth: width of line to plot
        @type linewidth: float

        @param shown: show number of used datasets in legend
        @type shown: bool

        @param kwargs: arguments for np.histogram function
        """

        #-check if Data object
        if isinstance(X,Data):
            x = X.data
        else:
            x = X

        if isinstance(x,np.ma.masked_array):
            x = x.data[~x.mask]

        x = x[~np.isnan(x)]

        #- REMOVE BINS ARGUMENT IF in kwargs, as global argument of class is used
        if 'bins' in kwargs.keys():
            bb = kwargs.pop('bins')

        if shown:
            show_legend=True

            if label == '':
                label = 'n='+str(sum(~np.isnan(x)))
            else:
                label= label + '(n='+str(sum(~np.isnan(x))) + ')'

        f,b = np.histogram(x,bins = self.bins,**kwargs)
        self.ax.plot(b[0:-1],f,color=color,linestyle=linestyle,linewidth=linewidth,label=label)
        if show_legend:
            self.ax.legend()


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


class ZonalPlot():
    def __init__(self,ax=None,dir='y'):
        '''
        @todo: still needs to take into account appropriately area weighting

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

#-----------------------------------------------------------------------

    def plot(self,x,xlim=None,timmean = False,show_ylabel=True):
        """
        plot zonal plot

        @param x: data to be plotted
        @type x: C{Data} object

        @param xlim: limits for the x-axis (e.g. values)
        @type xlim: tuple

        @param timmean: temporal mean calculation
        @type timmean: bool

        """

        #check if all latitudes are the same
        lu = x.lat.mean(axis=1)
        if any( abs(lu - x.lat[:,0]) > 1.E-5):
            print 'WARNING: latitudes are not unique!!!'
            print lu.shape,x.lat.shape
            print lu

            print x.lat[:,0]
            print x.lat[:,0] - lu

            raise ValueError, 'Cannot work with NOT unique latitude values!'

        if timmean:
            thex = x.timmean(return_object=True)
        else:
            thex = x

        if self.dir == 'y':
            dat = thex.get_zonal_mean() #no area weighting performed
        else:
            raise ValueError, 'Invalid option'

        #return dat
        if timmean:
            pass
        else:
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
                self.ax.plot(dat[i,:],x.lat[:,0],label='time='+str(i))
                self.ax.grid(b='on')

        self.ax.set_ylim(-90.,90.)

        if show_ylabel:
            self.ax.set_ylabel('latitude [deg]')
        else:
            self.ax.set_yticks([])

        if xlim != None:
            self.ax.set_xlim(xlim)

        self.ax.grid(b='on')


#-----------------------------------------------------------------------

class GlecklerPlot():
    """
    Class to generate a plot that to illustrate multi-model, multi-variable scores

    It was introdcued by Gleckler et al (2008)

    REFERENCES:
    * ﻿Gleckler, P.J., Taylor, K.E. & Doutriaux, C., 2008. Performance metrics for climate models. Journal of Geophysical Research, 113(D6). Available at: http://www.agu.org/pubs/crossref/2008/2007JD008972.shtml [Accessed February 29, 2012].

    EXAMPLE:
    G = GlecklerPlot()
    #register first models
    G.add_model('echam5'); G.add_model('mpi-esm')
    #then register variables
    G.add_variable('ta'); G.add_variable('P')
    #after that you can add values to be plotted; pos=1 mean that result is plotted in upper triangle
    G.add_data('ta','echam5',0.5,pos=1)
    G.add_data('P','echam5',0.25,pos=1)
    G.add_data('P','echam5',-0.25,pos=2)
    G.add_data('P','mpi-esm',-0.25,pos=1)
    G.plot() #do plot
    """

    def __init__(self,fig=None):
        """
        constructor of C{GlecklerPlot}

        @param fig: figure to which to plot to. If None, then a new figure will be generated
        @type fig: matplotlib figure
        """
        if fig == None:
            color='grey'
            fig = plt.figure(facecolor=color,edgecolor=color)
        self.fig = fig

        self.models = []
        self.variables   = []
        self.data = {} #store data for plot
        self.pos = {} #store position of plot

    def add_model(self,label):
        """
        register a model in the class
        @param label: string of the model
        @type label: str
        """
        s = label.replace(' ','_')
        if s not in self.models:
            self.models.append(s)

    def add_variable(self,label):
        """
        register a variable in the class
        @param label: string of variable
        @type label: str
        """
        self.variables.append(label)

    def __set_ax_prop(self,ax):
        """
        set axis properties of a subplot
        @param ax: subplot axis
        @type ax: matplotlib axis
        """
        ax.set_xticks([]); ax.set_yticks([])

    def __value2color(self,v):
        """
        return a color based on the value given
        the information on the colormap and its
        normalization is used for that purpose

        @param v: value of data
        @type v: float
        """
        return self.cmap(self.norm(v))

    def __plot_triangle(self,ax,value,pos='top'):
        """
        Plot a triangle and fill its color in accordance
        with the value given. Information on colormap
        will be obtained from class information

        @param ax: axis to plot to
        @type ax:  matplotlib axis object

        @param value: value to plot
        @type value: float

        @param pos: position of the triangle which will be plotted
                    top = plot an upper triangle, else = plot a lower triangle
        @type pos: str
        """
        if value == None:
            return

        color = self.__value2color(value)

        pmax = max(self.pos.values())

        if pmax > 4:
            raise ValueError, 'Only up to 4 observations supported!'

        if pmax > 2:
            #plot 4 triangles
            if pos == 'top':
                x = [0.,1.,0.5]
                y = [1.,1.,0.5]
            elif pos == 'bottom':
                x = [0.,0.5,1.]
                y = [0.,0.5,0.]
            elif pos == 'left':
                x = [0.,0.,0.5]
                y = [0.,1.,0.5]
            elif pos == 'right':
                x = [1.,0.5,1.]
                y = [0.,0.5,1.]
            else:
                raise ValueError, 'Invalid position for plot'

        else:
            #- plot only two triangles (diagonal)
            if pos == 'top':
                x = [0.,0.,1.]
                y = [0.,1.,1.]
            elif pos == 'bottom':
                x = [1.,1.,0.]
                y = [1.,0.,0.]
            else:
                raise ValueError, 'Invalid position for plot'

        xy = list(zip(x,y))
        p = Polygon(xy,edgecolor='white',linewidth=1,fill=True,linestyle='solid',facecolor=color)
        ax.add_patch(p)

#-----------------------------------------------------------------------

    def _normalize_data(self):
        """
        calculate for each observational data set
        the relative deviation from the average
        """
        pos = np.unique(self.pos.values())
        if hasattr(self,'_raw_data'):
            #data has been already normalized; take original data
            self.data = self._raw_data.copy()
        else:
            self._raw_data = self.data.copy() #preserve original calculated data

        for var in self.variables:
            for p in pos:
                xm = self._get_mean_value(p,var) #calculate multimodel mean
                for k in self.data:
                    if (self.pos[k] == p) & ('_' + var + '_' in k):
                        self.data[k] = (self.data[k] - xm) / xm #see Glecker et al, eq.2

#-----------------------------------------------------------------------

    def _get_mean_value(self,pos,var):
        """
        calculate mean value for a given observational dataset

        @param pos: position marker
        @type pos: int

        @param var: name of variable to analyze
        @type var: str
        """
        x = []
        for k in self.pos:
            if (self.pos[k] == pos) & ('_' + var + '_' in k):
                x.append(self.data[k])
        x = np.asarray(x)
        return x.mean()



#-----------------------------------------------------------------------

    def plot(self,cmap_name='RdBu_r',vmin=-1.,vmax=1.,nclasses=15,normalize=True,size=10):
        """
        plot Gleckler diagram

        @param cmap_name: name of the colormap to be used
        @type cmap_name: str

        @param vmin: lower boundary for plotting
        @type vmin: float

        @param vmax: upper boundary for plotting
        @type vmax: flaot

        @param nclasses: number of classes for colormap
        @type nclasses: int

        @param size: size of labels
        @type size: int

        @param normalize: normalize data relative to multimodel mean (affects self.data)
        @type normalize: bool
        """

        #~ todo implement normalization here per position

        if normalize:
            self._normalize_data()

        nm = len(self.models); nv = len(self.variables)

        #- colormap
        self.cmap = plt.cm.get_cmap(cmap_name, nclasses)
        self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        gs = gridspec.GridSpec(nm, nv, wspace=0.05,hspace=0.05,bottom=0.2) #generate grid for subplots

        cnt = 0; cnt_m = 0

#        model_list = self.models.sort()

        for model in self.models:
            cnt_m += 1; cnt_v  = 0
            for variable in self.variables:
                ax = self.fig.add_subplot(gs[cnt],frameon=True,aspect='equal',axisbg='grey')
                self.__set_ax_prop(ax)

                #labels
                if cnt_v == 0:
                    ax.set_ylabel(model,size=size,rotation='horizontal')
                if cnt_m == nm:
                    ax.set_xlabel(variable,size=size)

                self.__plot_triangle(ax,self.get_data(variable,model,1),pos='top')    #upper triangle
                self.__plot_triangle(ax,self.get_data(variable,model,2),pos='bottom') #lower triangle
                self.__plot_triangle(ax,self.get_data(variable,model,3),pos='left') #left triangle
                self.__plot_triangle(ax,self.get_data(variable,model,4),pos='right') #right triangle
                cnt += 1; cnt_v += 1

        #--- legend
        #- get positions of subplots to determine optimal position for legend
        def get_subplot_boundaries(g,f):
            x = g.get_grid_positions(f)
            b = x[0]; t = x[1]; l = x[2]; r = x[3]
            return l[0], r[-1], b[-1], t[0]

        left,right,bottom,top = get_subplot_boundaries(gs,self.fig)
        #draw legend
        self._draw_legend(left,right-left)

#-----------------------------------------------------------------------

    def get_data(self,v,m,p):
        """
        return data for a particular model and variable

        @param v: name of variable
        @type v: str

        @param m: model name
        @type m: str

        @param p: position
        @type p: int
        """
        r = None
        k = self.__gen_key(m,v,p)
        if k in self.data.keys():
            r = self.data[k]
        return r

#-----------------------------------------------------------------------

    def get_pos(self,v,m):
        r = None
        k = self.__gen_key(m,v)
        if k in self.pos.keys():
            r = self.pos[k]
        return r

#-----------------------------------------------------------------------

    def __gen_key(self,m,v,p):
        """
        generate a unique key for dictionaries
        comprised of model name, variable name and position

        @param v: name of variable
        @type v: str

        @param m: model name
        @type m: str

        @param p: position
        @type p: int
        """
        if m == None:
            return None
        if v == None:
            return None
        return m.replace(' ','_')+'_'+v.replace(' ','_')+'_'+str(p)

#-----------------------------------------------------------------------

    def add_data(self,v,m,x,pos=1):
        """
        add a data for plotting

        @param v: name of variable
        @type v: str

        @param m: model name
        @type m: str

        @param x: value to be plotted in Gleckler plot
        @type x: float

        @param pos: position where to plot data 1=top triangle, 2=lower triangle
        @type pos: int
        """

        if x != None:
            #- only use valid data
            if v in self.variables:
                if m in self.models:
                    self.data.update({ self.__gen_key(m,v,pos) :x})
                    self.pos.update({ self.__gen_key(m,v,pos) : pos})

#-----------------------------------------------------------------------

    def calc_index(self,x,y,model,variable,weights=None):
        """
        calculate model performance index
        the model and the variable need to have been registered
        already using add_model and add_variable

        @param x: reference data (observation)
        @type x: C{Data} object

        @param y: data to benchmark (e.g. model)
        @type y: C{Data} object

        @param model: model name
        @type model: str

        @param variable: variable name
        @type variable: str

        @param weights: weights to be applied to the data before index calculation; dedicated for spatial area weights
        @type weights: numpy array

        @return: returns performance index aggregated over time
        @rtype float
        """

        if weights == None:
            #set weights according to cell area
            if x.cell_area != None:
                weights = x._get_weighting_matrix()
            else:
                print 'WARNING: no weights when calculating performance index'
                weights = np.ones(x.data.shape)
        else: #weights are given
            if x.cell_area != None:
                print 'WARNING: cell weights are given, while cell_area available from data!!'

        from diagnostic import Diagnostic
        D = Diagnostic(x,y=y)
        e2 = D.calc_reichler_index(weights) #reichler performance index (might return a list if multiple times analyzed)
        if e2 == None:
            return None
        else:
            return np.nansum(e2) #temporal aggregation

#-----------------------------------------------------------------------

    def _draw_legend(self,left,width):
        """
        draw legend for Glecker plot. Requires information on
        the positioning of the colormap axis which can be obtained from

        left,right,bottom,top = get_subplot_boundaries(gs,self.fig)

        @param left: left position of axis
        @type left: float

        @param width: width of the axis to plot colorbar
        @type width: float
        """
        cax = self.fig.add_axes([left,0.05,width,0.05]) #left, bottom, width, height
        cb = mpl.colorbar.ColorbarBase(cax, cmap=self.cmap,
                                   norm=self.norm,
                                   orientation='horizontal')

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def __basemap_ancillary(m,latvalues = None, lonvalues = None,drawparallels=True,drawcountries=True):
    """
    routine to plot ancillary data like coastlines
    or meridians on a basemap plot

    @param m: map to add features to
    @type m: C{Basemap} object

    @param latvalues: latitude values for drawing grid (optional)
    @type latvalues: list or numpy array

    @param lonvalues: longitude values for drawing grid (optional)
    @type lonvalues: list or numpy array

    """

    if latvalues == None:
        latvalues=np.arange(-90.,120.,30.)
    if lonvalues == None:
        lonvalues= np.arange(-180.,180.,90.)
    if drawcountries:
        m.drawcountries()
    m.drawcoastlines()
    m.drawlsmask(lakes=True)
    m.drawmapboundary() # draw a line around the map region
    if drawparallels:
        m.drawparallels(latvalues,labels=[1, 0, 0, 0])
        m.drawmeridians(lonvalues,labels=[0, 0, 0, 1]) # draw meridians

#-----------------------------------------------------------------------

def pm_bar(x,y=None,pcolor='red',ncolor='blue',ax=None,**kwargs):
    """
    generate a nice looking barchart with different color for positive/negative numbers

    @param x: x-variable or variable to plot (if y is not given)
    @param y: y-variable (optional)
    @param ax: axis handle
    @return: returns handle axis
    """

    if ax is None:
        f = pl.figure(); ax = f.add_subplot(111)
    else:
        ax = ax

    if y == None:
        y = x*1.
        x = np.arange(len(y))
    else:
        pass

    yp = y*1.; yp[y<0.] = 0.
    yn = y*1.; yn[y>0.] = 0.

    #--- plot
    ax.bar(x,yp,color=pcolor,edgecolor='None',**kwargs)
    ax.bar(x,yn,color=ncolor,edgecolor='None',**kwargs)

    return ax



#-----------------------------------------------------------------------

def map_season(x,**kwargs):
    """
    generate a seasonal plot
    all arguments are parsed directly to map_plot function

    if kwargs contain a 'figure' argument, then this figure fill be used
    for plotting. Otherwise a new figure will be generated

    @param x: C{Data} object
    @type x : C{Data}

    @return: returns the figure where plot was done
    @rtype: C{figure}

    """

    nvals = len(x.data)
    if nvals == 12:
        year = True
    elif nvals == 4:
        year = False
    else:
        raise ValueError, 'Only data for 4-seasons or monthly data is supported!'

    #/// checks ///
    if x.data.ndim != 3:
        raise ValueError, 'only 3D data supported'

    #/// figure and axes
    if 'figure' in kwargs:
        f = kwargs['figure']
    else:
        f = pl.figure()

    if 'title' in kwargs:
        tit = kwargs.pop('title')
    else:
        tit = x.label

    #/// plot
    if year:
        labels=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    else:
        labels=['JFM','AMJ','JAS','OND']


    #/// in case that an overlay is provided, this needs to be processed for each timestep individually
    if 'overlay' in kwargs.keys():
        overlays = kwargs.pop('overlay')
    else:
        overlays = None

    for i in range(nvals):
        if year:
            ax = f.add_subplot(4,3,i+1)
            if i % 3 == 2:
                show_colorbar = True
            else:
                show_colorbar=False
        else:
            ax = f.add_subplot(2,2,i+1)
            if 'show_colorbar' in kwargs:
                show_colorbar  = kwargs.pop('show_colorbar')
            else:
                show_colorbar = True

        d = x.copy(); d.data = x.data[i,:,:]
        d.label = labels[i]

        if overlays == None:
            overlay = None
        else:
            overlay = overlays[i,:,:]

        map_plot(d,ax=ax,show_colorbar=show_colorbar,overlay = overlay, **kwargs); del d

    f.suptitle(tit,size=16)

    return f




#-----------------------------------------------------------------------

def map_plot(x,use_basemap=False,ax=None,cticks=None,region=None,nclasses=10,cmap_data='jet',
             title=None,regions_to_plot = None,logplot=False,logoffset=None,show_stat=False,
             f_kdtree=False,show_colorbar=True,latvalues=None,lonvalues=None,show_zonal=False,
             zonal_timmean=True,show_timeseries=False,scal_timeseries=1.,vmin_zonal=None,vmax_zonal=None,
             bluemarble = False, contours=False, overlay=None,titlefontsize=14,drawparallels=True,drawcountries=True,show_histogram=False,
             contourf = False, **kwargs):
    """
    produce a nice looking map plot

    @param drawparallels: option to draw parallels on the map
    @type drawparallels: bool

    @param x: data to plot
    @type x: C{Data} object

    @param use_basemap: specifies if Basemap should be used for plotting (=slow), otherwise a simple plot is generated (fast)
    @type use_basemap: bool

    @param ax: axis to plot to; if None, then new figure is generated
    @type ax: matplotlib axis

    @param cticks: ticks for the colorbar
    @type cticks: list of float values

    @param region: region that should be plotted. This is only used in case of Basemap maps
    @type region: C{Region}

    @param nclasses: number of classes for colormap
    @type nclasses: int

    @param cmap_data: colormap for data to be plotted
    @type cmap_data: str

    @param title: title of the plot
    @type title: str

    @param regions_to_plot: This variable might contain a list of regions
                            if the argument is given then each of the regions
                            is plotted as a rectangle into the map
    @type regions_to_plot: list of C{Region} objects

    @param logplot: show data as a logarithmic plot
    @type logplot: bool

    @param logoffset: offset that should be added to the data before performing
                      logarithmic plotting. Useful if negative data
    @type logoffset:  bool

    @param show_stat: show statistic of field in figure title. The mean and std correspond to the
                      SPATIAL mean and stdv of the temporal mean field. It is thus NOT the overall mean
                      and std.!
    @type show_stat: bool

    @param f_kdtree: use kdTree for interpolation of data to grid (might be slow, but might solve problem of stripes in plots)
    @type f_kdtree: bool

    @param latvalues: latitude values for drawing grid (optional)
    @type latvalues: list or numpy array

    @param lonvalues: longitude values for drawing grid (optional)
    @type lonvalues: list or numpy array

    @param vmin_zonal: minimum value for zonal plot
    @type vmin_zonal: float

    @param vmax_zonal: maximum value for zonal plot
    @type vmax_zonal: float

    @param bluemarble: load bluemarble as background image (does only work if no gridded data is shown!)
    @type bluemarble: bool

    @param contours: specifies if plot is done as contour plot
    @type contours: bool

    @param contourf: plot filled contours; is supposed to work only in combination with contours=True
    @type contourf: bool

    @param overlay: overlay for plot (e.g. significance)
    @type overlay: numpy array

    @param show_colorbar: specifies if colorbar should be plotted
    @type show_colorbar: bool

    @param show_zonal: specifies if zonal plots should be shown
    @type show_zonal: bool

    @param titlefontsize: fontsize of the figure title
    @type titlefontsize: int

    @param show_histogram: show a histogram below the map
    @type show_histogram: bool

    @param drawcountries: specifies if countries will be shown in map plot (default=TRUE)
    @type drawcountries: bool


    """

    #--- checks

    if overlay != None:

        if x.data.ndim == 2:
            if overlay.shape != x.data.shape:
                print overlay.shape, x.data.shape
                raise ValueError, 'Invalid geometry for overlay !'
        elif x.data.ndim == 3:
            if overlay.shape != x.data[0,:,:].shape:
                print overlay.shape, x.data.shape
                raise ValueError, 'Invalid geometry for overlay !'
        else:
            raise ValueError, 'Overlay for this geometry not supported!'

    #--- create new figure
    if ax == None:
        fig = plt.figure()

        #with timeseries plot?
        if show_timeseries:
            gs = gridspec.GridSpec(2, 1, wspace=0.05,hspace=0.05,bottom=0.2,height_ratios = [5,1])
            ax  = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
        else:
            ax = fig.add_subplot(111)
    else:
        fig = ax.figure


    #if cmap provided in kwargs, then remove it and set cmap_data
    kwargs1 = kwargs.copy()
    if 'cmap' in kwargs:
        cmap_data = kwargs1.pop('cmap')
    if ('levels' in kwargs) and (contours == False): #levels not needed
        dummy = kwargs1.pop('levels')

    #--- create colormap
    cmap = plt.cm.get_cmap(cmap_data, nclasses)

    #--- temporal mean fields as data to plot
    xm = x.timmean()

    #--- logscale plot ?
    if logplot:
        if logoffset == None:
            if xm.min() < 0.: logoffset = abs(xm.min())*1.01
            else: logoffset = 0.
        else: logoffset = logoffset
        print '     logoffset: ', logoffset
        xm = np.log10( xm + logoffset )

    #--- set projection parameters
    proj='robin'; lon_0=0.; lat_0=0.

    #--- plot using basemap
    if use_basemap:
        llcrnrlon=None; llcrnrlat=None; urcrnrlon=None; urcrnrlat=None

        """ if a region is specfied, then the plotting boundaries are set """
        if region !=None:
            if not hasattr(region,'lonmin'):
                print 'WARNING map boundaries can not be set, as region ' + region.label.upper() + ' has not lat/lon information'
            else:
                dlat = (region.latmax-region.latmin)*0.25; dlon = (region.lonmax-region.lonmin)*0.25
                di = 0. #with 0 it works; for other values problems may occur for negative lon!
                llcrnrlon=region.lonmin - di; llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di; urcrnrlat=region.latmax + di
                proj='tmerc' #use mercator projection at regional scale as robinson does not work!

        ############################################
        # generate Basemap map
        ############################################
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)

        if bluemarble:
            m1.bluemarble()

        if f_kdtree:
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

            #~ omask[dist != 0.] = True

            xm1 = xm.copy(); xm1.shape = (-1)
            Z[idx]   = xm1 #assign data and reshape it and set generate masked array
            #~ Z[dist != 0.] = np.nan

            #~ print Z.shape, omask.shape
            Z[omask] = np.nan
            Z = np.reshape(Z,shape0); Z = np.ma.array(Z,mask=np.isnan(Z))

        else: #f_kdtree --> not kdtree

            #/// in the following, we check if the longitudes are in the right order
            #    to allow for an appropirate plotting. Basemap assumes ascending order
            #    of longitudes. For global projections, problems occur if lon are given
            #    as 0...360 deg. It needs to be shifted then. Otherwise the data will produce
            #    strange stripes when plotting.
            #
            #    REFERENCES:
            #    * http://pl.digipedia.org/usenet/thread/15998/16891/

            if x._lon360: #if lon 0 ... 360, then shift data
                tmp_lon = x._get_unique_lon() #get unique longitudes
                tmplon1 = tmp_lon.copy()
                Z, tmp_lon = shiftgrid(180, xm, tmp_lon, start=False)
                if overlay != None:
                    overlay, nope = shiftgrid(180, overlay, tmplon1, start=False)
                lon, lat = np.meshgrid(tmp_lon, np.arange(Z.shape[0]))
                lat = x.lat
            else:
                print '*** WARNING: not lon360 not validated yet, try KDTREE option if stripes in plot ***'
                lon = x.lon; lat=x.lat
                Z = xm

            X, Y = m1(lon, lat)

        #here is still a problem in the plotting over land; masking does not work properly,
        #while the data as such is o.k.!
        #~ im1=m1.pcolormesh(xmap,ymap,xm,cmap=cmap,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)

        if not bluemarble:
            if contours:
                if 'levels' in kwargs1.keys():
                    levels = kwargs1.pop('levels')
                else:
                    raise ValueError, 'When plotting with contours, you need to specify the levels option (see contour documnetation)'
                if contourf:
                    im1=m1.contourf(X,Y,Z,levels,cmap=cmap,**kwargs1)
                else:
                    im1=m1.contour(X,Y,Z,levels,cmap=cmap,**kwargs1)
                    ax.clabel(im1, inline=1, fontsize=10) #contour label


            else:
                im1=m1.pcolormesh(X,Y,Z,cmap=cmap,**kwargs1) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)

            if overlay != None:
                #print 'Doing overlay plot! ...'
                xcoordnew=X[overlay]; ycoordnew=Y[overlay]
                m1.scatter(xcoordnew,ycoordnew,marker='x',s=50,c='white',edgecolors='white',linewidth=1)
                #todo: there is still a problem that the coordinates are not properly aligned with the grid cells!!!

            #/// ANCILLARY PLOT FOR BASEMAP ///
            __basemap_ancillary(m1,latvalues=latvalues,lonvalues=lonvalues,drawparallels=drawparallels,drawcountries=drawcountries)

    else: #use_basemap = False
        #- normal plots
        if contours:
            if contourf:
                im1 = ax.contourf(xm,cmap=cmap,**kwargs1)
            else:
                im1 = ax.contour(xm,cmap=cmap,**kwargs1)
                ax.clabel(im1, inline=1, fontsize=10) #contour label
        else:
            im1=ax.imshow(xm,cmap=cmap,interpolation='nearest', **kwargs1)

        #--- overlay
        if overlay != None:
            ny,nx = xm.shape
            xnew = arange(nx); ynew = arange(ny)
            XN,YN = np.meshgrid(xnew,ynew)
            del xnew,ynew
            xcoordnew=XN[overlay]; ycoordnew=YN[overlay]
            pl.scatter(xcoordnew,ycoordnew,marker='.',s=1,c='white',edgecolors='white',alpha=1.)

        ax.set_xticks([]); ax.set_yticks([])

    #set legend aligned with plot (nice looking)
    if show_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size="3%", pad=0.1, axes_class=maxes.Axes)
        ax.figure.add_axes(cax)
        vmin = im1.get_clim()[0]; vmax = im1.get_clim()[1]
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cb   = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,ticks=cticks)


    #--- add a histogram below the map plot
    if show_histogram:
        add_histogram(ax,x)

    #Zonal plot
    if show_zonal:
        add_zonal_plot(ax,x,timmean=zonal_timmean,vmin=vmin_zonal,vmax=vmax_zonal) #,vmin=im1.get_clim()[0],vmax=im1.get_clim()[1])






    def _add_region(m,r,color='red'):
        """
        plot region r on top of basemap map m

        @param m: map
        @type m: C{Basemap} object

        @param r: region to plot
        @type r: C{Region}

        @param color: color to plot region
        @type color: str
        """
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
        title = x._get_label()
    else:
        pass

    #--- show field statistics in title ?
    # calculates first the temporal mean and results shown then
    # are the means and std of the temporal mean fields which are weighted
    # appropriately according to the cell area
    if show_stat:
        tmp_xm = x.copy()
        tmp_xm.data = xm.reshape(1,xm.shape[0],xm.shape[1]) #make a 3D object thus fldmean works appropriately
        me = tmp_xm.fldmean(); st=tmp_xm.fldstd()
        assert(len(me) == 1)
        assert(len(st) == 1)
        me = me[0]; st=st[0]
        title = title + '\n ($' + str(round(me,2))  + ' \pm ' + str(round(st,2)) + '$' + ')'


    ax.set_title(title,size=titlefontsize)


    #/// show timeseries? ///
    if show_timeseries:
        ax2.plot(plt.num2date(x.time),x.fldmean())
        ax2.grid()
        ax2.set_ylim(im1.get_clim()[0]*scal_timeseries,im1.get_clim()[1]*scal_timeseries)
        ti = ax2.get_yticks(); n=len(ti) / 2
        ax2.set_yticks([ti[0],ti[n],ti[-1]])
        ax2.set_ylabel(x._get_unit())

    return fig

#-----------------------------------------------------------------------------------------------------------------------

def add_histogram(ax,x):
    """
    add a histogram to an axis

    @param ax: axis to add the histogram plot to
    @type ax: axis

    @param x: Data that should be vizualized as histogram
    @type x: Data
    """

    divider = make_axes_locatable(ax)
    #zax     = divider.new_vertical("30%", pad=0.1, axes_class=maxes.Axes,pack_start=True)
    zax = divider.append_axes("bottom","30%",pad=0.1)


    ax.figure.add_axes(zax,axisbg=ax.figure.get_facecolor())

    H = HistogrammPlot(ax=zax) #bins ???
    H.plot(x) #plot options ????


#-----------------------------------------------------------------------------------------------------------------------

def add_zonal_plot(ax,x,timmean=True,vmin=None,vmax=None):
    """
    add a zonal plot to the axis

    @param ax: axis where zonal plot should be added to
    @type ax: axis

    @param x: data to plot
    @type x: C{Data} object

    @param timmean: temporal mean for zonal plot?
    @type timmean: bool

    @param vmin: minimum value for y-axis
    @type vmin: float

    @param vmax: maximum value for y-axis
    @type vmax: float


    """

    divider = make_axes_locatable(ax)
    zax     = divider.new_horizontal("15%", pad=0.1, axes_class=maxes.Axes,pack_start=True)
    ax.figure.add_axes(zax,axisbg=ax.figure.get_facecolor())

    ZP = ZonalPlot(ax=zax,dir='y')

    if x.data.ndim == 2:
        weights = np.ones(x.data.shape)
    elif x.data.ndim == 3:
        nt,ny,nx = x.data.shape
        weights = np.ones((ny,nx))
    weights = np.ma.array(weights,mask = weights != weights)
    ZP.plot(x,timmean=timmean,show_ylabel=False)


    #- set limits
    if ((vmin == None) & (vmax == None)):
        vmin = zax.get_xlim()[0]
        vmax = zax.get_xlim()[1]
        #symmetry if neg. and posistive limits
        if (vmin < 0.) & (vmax>0.):
            val = max(abs(vmin),abs(vmax))
            vmin = -val; vmax = val

    if vmin == None:
        vmin = zax.get_xlim()[0]
    if vmax == None:
        vmax = zax.get_xlim()[1]
    zax.set_xlim(vmin,vmax)

    #set only first and last label
    zax.set_xticks([vmin,vmax])
    zax.plot([0,0],zax.get_ylim(),linestyle='-',color='grey')

    for tick in zax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)

    return zax

#-----------------------------------------------------------------------

def add_nice_legend(ax,im,cmap,cticks=None,dummy=False,fontsize=8):
    """
    add a nice looking legend

    @param ax: major axis with plot
    @type ax: matpltlib axis object

    @param im: result from command like imshow
    @param im: matplotlib im object (???)

    @param dummy: add colorbar axis as a dummy axis which is not visible
                  this is useful if you have multiple subplots which should
                  have the same size. Adding a colorbar will slightly change the size
    @type dummy: bool

    @param fontsize: fontsize for colorbar ticks
    @type fontsize: int

    @param cmap: colormap
    @type cmap: str or colorbar class

    @param cticks: colorbar ticks; if None, then default setup is used
    @type cticks: list


    #todo: add option to add units
    """

    #set legend aligned with plot (nice looking)
    divider = make_axes_locatable(ax)
    #cax = divider.new_horizontal("5%", pad=0.05, axes_class=maxes.Axes)

    cax = divider.append_axes("right","5%",pad=0.05)


    ax.figure.add_axes(cax,axisbg=ax.figure.get_facecolor())
    if dummy:
        cax.set_xticks([])
        cax.set_yticks([])
        cax.set_frame_on(False)
    else:
        norm = mpl.colors.Normalize(vmin=im.get_clim()[0], vmax=im.get_clim()[1])
        cb   = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,ticks=cticks)

        for t in cb.ax.get_yticklabels():
            t.set_fontsize(fontsize)


#-----------------------------------------------------------------------

def hov_difference(x,y,climits=None,dlimits=None,data_cmap='jet',nclasses=15,cticks=None,cticks_dif=None,ax1=None,ax2=None,ax3=None,rescaley=6,grid=True,rescalex=1,**kwargs):
    """
    class to plot hovmoeller diagrams of two datasets
    and their difference

    x,y two Data structures

    axextra: plot difference on separate axis
    """

    if climits == None:
        sys.exit('Please specify climits for hovmoeller')
    if dlimits == None:
        sys.exit('Please specify dlimits for hovmoeller')

    fig = plt.figure()
    if ax1 == None:
        ax1 = fig.add_subplot(311)
    if ax2 == None:
        ax2 = fig.add_subplot(312)
    if ax3 == None:
        ax3 = fig.add_subplot(313)

    #set all invalid data to NAN
    xdata = x.data; ydata = y.data

    hov1 = hovmoeller(num2date(x.time),xdata,rescaley=rescaley,lat=x.lat,rescalex=rescalex)
    hov2 = hovmoeller(num2date(y.time),ydata,rescaley=rescaley,lat=y.lat,rescalex=rescalex)

    hov1.time_to_lat(**kwargs); hov2.time_to_lat(**kwargs)

    cmap = plt.cm.get_cmap(data_cmap, nclasses)

    hov1.plot(title=x._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax1,showcolorbar=False,climits=climits,grid=grid)
    hov2.plot(title=y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax2,showcolorbar=False,climits=climits,grid=grid)

    add_nice_legend(ax1,hov1.im,cmap,cticks=cticks); add_nice_legend(ax2,hov2.im,cmap,cticks=cticks)

    if x.data.shape == y.data.shape:
        hov3 = hovmoeller(num2date(y.time),x.data - y.data,rescaley=rescaley,lat=y.lat,rescalex=rescalex)
        hov3.time_to_lat(**kwargs)
        cmap_diff = plt.cm.get_cmap('RdBu', nclasses)
        hov3.plot(title=x._get_label() + ' - ' + y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap_diff,ax=ax3,showcolorbar=False,climits=dlimits,grid=grid)
        add_nice_legend(ax3,hov3.im,cmap_diff,cticks=cticks_dif)
    else:
        msg = 'Difference plot not possible as data has different shape'
        ax3.text(0.5, 0.5,msg,
             horizontalalignment='center',
             verticalalignment='center') #,
             #transform = ax.transAxes)
        ax3.set_xticks([]); ax3.set_yticks([])

    return fig,hov1,hov2


#-----------------------------------------------------------------------


def map_difference(x,y,dmin=None,dmax=None,use_basemap=False,ax=None,title=None,cticks=None,region=None,nclasses=10,cmap_data='jet',cmap_difference = 'RdBu_r',rmin=-1.,rmax=1., absthres=None, show_stat=True,show_zonal=True,zonal_timmean=False, **kwargs):
    """
    Given two datasets, this map generates a map plot of each dataset as
    well as of the difference of the two datasets

    @param x: first dataset
    @type x: C{Data} object

    @param y: second dataset
    @type y: C{Data} object

    @param dmin: minimum value of difference map
    @type dmin: float

    @param dmax: maximum value of difference map
    @type dmax: float

    @param use_basemap: flag if Basemap should be used for plotting
    @type use_basemap: bool

    @param ax: axis to plot to; if None, then new figure is generated
    @type ax: matplotlib axis

    @param title: title of the plot
    @type title: str

    @param cticks: ticks for the colorbar
    @type cticks: list of float values

    @param region: region that should be plotted. This is only used in case of Basemap maps
    @type region: C{Region}

    @param nclasses: number of classes for colormap
    @type nclasses: int

    @param cmap_data: colormap for data to be plotted
    @type cmap_data: str

    @param cmap_difference: colormap for difference map to be plotted
    @type cmap_difference: str

    @param rmin: minimum value for data plot
    @type rmin: float

    @param rmax: maximum value for data plot
    @type rmax: float

    @param absthres: threshold that will be estimated based on
                     absolute difference and will be applied to
                     relative difference maps
    @type absthres: float

    @param show_zonal: plot zonal statistic plot
    @type show_zonal: bool

    """


    if 'cticks_diff' in kwargs:
        cticks_diff = kwargs.pop('cticks_diff')
    else:
        cticks_diff = None

    fig = plt.figure()

    ax1 = fig.add_subplot(221); ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223); ax4 = fig.add_subplot(224)

    #--- get colormap
    cmap = plt.cm.get_cmap(cmap_data, nclasses)

    #- temporal mean fields
    xm = x.timmean(); ym = y.timmean()

    proj='robin'; lon_0=0.; lat_0=0.

    #- plot first dataset
    map_plot(x,use_basemap=use_basemap,ax=ax1,cticks=cticks,region=region,nclasses=nclasses,cmap_data=cmap_data, title=title,show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean, **kwargs)

    #- plot second dataset
    map_plot(y,use_basemap=use_basemap,ax=ax2,cticks=cticks,region=region,nclasses=nclasses,cmap_data=cmap_data, title=title,show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean,  **kwargs)

    #-first minus second dataset
    adif = x.sub(y) #absolute difference #todo where to get std of seasonal means !!!! needs to be realized before beeing able to use significance ????

    map_plot(adif,use_basemap=use_basemap,ax=ax3,vmin=dmin,vmax=dmax,cticks=cticks_diff,region=region,nclasses=nclasses,cmap_data=cmap_difference, title='absolute difference [' + x.unit + ']',show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean)

    #- relative error
    rdat = adif.div(x) #y.div(x).subc(1.) #relative data
    if absthres != None:
        mask = abs(x.timmean()) < absthres
        rdat._apply_mask(~mask)

    map_plot(rdat,use_basemap=use_basemap,ax=ax4,vmin=rmin,vmax=rmax,title='relative difference',cticks=[-1.,-0.75,-0.5,-0.25,0.,0.25,0.5,0.75,1.],region=region ,nclasses=8,cmap_data=cmap_difference,show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean)


    return fig

#-----------------------------------------------------------------------

def plot_hovmoeller(x,rescaley=10,rescalex=1,monthsamp=24,dlat=1.,cmap=None,ax=None,climits=None,xtickrotation=0,showxticks=True,title=None,cticks=None,ylabel=None):
    """
    plot hovmoeller plots given a C{Data} object

    @param x: C{Data} object
    @type x: C{Data} object

    @param rescalex: rescaling parameter for x-axis for hovmoeller plot
    @type rescalex: int

    @param rescaley: rescaling parameter for y-axis for hovmoeller plot
    @type rescaley: int


    """

    if climits == None:
        raise ValueError, 'climits need to be given!'

    if cmap == None:
        cmap = plt.cm.get_cmap('RdBu', 15)
    if ax == None:
        f = plt.figure()
        ax = f.add_subplot(111)

    if title == None:
        tit=x._get_label()
    else:
        tit = title

    if showxticks:
        xlab = 'time'
    else:
        xlab=None

    if ylabel == None:
        ylabel = 'lat'



    h = hovmoeller(pl.num2date(x.time),x.data,rescaley=rescaley,lat=x.lat,rescalex=rescalex)
    h.time_to_lat(dlat=dlat,monthly = True, yearonly = True,monthsamp=monthsamp)
    h.plot(title=tit,xlabel=xlab,ylabel=ylabel,origin='lower',xtickrotation=xtickrotation,cmap=cmap,ax=ax,showcolorbar=False,climits=climits,grid=False,showxticks=showxticks)

    add_nice_legend(ax,h.im,cmap,cticks=cticks)

    return h



