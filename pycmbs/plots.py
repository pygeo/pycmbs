# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""


"""
Module that contains relevant classes for diagnostic plots

@todo: implement writing of statistics to an ASCII file as export
@todo: implement taylor plots
@todo: faster implementation of Basemap plots. For large number of grid cells, the current KTree implementation is by far too slow!

"""
import matplotlib as mpl
mpl.use('agg')

from pycmbs.data import *
from pycmbs.hov import *

from matplotlib import pyplot as plt
from matplotlib import pylab

from matplotlib.patches import Polygon
import matplotlib.path as mpath
from matplotlib.collections import PatchCollection

from mpl_toolkits.basemap import Basemap,shiftgrid
from scipy import stats

import numpy as np
import matplotlib.gridspec as grd

from matplotlib.patches import Circle
import matplotlib.patches as mpatches
import sys
import copy
import numpy as np


from scipy.spatial import cKDTree as KDTree #import the C version of KDTree (faster)
from matplotlib.ticker import MaxNLocator

import matplotlib.gridspec as gridspec

#http://old.nabble.com/manual-placement-of-a-colorbar-td28112662.html
from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes
import datetime



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


def rotate_ticks(ax,angle):
    """
    Rotate ticks of an axis
    @param ax:
    @param angle:
    @return:
    """
    ts = ax.xaxis.get_major_ticks()
    for t in ts:
        t.label.set_rotation(angle)


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class CorrelationAnalysis(object):
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

            if ax is None:
                f = plt.figure()
                self.ax = f.add_subplot(111)
            else:
                self.ax = ax

#-----------------------------------------------------------------------

        def do_analysis(self):
            """
            perform correlation analysisglobal

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

            #print 'RMSE: ', rmse
            #print 'R   : ', r
            #print 'N   : ', n

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class HovmoellerPlot(object):
    def __init__(self, D, rescaley=10, rescalex=10, yticksampling=1,
                 monthly=False, ax=None, figsize=(10, 5)):
        """
        D : Data object
            Data object that should be used for plotting the Hovmoeller
            diagram


        if the argument lat is provided it is assumed that lat/lon are 2D matrices
        In this case the value is expected to be a 3D variables as
        value(time,ny,nx)
        """
        if ax is None:
            self.figure = pl.figure(figsize=figsize)
            self.ax = self.figure.add_subplot(111)
        else:
            self.ax = ax
            self.figure = self.ax.figure

        self.hov = hovmoeller(D.num2date(D.time), None,
                              rescaley=rescaley, rescalex=rescalex)
        #self.hov.time_to_lat(dlat=dlat,yticksampling=yticksampling,monthly=monthly)
        self.x = D

    def plot(self, title=None, climits=None, showxticks=True,
             showcolorbar=True, cmap='jet', xtickrotation=90,
             ylim=None):
        if climits is None:
            raise ValueError('CLIMITS needs to be specified!')
        self.hov.plot(input=self.x, ax=self.ax, title=title,
                      ylabel='lat', xlabel='days', origin='lower',
                      xtickrotation=xtickrotation, climits=climits,
                      showxticks=showxticks, showcolorbar=showcolorbar,
                      cmap=cmap)
        if ylim is not None:
            self.ax.set_ylim(ylim)
        if title is None:
            self.ax.set_title(self.x.label)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class ReichlerPlot(object):
    """
    class for Reichler plot generation

    @todo: add example how to use Reichler plotting
    @todo: provide references Glecher and Reichler + Kim
    """
    def __init__(self, ax=None):
        """
        constructor for Reichler plot

        @param ax: axis to plot data to; if None, new figure will be generated
        @type ax: matplotlib axis
        """

        if ax is None:
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

        if self.e2[0] is None:
            return self.ax.figure



        #- normalize results (relative modle performance)
        self._normalize()
        x = np.arange(len(self.e_norm))
        y1 = self.e_norm*100.; y2 = self.e_norm*100.

        y1 = np.ma.array(y1,mask=y1<0.) #positive values only
        y2 = np.ma.array(y2,mask=y2>=0.) #negative values only

        #print 'Reichler data for plotting: ', y1,y2
        #print 'Original Reichler data:'
        #print self.e2
        #print self.e2_norm

        if 'color' in kwargs.keys():
            upcolor  =kwargs.pop('color')
            lowcolor =kwargs.pop('color')
        else:
            upcolor='red'
            lowcolor='blue'

        self.ax.bar(x,y1,color=upcolor ,edgecolor='None',**kwargs)
        self.ax.bar(x,y2,color=lowcolor,edgecolor='None',**kwargs)

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
        """
        do a very simple plot of diagnostics
        """
        for i in np.arange(len(self.e2)):
            self.ax.plot(self.e2[i],'o',label=self.labels[i])

#-----------------------------------------------------------------------

    def circle_plot(self):
        """
        nice looking plot of Reichler index
        """
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
        E = []

        for e2 in self.e2:
            if len(e2) != n:
                print 'WARNING: non consistent length in error statistics!!!'
            E.append(np.nansum(np.sqrt(e2))) #temporal aggregation

        E = np.asarray(E);  EM = E.mean() #take square root, as e2 is still the squared error!
        self.e_norm =  (E - EM) / EM #see Glecker et al, eq.2

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class ScatterPlot(object):
    """
    Class for generation of scatterplots
    """
    def __init__(self, x, ax=None, ticksize=10, normalize_data=False, show_xlabel=True, figsize=None):
        """
        constructor of class C{ScatterPlot}

        @param x: Variable that will be used as the x-variable
        @type x: C{Data} object

        @param normalize_data: if True, then the dataseries is normalized internally so that it has zero mean and a std of 1
        @type normalize_data: bool

        @param show_xlabel: show xlabel in plot
        @type show_xlabel: bool

        @param ticksize: ticksize of labels
        @type ticksize: int
        """

        self.show_xlabel = show_xlabel

        if ax is None:
            f = plt.figure(figsize=figsize)
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

        self.figure = self.ax.figure
        self.x = x
        self.lines = []
        self.labels = []
        self.ticksize = ticksize
        self.normalize = normalize_data

#-----------------------------------------------------------------------

    def __normalize_data(self, x):
        """
        normmalize timeseries
        """
        return (x-x.mean()) / x.std()

#-----------------------------------------------------------------------

    def plot(self, y, regress=True, fldmean=True, hexbin=False, **kwargs):
        """
        add a dataset to the scatterplot and plot
        it. It also allows to perform a regression analysis

        @param y: data to be plotted on y-axis
        @type y: C{Data} object

        @param regress: Perform linear regression analysis
        @type regress: bool

        @param fldmean: show in scatterplot fldmean() values, else datapairs are copnstructed for each grid cell
        @type fldmean: bool
        """

        if y.label is None:
            label=''
        else:
            label=y.label

        if fldmean:
            xdat = self.x.fldmean()
            ydat = y.fldmean()
        else:
            if self.x.data.shape != y.data.shape:
                print self.x.data.shape
                print self.y.data.shape
                raise ValueError, 'Invalid geometry between X and Y. fldmean=True option therefore not possible!'
            else:
                xdat = self.x.data.flatten()
                ydat = y.data.flatten()

        assert(isinstance(xdat,np.ma.core.MaskedArray))
        assert(isinstance(ydat,np.ma.core.MaskedArray))

        #--- mask invalid data
        msk = np.isnan(xdat) | np.isnan(ydat)
        xdat = np.ma.masked_where(msk, xdat)
        ydat = np.ma.masked_where(msk, ydat)

        if self.normalize:
            xdat = self.__normalize_data(xdat)
            ydat = self.__normalize_data(ydat)

        #- calculate linear regression
        if regress:
            slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(xdat, ydat)
            #~ print 'r_value: ', r_value
            #~ print xdat
            #~ print ydat
            nval = (~(xdat-ydat).mask).sum()  # number of valid datasets used for comparison

            assert(isinstance(xdat,np.ma.core.MaskedArray))
            assert(isinstance(ydat,np.ma.core.MaskedArray))

            rms_error = calc_rms_error(xdat,ydat)
            bias,c_rms = calc_centered_rms_error(xdat,ydat)
            std_error = np.std(xdat-ydat)

            if p_value < 0.01:
                spvalue = 'p < 0.01'
            else:
                spvalue = 'p=' + str(round(p_value,2))

            if r_value is None:
                label = ''
            else:
                label = '\n' + label + '\nr=' + str(round(r_value,2)) + ', ' + spvalue + ', ' + 'rmsd: ' + str(rms_error) + ', N=' + str(int(nval)) + '\n' + 'y=' + str(slope) + 'x+' + str(intercept) + ''

        # actual plot
        if hexbin:
            l = self.ax.hexbin(xdat, ydat,**kwargs)
        else:
            l = self.ax.plot(xdat, ydat,'.', label=label, **kwargs)[0]

        if hexbin:
            pass
        else:
            self.lines.append(l); self.labels.append(label)

        if regress:
            if r_value is not None:
                if hexbin:
                    self.ax.plot(xdat,xdat*slope+intercept,'--')
                else:
                    self.ax.plot(xdat,xdat*slope+intercept,'--', color=l.get_color())

        if self.show_xlabel:
            self.ax.set_xlabel(self.x._get_label(),size=self.ticksize )
        self.ax.set_ylabel(y._get_unit(),size=self.ticksize)

        self._change_ticklabels()

        if regress:
            return r_value,p_value,rms_error,c_rms, std_error, nval, np.ma.std(xdat.flatten()), np.ma.std(ydat.flatten())
        else:
            return None

#-----------------------------------------------------------------------

    def _change_ticklabels(self):
        for tick in self.ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(self.ticksize)
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)

#-----------------------------------------------------------------------

    def legend(self,size=8.):
        """
        plot legend
        """
        self.ax.legend(self.lines,self.labels,prop={'size':size})

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class LinePlot(object):
    """
    class for a pyCMBS Line Plot
    This class is usefull for plotting timeseries
    """

    def __init__(self, ax=None, regress=False, title=None,
                 show_xlabel=True, show_ylabel=True, ticksize=10,
                 normx=1., show_equation=True, xtickrotation=90,
                 figsize=(8, 7)):
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

        if ax is None:
            f = plt.figure(figsize=figsize)
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
        self.ticksize = ticksize
        self.normx=normx
        self.show_equation = show_equation
        self.xtickrotation = xtickrotation

#-----------------------------------------------------------------------

    def legend(self, prop={'size':8}, **kwargs):
        """
        plot legend
        """
        self.ax.legend(self.lines, self.labels, prop=prop, **kwargs)

#-----------------------------------------------------------------------

    def _change_ticklabels(self, ax=None):
        if ax is None:
            ax = self.ax

        # yaxis
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)
        # xaxis
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)
            tick.label.set_rotation(self.xtickrotation)

#-----------------------------------------------------------------------

    def plot(self, x, ax=None, vmin=None, vmax=None, label = None, norm_std = False, set_ytickcolor=True, **kwargs):
        """
        plot LinePlot data. If a spatial field is provided, this is aggregated
        using the fldmean() function of C{Data}

        @param x: data to be plotted
        @type x: C{Data}

        @param ax: axis to plot to. If None, then a new figure is generated.
                   be aware that problems might occur when using axes that were generated with
                   the twinx() command. Problems can be avoided by generating the figure axcopy=ax.twinx() just
                   immediately BEFORE calling the plot routine.
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

            if ax is None:
                ax = self.ax
                set_axiscolor=False
            else:
                ax = ax
                set_axiscolor=True

            if x.ndim == 1: #if a vector already provided
                y = x.data * 1.
            else:
                y = x.fldmean() #... otherwise use fldmean() to get timeseries

            if norm_std:
                y /= y.std()

            if label is None:
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

            p = ax.plot(x.date, y , label=label, **kwargs)[0]
            self.lines.append(p)
            if self.regress:
                ax.plot(x.date,x.time*slope+intercept,'--',color=p.get_color()) #plot regression line

            if self.show_ylabel:
                ax.set_ylabel(x._get_unit(),size=self.ticksize)
            if self.show_xlabel:
                ax.set_xlabel('time',size=self.ticksize)

            if self.title != None:
                ax.set_title(self.title,size=self.ticksize)

            if vmin != None and vmax != None:
                ax.set_ylim(vmin,vmax)

            if set_ytickcolor:
                for tl in ax.get_yticklabels():
                    tl.set_color(p.get_color())

            #self._change_ticklabels(ax)


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class GlobalMeanPlot(object):
    """
    plots timeseries of global mean field
    """

    def __init__(self,ax=None,climatology=True,ax1=None):
        """
        @param ax: specifies axis to plot the data to
        @type ax: axis

        @param ax1: specifies axis for second plopt (only used when climatology==True)
        @type ax1: axis

        @param climatology: specifies if a second plot for a climatological mean value shall be generated
        @type climatology: bool
        """

        if climatology:
            nplots = 2
        else:
            nplots = 1
        self.climatology = climatology
        self.nplots=nplots

        if ax is None:
            #--- create new figure if needed
            f = plt.figure()
            self.ax = f.add_subplot(nplots,1,1)
            if self.climatology:
                if ax1 is None:
                    self.ax1 = self.ax.figure.add_subplot(nplots,1,2)
                else:
                    self.ax1 = ax1
        else: #figure already existing
            self.ax = ax
            if self.climatology:
                if ax1 is None:
                    raise ValueError, 'If option climatology is chosen for GlobalMeanPlot and axis is provided, then axis1 needs also be provided!'
                else:
                    self.ax1 = ax1

        self.labels=[]; self.plots=[]
        self.pdata={}
        self.pdata_clim={}

#-----------------------------------------------------------------------------------------------------------------------

    def plot_mean_result(self,dt=0.,colors=None,plot_clim=False):
        """
        plot mean result

        dt: difference for coregistering in fractions of days

        @param plot_clim: if True, then the climatology will be plotted instead of the actual data
        @type plot_clim: bool

        @return:returns figure handle
        """

        def coregister(x, x1, y1, d):
            """
            This routine provides the functionality to coregister two timeseries of data
            x : reference x vector
            x1: x-value of data
            y1: y-value of data
            d : difference in t (threshold)
            """

            o=[]
            for xx in x:
                m1=abs(x1-xx)<d
                m2=~np.isnan(y1)
                m = m1 & m2
                if sum(m)>0:
                    o.append(np.mean(y1[m]))
                else:
                    o.append(np.nan)
            o=np.asarray(o)
            return o

        if hasattr(self,'pdata'):
            pass
        else:
            raise ValueError('Can not plot mean results for GlobalMeanPlot! Missing data!')

        f = plt.figure()
        ax = f.add_subplot(111)
        if plot_clim:
            pdata = self.pdata_clim
            print('GlobalMeanPlot climdata: %s ' % pdata)
        else:
            pdata = self.pdata

        groups = pdata.keys()
        for g in groups:
            dat = pdata[g] #this gives a list, where each entry is a dictionary of ['time','data','unit']
            n = 0
            for i in xrange(len(dat)):
                if i == 0:
                    if plot_clim:
                        tref = dat[i]['time']  # climatology 1...12
                    else:
                        tref = pylab.date2num(dat[i]['time'])  # reference time vector
                    y = dat[i]['data']*1.
                    ys = y*y
                    n = np.ones(len(y))*1.
                else:
                    # interpolate results to reference time period
                    if plot_clim:
                        t1 = dat[i]['time']
                    else:
                        t1 = pylab.date2num(dat[i]['time'])
                    y1 = dat[i]['data']*1.
                    yo = coregister(tref, t1, y1, dt)

                    m = ~np.isnan(yo)

                    y[m]=y[m]+yo[m]
                    ys[m] = ys[m] + yo[m]*yo[m]
                    n[m] += 1.

                    if plot_clim:
                        pass
                    del m


            if len(n) > 0:
                n = map(float,n)
                ym = y / n
                ys /=  n #squared values
                std_data = np.sqrt(ys - ym*ym)

                color=None
                if colors is not None:
                    if g in colors.keys():
                        color=colors[g]

                if plot_clim:
                    tval = tref
                else:
                    tval = pylab.num2date(tref)

                ax.fill_between(tval,ym-std_data,ym+std_data,color=color,alpha=0.5)
                ax.plot(tval,ym,label=g+'$\pm 1\sigma$',color=color)
                ax.set_ylabel(dat[0]['unit'])
                ax.set_xlabel('months')

        ax.legend()
        ax.grid()

        return f

#-----------------------------------------------------------------------------------------------------------------------

    def plot(self, D1, color=None, linewidth=1, show_std=False,
                label=None, linestyle='-', mask=None, group='A',
                stat_type='mean'):
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

        @param group: specifies a group that will be used to combine plots of the same type. The group is used as a key for a dictionary that stores the results
        @type group: str

        @param stat_type: specifies which statistic shall be plotted ['mean','sum'], either area weighted mean or sum
        @type stat_type: str

        """

        if stat_type not in ['mean','sum']:
            raise ValueError, 'Invalid stat_type in GlobalMean plot'

        if (label is None) and (D1.label in self.labels):
            #print 'Label already existing: ', D.label, ' skipping analysis'
            return
        elif ((label is not None) and (label in self.labels)):
            #print 'Label already existing: ', label, ' skipping analysis'
            return

        # ensure to work with a data object
        if 'tuple' in str(type(D1)): #only a vector is provided as data together with time (time,data)
            D = D1[2]  # (time,data,orgdata)
        else:
            D = D1

        if D.data.ndim != 3:
            raise ValueError('Global plot only supported for 3D data')

        if mask is not None:
            D._apply_mask(mask)

        if stat_type == 'mean':
            # mean field
            m = D.fldmean(return_data=True)  # mean
        elif stat_type == 'sum':
            m = D.areasum(return_data=True)  # area weighted sum
        else:
            raise ValueError, 'Unsupported stat_type: ' + stat_type
        mdata = m.data.flatten()

        #time
        #thlp = m.num2date(m.time) #this gives a datetime object based on the calendar of the dataset first and this is then converted to the standard matplotlib datetime object
        #t = [datetime.datetime(x.year,x.month,x.day,x.hour,x.minute,x.second) for x in thlp] #datetime object
        t = m.date
        #del thlp

        #--- plot generation ---
        if color is None:
            p = self.ax.plot(t,mdata,linewidth=linewidth,linestyle=linestyle)
        else:
            p = self.ax.plot(t,mdata,color=color,linewidth=linewidth,linestyle=linestyle)

        if group in self.pdata.keys():
            vdata = self.pdata[group]
        else:
            vdata=[]
        vdata.append({'time':t, 'data':mdata,'unit':m._get_unit()})

        #print 'vdata: ', t
        self.pdata.update({group:vdata}) #store results for current group
        del vdata

        if show_std & (stat_type == 'mean'):
            s = D.fldstd (return_data=True) #std
            sdata = s.data.flatten()
            self.ax.fill_between(t,mdata-sdata,y2=mdata+sdata,color=p[0].get_color(),alpha=0.5)

        # plot climatology if desired
        if self.climatology:
            if hasattr(D,'time_cycle'):
                tmp = D.get_climatology(return_object=True)
                tmp.adjust_time(year=1700, day=15)
                tmp.timsort()
                m = tmp.fldmean(return_data=False).flatten()

                self.ax1.plot(np.arange(1,13), m, linestyle=linestyle)
                self.ax1.set_xlim(0.,13.)
                self.ax1.set_ylabel(tmp._get_unit())
                self.ax1.set_xlabel('months')

                self.ax1.grid()

                # store values for aggregated plot
                if group in self.pdata_clim.keys():
                    vdata = self.pdata_clim[group]
                else:
                    vdata = []
                vdata.append({'time':np.arange(1,13),'data':m,'unit':tmp._get_unit()})
                self.pdata_clim.update({group:vdata}) #store results for current group
                del vdata, tmp
            else:
                print 'WARNING: Global mean plot can not be generated due to missing time_cycle!'


        # store information for legend
        self.plots.append(p[0])
        if label is None:
            self.labels.append(D.label)
        else:
            self.labels.append(label)

        # labels
        self.ax.set_ylabel(D._get_unit())

       # LEGEND always below the figure
        if self.nplots == 1:
            lax = self.ax
            loff = 0.2
        else:
            lax = self.ax1
            loff = 0.4
        box = lax.get_position()

        lax.figure.subplots_adjust(bottom=loff)  # make space on bottom for legend
        lax.legend(self.plots,self.labels,loc='upper center', bbox_to_anchor=(0.5, -loff),fancybox=True, shadow=True, ncol=3,prop={'size':8})

        self.ax.grid()

#-----------------------------------------------------------------------

class HistogrammPlot(object):
    """
    class to plot histograms based on C{Data} objects
    """
    def __init__(self, ax=None, bins=10, normalize=False, percent=True):
        """
        @param ax: axis to plot data to. If not specified, then a new figure is created
        @type: ax: axis

        @param bins: bins for histogram calculation, either int or a list
        @type bins: int or list or array

        @param normalize: specifies if data should be normalized relative to the sample size
        @type normalize: bool

        @param percent: resulting frequencies in percent (applies only if normalize=True)
        @type percent: bool
        """

        #- Figure init
        if ax is None:
            self.figure = pl.figure()
            self.ax = self.figure.add_subplot(111)
        else:
            self.ax = ax
            self.figure = self.ax.figure

        self.bins = bins
        self.normalize = normalize
        self.percent = percent

    def plot(self, X, color='black', linestyle='-', linewidth=1., label='', shown=False, show_legend=False, **kwargs):
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

        #--- calculate frequency distribution
        f,b = np.histogram(x,bins = self.bins,**kwargs)
        if self.normalize:
            f /= float(sum(f))
            if self.percent:
                f *= 100.

        self.ax.plot(b[0:-1], f, color=color, linestyle=linestyle,
                    linewidth=linewidth, label=label)

        if show_legend:
            self.ax.legend()


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


class ZonalPlot(object):
    def __init__(self, ax=None, dir='y'):
        """
        @todo: still needs to take into account appropriately area weighting

        constructor for zonal plot

        dir - specifies direction for aggregation: y = zonal, x = meridional aggregation

        CAUTION: the class simply aggregates x/y. Thus the user needs to ensure, that the data is projected
        in a way that all lon/lat are on the same row/col
        """

        # directionalities
        if dir == 'y': # zonal plot
            self.dir = 'y'
        elif dir == 'x':
            self.dir = 'x'
        else:
            raise ValueError('Invalid value for agregation direction (ZonalPlot): ' % dir)

        # set axis
        if ax is None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

#-----------------------------------------------------------------------

    def plot(self, x, xlim=None, timmean=False, show_ylabel=True, label=''):
        """
        plot zonal plot

        @param x: data to be plotted
        @type x: C{Data} object

        @param xlim: limits for the x-axis (e.g. values)
        @type xlim: tuple

        @param timmean: temporal mean calculation
        @type timmean: bool

        """

        if not x._latitudecheckok:
            print('WARNING: can not do zonal plot as not regular latitudes!')
            return

        # check if all latitudes are the same
        lu = x.lat.mean(axis=1)
        if any( abs(lu - x.lat[:, 0]) > 1.E-5):
            print 'WARNING: latitudes are not unique!!!'
            print lu.shape,x.lat.shape
            print lu

            print x.lat[:, 0]
            print x.lat[:, 0] - lu
            raise ValueError('Cannot work with NOT unique latitude values!')

        if timmean:
            thex = x.timmean(return_object=True)
        else:
            thex = x

        if self.dir == 'y':
            dat = thex.get_zonal_mean()
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

        # plot zonal statistics
        if dat.ndim == 1:
            self.ax.plot(dat,x.lat[:,0], label=label)
        elif dat.ndim == 2:
            for i in range(len(dat)):
                self.ax.plot(dat[i,:], x.lat[:,0], label='time='+str(i))
                self.ax.grid(b='on')

        self.ax.set_ylim(-90., 90.)

        if show_ylabel:
            self.ax.set_ylabel('latitude [deg]')
        else:
            self.ax.set_yticks([])

        if xlim != None:
            self.ax.set_xlim(xlim)

        self.ax.grid(b='on')


#-----------------------------------------------------------------------

class GlecklerPlot(object):
    """
    Class to generate a Portraet plot to illustrate multi-model,
    multi-variable scores.

    It was introdcued by Gleckler et al (2008)

    References
    ----------
    [1] ﻿Gleckler, P.J., Taylor, K.E. & Doutriaux, C., 2008.
             Performance metrics for climate models. Journal of
             Geophysical Research, 113(D6).

    Example
    -------
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

    def __init__(self, fig=None, figsize=(8,6)):
        """
        constructor of C{GlecklerPlot}

        @param fig: figure to which to plot to. If None, then a new figure will be generated
        @type fig: matplotlib figure
        """
        if fig is None:
            color = 'grey'
            fig = plt.figure(facecolor=color,edgecolor=color,figsize=figsize)
        self.fig = fig

        self.models = []
        self.variables = []
        self.data = {}  # store data for plot
        self.pos = {}  # store position of plot
        self.ndigits = 2  # number of digits for number plotting


    def add_model(self, label):
        """
        register a model in the class
        @param label: string of the model
        @type label: str
        """
        s = label.replace(' ','_')
        if s not in self.models:
            self.models.append(s)

    def _get_model_sequence(self):
        seq = {}
        cnt = 1
        for m in self.models:
            seq.update({m:cnt})
            cnt += 1
        return seq

    def _model2short_label(self, m):
        """
        given a model key a short label for plotting is returned
        """
        return str(self._get_model_sequence()[m])

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
        ax.set_xticks([])
        ax.set_yticks([])

    def __value2color(self,v):
        """
        return a color based on the value given
        the information on the colormap and its
        normalization is used for that purpose

        @param v: value of data
        @type v: float
        """
        if np.isscalar(v):
            r = self.cmap(self.norm(np.asarray([v])))
        else:
            r =  self.cmap(self.norm(v))
        return r.flatten()

    def __plot_triangle(self, ax, value, pos='top'):
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
        if value is None:
            return

        color = self.__value2color(value)
        pmax = int(max(self.pos.values()))

        if pmax > 4:
            raise ValueError('Only up to 4 observations supported!')

        if pmax > 2:
            # plot 4 triangles
            if pos == 'top':
                x = [0.,1.,0.5]
                y = [1.,1.,0.5]
                tpos = (0.5,0.75)
            elif pos == 'bottom':
                x = [0.,0.5,1.]
                y = [0.,0.5,0.]
                tpos = (0.5,0.25)
            elif pos == 'left':
                x = [0.,0.,0.5]
                y = [0.,1.,0.5]
                tpos = (0.25,0.5)
            elif pos == 'right':
                x = [1.,0.5,1.]
                y = [0.,0.5,1.]
                tpos = (0.75,0.5)
            else:
                print pos
                raise ValueError('Invalid position for plot')
        else:
            # plot only two triangles (diagonal)
            if pos == 'top':
                x = [0.,0.,1.]
                y = [0.,1.,1.]
                tpos = (0.25,0.75)
            elif pos == 'bottom':
                x = [1.,1.,0.]
                y = [1.,0.,0.]
                tpos = (0.75, 0.25)
            else:
                print 'Positions: ', self.pos.values()
                raise ValueError, 'Invalid position for plot: ' + str(pos) + ' pmax: ' + str(pmax)

        xy = list(zip(x,y))
        p = Polygon(xy,edgecolor='white',linewidth=1,fill=True,linestyle='solid',facecolor=color)
        ax.add_patch(p)

        #--- add value as text if required
        if self.show_value:
            if self.labelthreshold is None:
                labelcolor=self.labelcolor
            else:
                if np.abs(value) >= self.labelthreshold:
                    labelcolor = self.labelcolor
                else:
                    labelcolor = 'black'
            ax.text(tpos[0],tpos[1],str(np.round(value,self.ndigits)),fontdict={'size':6,'color':labelcolor},horizontalalignment='center',verticalalignment='center')


#-----------------------------------------------------------------------

    def _normalize_data(self, method='median'):
        """
        calculate for each observational data set
        the relative deviation from the average

        @param method: specifies which method should be used for calculation of the "mean" model. The paper of Gleckler et al. (2008)
                       uses the median, that's why it is the default. another option is to use the 'mean'
                       ['median','mean']
        @type method: str
        """
        pos = np.unique(self.pos.values())
        if hasattr(self,'_raw_data'):
            # data has been already normalized; take original data
            self.data = copy.deepcopy(self._raw_data)
        else:
            # preserve original calculated data
            self._raw_data = copy.deepcopy(self.data)

        for var in self.variables:
            for p in pos:
                # calculate multimodel mean/median
                xm = self._get_mean_value(p,var,method=method)
                for k in self.data:
                    if (self.pos[k] == p) & ('_' + var + '_' in k):
                        self.data[k] = (self.data[k] - xm) / xm #see Glecker et al, eq.2

#-----------------------------------------------------------------------

    def __draw_ranking_scatter(self, p1, p2, var, ax=None, marker='o', color='red', show_text=False):
        """
        ranking scatterplot between two positions
        """

        def _pos2label(p):
            if p == 1:
                return 'top'
            elif p == 2:
                return 'bottom'
            elif p == 3:
                return 'left'
            elif p == 4:
                return 'right'

        r1=self._get_model_ranking(p1, var)  # ranking for first position
        r2=self._get_model_ranking(p2, var)  # ranking for second position

        if len(r1) == 0:
            return None
        if len(r2) == 0:
            return None

        # generate two array where each specifies the rank of a particular model
        x = np.arange(len(r1))+1    #pos1
        y = np.ones_like(x)*np.nan  #pos 2
        for i in xrange(len(r1)):
            for j in xrange(len(r1)):
                if r1[i] == r2[j]:
                    y[i] = j+1

        # Spearman correlation coefficient based on ranks
        spear = np.corrcoef(np.asarray(x), np.asarray(y))[0][1]

        #--- generate plot ---
        if ax is None:
            f = pl.figure()
            ax = f.add_subplot(1,2,1)
        else:
            f = ax.figure

        ax.plot(x,y,marker=marker,color=color,label=_pos2label(p1) + ' vs. ' + _pos2label(p2) + ' ($r_s$=' + str(round(spear,2)) + ')' ,linestyle='None')
        if show_text:
            for i in xrange(len(x)):
                xy = (x[i],y[i])
                label = self._model2short_label(r1[i])
                ax.annotate(label,xy,color=color)

        return ax

#-----------------------------------------------------------------------

    def _draw_error_scatter(self, p1, p2, var, color='red', marker='*', ax=None):
        """
        draw scatterplot of errors for a combination
        of two different observations p1,p2
        """

        def _pos2label(p):
            if p == 1:
                return 'top'
            elif p == 2:
                return 'bottom'
            elif p == 3:
                return 'left'
            elif p == 4:
                return 'right'

        if ax is None:
            f = pl.figure()
            ax = f.add_subplot(111)
        d1 = self._get_model_error(p1,var) #this gives a dict with (key,value)
        d2 = self._get_model_error(p2,var)

        if len(d1) == 0:
            return None
        if len(d2) == 0:
            return None

        x = []
        y=[]
        labels=[]
        for k in d1.keys():
            labels.append(k)
            x.append(d1[k])
            #find corresponding model in second observation
            for l in d2.keys():
                if l == k:
                    y.append(d2[k])
        x = np.asarray(x)
        y = np.asarray(y)

        if len(x) != len(y):
            print(d1)
            print(d2)
            raise ValueError('The length of the arrays are not the same !')

        slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(x,y)
        if r_value is None:
            return ax
        else:
            ax.plot(x,y,marker=marker,color=color,label=_pos2label(p1) + ' vs. ' + _pos2label(p2) + ' ($r$=' + str(round(r_value, 2)) + ')' ,linestyle='None')
        return ax

    def plot_model_error(self, var):
        """
        plots a scatterplot of errors of different models
        for a particular variable

        Parameters
        ----------
        var : str
            variable to be investigated
        """

        fig = pl.figure()
        gs = gridspec.GridSpec(1, 2, wspace=0.05, hspace=0.05, bottom=0.2, width_ratios = [3,1])
        ax = fig.add_subplot(gs[0])

        # 1 vs. 2
        self._draw_error_scatter(1, 2, var, color='red', marker='o',ax=ax)

        # 1 vs. 3
        self._draw_error_scatter(1, 3, var, color='green', marker='*',ax=ax)

        # 1 vs. 4
        self._draw_error_scatter(1, 4, var, color='blue', marker='^',ax=ax)

        # 2 vs. 3
        self._draw_error_scatter(2, 3, var, color='grey', marker='x', ax=ax)

        # 2 vs 4
        self._draw_error_scatter(2, 4, var, color='m', marker='+', ax=ax)

        # 3 vs 4
        self._draw_error_scatter(3, 4, var, color='c', marker='h', ax=ax)

        if ax is not None:
            ax.legend(prop={'size':8}, ncol=1, fancybox=True, loc='upper left')
            ax.set_xlabel('$\epsilon$ (observation X)')
            ax.set_ylabel('$\epsilon$ (observation Y)')

            xmi,xma = ax.get_xlim()
            ymi,yma = ax.get_ylim()

            ax.set_ylim(min(xmi, ymi),max(xma, yma))
            ax.set_xlim(min(xmi, ymi),max(xma, yma))
            ax.grid()
            ax.set_title('Comparison of model errors: ' + var.upper())
            ax.plot(ax.get_xlim(),ax.get_xlim(),'k--') #1:1 line
        return fig

    def write_ranking_table(self, var, filename, fmt='latex'):
        """
        write results of model ranking to an ASCII table

        Parameters
        ----------
        var : str
            name of variable to analyze
        filename : str
            name of file to write table to
        fmt : str
            specifies output format for table ['latex','markdown']
        """

        if fmt not in ['latex', 'markdown']:
            raise ValueError('Invalid output format')

        if fmt == 'latex':
            sep = ' & '
            eol = ' \\\  \n'
            sol = '            '
        elif fmt == 'markdown':
            sep = ' | '
            eol = ' | \n'
            sol = ''
        else:
            raise ValueError('Unrecognized output format')

        if os.path.exists(filename):
            os.remove(filename)

        def _rnk2str(r):
            if r <= 3:
                return '{\\bf ' + str(r) + '}'
            else:
                return str(r)

        def _get_model_rank(k, x):
            """
            returns rank of model given a key k and an ordered list x

            Returns
            --------
            rank of key in list; NOTE that this is NOT the index!
            """
            if k not in x:
                if len(x)==0:  # no observational data available
                    return '-'
                else:
                    # here observations are there, but model key is not in
                    return '-'
            for i in xrange(len(x)):
                if x[i] == k:
                    return i+1
            raise ValueError('Fatal error: no solution found!')

        # get for each obs. dataset a list of models keys which are ordered
        r1 = self._get_model_ranking(1, var)  # returns model keys
        r2 = self._get_model_ranking(2, var)
        r3 = self._get_model_ranking(3, var)
        r4 = self._get_model_ranking(4, var)

        # now write a table with different columns for each dataset
        o = open(filename, 'w')
        if fmt == 'latex':
            o.write('        \\begin{tabular}{lcccc} \n')
            o.write(sol + '\\hline \n')
            s = sol + 'model' + sep + 'obs-top' + sep + 'obs-bott' + sep + 'obs-left' + sep + 'obs-right' + eol
            o.write(s)
            o.write(sol + '\\hline \n')
        for m in self.models:
            rnk1 = _get_model_rank(m, r1)
            rnk2 = _get_model_rank(m, r2)
            rnk3 = _get_model_rank(m, r3)
            rnk4 = _get_model_rank(m, r4)

            rnk1 = _rnk2str(rnk1)
            rnk2 = _rnk2str(rnk2)
            rnk3 = _rnk2str(rnk3)
            rnk4 = _rnk2str(rnk4)

            s = sol + m.replace('_','-') + sep + rnk1 + sep + rnk2 + sep + rnk3 + sep + rnk4 + eol
            o.write(s)
        if fmt == 'latex':
            o.write(sol + '\\hline \n')
            o.write('        \end{tabular}')
        o.close()

    def plot_model_ranking(self, var, show_text=False):
        """
        plots a model ranking scatterplot, indicating
        if models have similar ranks between different observational
        datasets

        in case that a certain combination does not exist, the
        corresponding plot is simply not generated.

        Parameters
        ----------
        var : str
            name of variable to analyze
        show_text : bool
            annotate plot using text for models as labels
        """

        tmp=self._get_model_ranking(1, var)

        fig = pl.figure()
        gs = gridspec.GridSpec(1, 2, wspace=0.05, hspace=0.05, bottom=0.2, width_ratios = [3,1])
        ax = fig.add_subplot(gs[0])

        # 1 vs. 2
        self.__draw_ranking_scatter(1,2,var,color='red',marker='o',show_text=show_text,ax=ax)
        # 1 vs. 3
        self.__draw_ranking_scatter(1,3,var,color='green',marker='*',ax=ax,show_text=show_text)
        # 1 vs. 4
        self.__draw_ranking_scatter(1,4,var,color='blue',marker='^',ax=ax,show_text=show_text)
        # 2 vs. 3
        self.__draw_ranking_scatter(2,3,var,color='grey',marker='x',ax=ax,show_text=show_text)
        # 2 vs 4
        self.__draw_ranking_scatter(2,4,var,color='m',marker='+',ax=ax,show_text=show_text)
        # 3 vs 4
        self.__draw_ranking_scatter(3,4,var,color='c',marker='h',ax=ax,show_text=show_text)

        if ax is not None:
            ax.legend(prop={'size':8}, ncol=1,fancybox=True, loc='upper left')
            ax.set_xlabel('rank(observation X)')
            ax.set_ylabel('rank(observation Y)')
            ax.set_ylim(ymin=0, ymax=len(tmp)+1)
            ax.set_xlim(xmin=0, xmax=len(tmp)+1)
            ax.grid()
            ax.set_title('Comparison of model ranking: ' + var.upper())
            ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--')  # 1:1 line

        ax2 = fig.add_subplot(gs[1])

        dy = 0.1
        yoff = dy
        for k in tmp:
            ax2.text(0.1, yoff, self._model2short_label(k) + ': ' + k)
            yoff += dy
        ax2.set_ylim(0., yoff)
        ax2.set_xticks([])
        ax2.set_yticks([])

        return fig

#-----------------------------------------------------------------------

    def _get_model_error(self,pos,var):
        """
        get error of each model for a certain variable and observation
        NOTE: to obtain a relative model ranking, one needs to normalize the data before, otherwise the absolute values
              are used!
        """
        x = []; keys=[]
        for k in self.pos:
            if (self.pos[k] == pos) & ('_' + var + '_' in k):
                x.append(self.data[k])
                keys.append(k[:k.index(var)-1]) #model name

        x    = np.asarray(x)
        keys = np.asarray(keys)

        return dict(zip(keys, x))


#-----------------------------------------------------------------------


    def _get_model_ranking(self, pos, var):
        """
        get ranking of each model for a certain variable and observation
        NOTE: to obtain a relative model ranking, one needs to
        normalize the data before, otherwise the absolute values
        are used!
        """
        x = []
        keys=[]
        for k in self.pos:
            if (self.pos[k] == pos) & ('_' + var + '_' in k):
                x.append(self.data[k])
                keys.append(k[:k.index(var)-1])  # model name

        x    = np.asarray(x)
        keys = np.asarray(keys)
        idx  = x.argsort()
        rnk  = np.arange(len(x))

        return keys[idx] #return list with keys which give ranked sequence


#-----------------------------------------------------------------------

    def _get_mean_value(self,pos,var,method='median'):
        """
        calculate mean value for a given observational dataset

        @param pos: position marker
        @type pos: int

        @param var: name of variable to analyze
        @type var: str

        @param method: specifies which method should be used for calculation of the "mean" model. The paper of Gleckler et al. (2008)
                       uses the median, that's why it is the default. another option is to use the 'mean'
                       ['median','mean']
        @type method: str
        """
        x = []
        for k in self.pos:
            if (self.pos[k] == pos) & ('_' + var + '_' in k):
                x.append(self.data[k])
        x = np.asarray(x)

        if method == 'median':
            return np.median(x)   #todo unittest for this!
        elif method == 'mean':
            return x.mean()
        else:
            raise ValueError, 'Invalid option in _get_mean_value() ' + method






#-----------------------------------------------------------------------

    def plot(self, cmap_name='RdBu_r', vmin=-1.0, vmax=1.0, nclasses=15,
             normalize=True, size=10, method='median', title=None,
             show_value=False, logscale=False, labelcolor='black',
             labelthreshold=None, cmap=None, norm=None,
             colorbar_boundaries=None, show_colorbar=True,
             autoscale=True, ticks=None):
        """
        plot Gleckler diagram

        @param cmap_name: name of the colormap to be used
        @type cmap_name: str

        @param vmin: lower boundary for plotting
        @type vmin: float

        @param vmax: upper boundary for plotting
        @type vmax: float

        @param nclasses: number of classes for colormap
        @type nclasses: int

        @param size: size of labels
        @type size: int

        @param normalize: normalize data relative to multimodel mean (affects self.data)
        @type normalize: bool

        @param method: specifies which method should be used for calculation of the "mean" model. The paper of Gleckler et al. (2008)
                       uses the median, that's why it is the default. another option is to use the 'mean'
                       ['median','mean']
        @type method: str

        @param logscale: reformats labels of colorbar to log with base 10. NOTE, this does NOT change the data! The data needs to be added in log10() scale already in function add_data() !!!
        @type logscale: bool

        @param labelcolor: color of text labels
        @type labelcolor: str

        @param labelthreshold: allows for plotting of labels in different colors. if abs(data)>= labelthreshold, then the label is plotted in labelcolor, else it is plotted in black
        @type labelthreshold: float

        @param show_colorbar: show colorbar plot
        @type show_colorbar: bool

        @param autoscale: autoscale figure size
        @type autoscale: bool

        """


        self.show_value=show_value
        self.labelcolor=labelcolor
        self.labelthreshold=labelthreshold
        self.bounds = colorbar_boundaries

        if normalize:
            self._normalize_data(method=method)

        nm = len(self.models); nv = len(self.variables)
        if nm == 0:
            ax = self.fig.add_subplot(111,frameon=True,aspect='equal',axisbg='grey')
            ax.text(0.5, 0.5,'no plot possible (missing model data)',
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform = ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            print '    Gleckler plot can not be generated, as no model data available!'
            return


        #- colormap
        if cmap is None:
            self.cmap = plt.cm.get_cmap(cmap_name, nclasses)
            self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        else: #use provided colormap
            self.cmap = cmap
            if norm is None:
                raise ValueError, 'If a colormap is provided, the NORM function shall also be provided!'
            else:
                self.norm = norm
                #self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        gs = gridspec.GridSpec(nm, nv, wspace=0.05,hspace=0.05,bottom=0.2) #generate grid for subplots

        cnt = 0
        cnt_m = 0

#        model_list = self.models.sort()

        for model in self.models:
            cnt_m += 1
            cnt_v = 0
            for variable in self.variables:
                ax = self.fig.add_subplot(gs[cnt], frameon=True, aspect='equal', axisbg='grey')
                self.__set_ax_prop(ax)

                # labels
                if cnt_v == 0:
                    ax.set_ylabel(model,size=size,rotation='horizontal',horizontalalignment='right') #set labels for models
                if cnt_m == nm:
                    ax.set_xlabel(variable,size=size)

                self.__plot_triangle(ax,self.get_data(variable,model,1), pos='top')    # upper triangle
                self.__plot_triangle(ax,self.get_data(variable,model,2), pos='bottom') # lower triangle
                #~ if variable == 'albedo':
                    #~ print 'POS 2: ', self.get_data(variable,model,2)
                self.__plot_triangle(ax,self.get_data(variable,model,3), pos='left')   # left triangle
                self.__plot_triangle(ax,self.get_data(variable,model,4), pos='right')  # right triangle
                cnt += 1
                cnt_v += 1

        #--- legend
        #- get positions of subplots to determine optimal position for legend
        def get_subplot_boundaries(g,f):
            x = g.get_grid_positions(f)
            b = x[0]
            t = x[1]
            l = x[2]
            r = x[3]
            return l[0], r[-1], b[-1], t[0]

        if autoscale:
            self.fig.set_size_inches(3.+0.5*nv,4.+0.5*nm) #important to have this *before* drawring the colorbar

        left,right,bottom,top = get_subplot_boundaries(gs,self.fig)
        # draw legend
        c=1.
        width=(right-left)*c
        if show_colorbar:
            self._draw_colorbar(left,width,logscale=logscale,ticks=ticks)

        if title is not None:
            self.fig.suptitle(title)

        return self.fig

#-----------------------------------------------------------------------

    def get_data(self, v, m, p):
        """
        return data for a particular model and variable

        v : str
            name of variable
        m : str
            model name
        p : int
            position
        """
        r = None
        k = self.__gen_key(m,v,p)
        if k in self.data.keys():
            r = self.data[k]
        else:
            #print 'NOT FOUND!!', k
            pass
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
        if m is None:
            return None
        if v is None:
            return None
        return m.replace(' ','_')+'_'+v.replace(' ','_')+'_'+str(int(p))

#-----------------------------------------------------------------------

    def add_data(self, v, m, x, pos=1):
        """
        add a data for plotting

        Parameters
        ----------
        v : str
            name of variable
        m : str
            model name
        x : float
            value to be plotted in Gleckler plot
        pos : int
            position where to plot data 1=top triangle, 2=lower triangle
        """
        if x is not None:
            if v in self.variables:
                if m in self.models:
                    self.data.update({ self.__gen_key(m,v,pos) : x})
                    self.pos.update({ self.__gen_key(m,v,pos) : pos})
                else:
                    pass
            else:
                pass
        else:
            pass

#-----------------------------------------------------------------------

    def calc_index(self, x, y, model, variable, weights=None, time_weighting=True):
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

        @param weights: weights to be applied to the data before index
        calculation; dedicated for spatial area weights
        @type weights: numpy array

        @return: returns performance index aggregated over time
        (NOTE: this is still E**2 !!!)
        @rtype float

        """

        if weights is None:
            #set weights according to cell area
            if x.cell_area is not None:
                weights = x._get_weighting_matrix()
            else:
                print('WARNING: no weights when calculating performance index')
                weights = np.ones(x.data.shape)
        else: #weights are given
            if x.cell_area is not None:
                print('WARNING: cell weights are given, while cell_area available from data!!')

        from diagnostic import Diagnostic
        D = Diagnostic(x, y=y)
        # reichler performance index
        # (might return a list if multiple times analyzed)
        e2 = D.calc_reichler_index(weights)  # returns [time, index]
        if e2 is not None:
            e2 = np.ma.masked_where(np.isnan(e2), e2)
        # Note that the returned value is E**2 for each timestep!
        # When doing the normalization, one needs to take the sqrt(E++2)
        # to obtain the actual RMSE

        if e2 is None:
            return None
        else:
            if time_weighting:
                days = np.asarray(x._days_per_month())  # number of days for each month
                wt = days / float(days.sum())  # weights for time
            else:
                wt = np.ones(len(e2)) / float(len(e2))

            return np.sqrt((e2*wt).sum())

#-----------------------------------------------------------------------

    def _draw_colorbar(self, left, width, logscale=False, ticks=None):
        """
        draw legend for Glecker plot. Requires information on
        the positioning of the colormap axis which can be obtained from

        left,right,bottom,top = get_subplot_boundaries(gs,self.fig)

        left : float
            left position of axis
        width : float
            width of the axis to plot colorbar
        """

        if logscale:
            #http://www.mailinglistarchive.com/html/matplotlib-users@lists.sourceforge.net/2010-06/msg00311.html
            from matplotlib.ticker import LogLocator, LogFormatter, FormatStrFormatter
            l_f = FormatStrFormatter('$10^{%d}$')
        else:
            l_f = None

        # left, bottom, width, height
        cax = self.fig.add_axes([left, 0.05, width, 0.05])
        if self.bounds is not None:
            extend = 'both'
        else:
            extend = 'neither'
        cb = mpl.colorbar.ColorbarBase(cax, cmap=self.cmap,
                                   norm=self.norm,
                                   orientation='horizontal',
                                   format=l_f,boundaries=self.bounds,
                                   extend=extend, ticks=ticks)

    def _draw_legend(self, labels, title=None):
        """
        draw automatic legend for Gleckler plot

        EXAMPLE:
        fl=G._draw_legend({1:'top',2:'bott',3:'left',4:'right'})

        CAUTION: when saving the figure, do NOT use bbox_inches='tight', as this might cut the labels!

        @param labels: dictionary as {position:'label'}; e.g. {1:'label1',2:'label2',3:'label3',4:'label4'}
        """

        if len(self.pos) < 1:
            print 'Legend can not be plotted for Gleckler, as no data available!'
            return

        pmax = max(self.pos.values())

        #generate separate figure for legend
        f=plt.figure()
        ax=f.add_subplot(111, frameon=True, aspect='equal', axisbg='grey')
        f.subplots_adjust(bottom=0.25, top=0.75, left=0.25, right=0.75)

        for k in labels.keys():
            if k == 1:
                pos = 'top'
            elif k == 2:
                pos = 'bottom'
            elif k == 3:
                pos = 'left'
            elif k == 4:
                pos = 'right'
            else:
                raise ValueError, 'Can not draw Gleckler legend! Invalid position value! ' + str(k)

            oldval = self.show_value
            self.show_value=False
            self.__plot_triangle(ax, np.random.random(), pos=pos)
            self.show_value=oldval
            ax.set_xticks([])
            ax.set_yticks([])

        fontsize=16
        linewidth=3

        for k in labels.keys():
            if k == 1: #top
                ax.annotate(labels[k], xy=(0.5, 0.9), xycoords='axes fraction', xytext=(0., 1.2), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="angle3,angleA=0,angleB=-90",linewidth=linewidth),horizontalalignment='left',size=fontsize)
            elif k == 2:
                ax.annotate(labels[k], xy=(0.5, 0.1), xycoords='axes fraction', xytext=(0., -0.3), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="angle3,angleA=0,angleB=-90",linewidth=linewidth),horizontalalignment='left',size=fontsize)
            elif k == 3:
                ax.annotate(labels[k], xy=(0.1, 0.5), xycoords='axes fraction', xytext=(-0.6,0.2), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="angle3,angleA=0,angleB=-90",linewidth=linewidth),horizontalalignment='left',size=fontsize)
            elif k == 4:
                ax.annotate(labels[k], xy=(0.9, 0.5), xycoords='axes fraction', xytext=(1.1,0.8), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",connectionstyle="angle3,angleA=0,angleB=-90",linewidth=linewidth),horizontalalignment='left',size=fontsize)

        if title is not None:
            f.suptitle(title, size=fontsize)

        return f

#-----------------------------------------------------------------------

def __basemap_ancillary(m,latvalues = None, lonvalues = None,drawparallels=True,drawcountries=True,land_color=0.8):
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

    if latvalues is None:
        latvalues=np.arange(-90.,120.,30.)
    if lonvalues is None:
        lonvalues= np.arange(-180.,180.,90.)
    if drawcountries:
        m.drawcountries()
    m.drawcoastlines()
    m.drawlsmask(lakes=True,land_color=land_color)
    m.drawmapboundary() # draw a line around the map region
    if drawparallels:
        m.drawparallels(latvalues,labels=[1, 0, 0, 0])
        m.drawmeridians(lonvalues,labels=[0, 0, 0, 1]) # draw meridians

#-----------------------------------------------------------------------

def pm_bar(x, y=None, pcolor='red', ncolor='blue', ax=None, **kwargs):
    """
    generate a nice looking barchart with different color for positive/negative numbers

    @param x: x-variable or variable to plot (if y is not given)
    @param y: y-variable (optional)
    @param ax: axis handle
    @return: returns handle axis
    """

    if ax is None:
        f = pl.figure()
        ax = f.add_subplot(111)
    else:
        ax = ax

    if y is None:
        y = x*1.
        x = np.arange(len(y))
    else:
        pass

    yp = y*1.
    yp[y<0.] = 0.
    yn = y*1.
    yn[y>0.] = 0.

    #--- plot
    ax.bar(x, yp, color=pcolor, edgecolor='None', **kwargs)
    ax.bar(x, yn, color=ncolor, edgecolor='None', **kwargs)

    return ax

#-----------------------------------------------------------------------

def map_season(x, figsize=(8,6), **kwargs):
    """
    generate a seasonal plot
    all arguments are parsed directly to map_plot function

    if kwargs contain a 'figure' argument, then this figure fill be used
    for plotting. Otherwise a new figure will be generated

    Note, that it is not checked, if the seasonal mean values were
    precalculated correctly.
    It is ASSUMED that the seasonal means are calculated using cdo, which leads to

    a) yseasmean --> DJF,MAM,JJA,SON where the timestamp in the file corresponds to the LAST valid data used
    b) ymonmean --> monthly means

    Parameters
    ----------
    x : Data
        data to be plotted
    figsize : tuple
        specifies size of figure (see pyplot documentation)

    Returns
    -------
    returns a figure handler
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

    if 'vmin' not in kwargs.keys():
        raise ValueError, 'vmin argument is obligatory for map_seasons()'
    if 'vmax' not in kwargs.keys():
        raise ValueError, 'vmax argument is obligatory for map_seasons()'

    if kwargs['vmin'] is None:
        raise ValueError, 'vmin MUST NOT be None!'
    if kwargs['vmax'] is None:
        raise ValueError, 'vmax MUST NOT be None!'

    #/// figure and axes
    if 'figure' in kwargs:
        f = kwargs['figure']
    else:
        f = pl.figure(figsize=figsize)

    if 'title' in kwargs:
        tit = kwargs.pop('title')
    else:
        tit = x.label

    if 'drawparallels' in kwargs:
        drawparallels = kwargs.pop('drawparallels')
    else:
        drawparallels = False

    if 'savefile' in kwargs:
        savefile = kwargs.pop('savefile')
        if '.nc' in savefile:
            savefile=savefile[:-3]
    else:
        savefile = None

    # plot
    if year:
        labels=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    else:
        labels=['DJF','MAM','JJA','SON']

    # check dates
    if year:
        mo = 1
        for t in x.time:
            if x.num2date(t).month != mo:
                print x.num2date(t), mo
                raise ValueError, 'Invalid monthly sequence! Can not plot results!'
            mo +=1

    #/// in case that an overlay is provided, this needs to be processed for each timestep individually
    if 'overlay' in kwargs.keys():
        overlays = kwargs.pop('overlay')
    else:
        overlays = None

    for i in range(nvals):
        if year:
            ax = f.add_subplot(4,3,i+1)
            #if i % 3 == 2:
            if i > 8:
                show_colorbar = True
            else:
                show_colorbar=False
        else:
            ax = f.add_subplot(2,2,i+1)
            if 'show_colorbar' in kwargs:
                show_colorbar = kwargs.pop('show_colorbar')
            else:
                show_colorbar = True

        d = x.copy()
        d.data = x.data[i,:,:]
        d.label = labels[i]

        if overlays is None:
            overlay = None
        else:
            overlay = overlays[i,:,:]

        if savefile is not None:
            tmpoutname=savefile + '_' + labels[i]
        else:
            tmpoutname = None

        map_plot(d, ax=ax, show_colorbar=show_colorbar, overlay=overlay,
                 savefile=tmpoutname, colorbar_orientation='horizontal',
                 drawparallels=drawparallels, **kwargs)
        del d
    f.suptitle(tit,size=16)
    return f

#-----------------------------------------------------------------------


def map_plot(x, use_basemap=False, ax=None, cticks=None, region=None,
             nclasses=10, cmap_data='jet',
             title=None, regions_to_plot=None, logplot=False,
             logoffset=None, show_stat=False,
             f_kdtree=False, show_colorbar=True, latvalues=None,
             lonvalues=None, show_zonal=False,
             zonal_timmean=True, show_timeseries=False,
             scal_timeseries=1., vmin_zonal=None, vmax_zonal=None,
             bluemarble = False, contours=False, overlay=None,
             titlefontsize=14, drawparallels=True, drawcountries=True,
             show_histogram=False,
             contourf = False, land_color=(0.8,0.8,0.8),
             regionlinewidth=1, bins=10,
             colorbar_orientation='vertical', stat_type='mean',
             cax_rotation=0., cticklabels=None, proj='robin',
             plot_method='colormesh', boundinglat=60.,
             savefile=None, lon_0=0., lat_0=0., savegraphicfile=None,
             **kwargs):
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

    @param bins: bins for histogram
    @type bins: int

    @param colorbar_orientation: specifies if colorbar shall be vertical or horizontal aligned ['horizontal','vertical']
    @type colorbar_orientation: str

    @param stat_type: specifies if mean or median shall be used for statistics ['mean','median']
    @type stat_type: str

    @param cax_rotation: rotation of labels for colorbar axis
    @type cax_rotation: float

    @param cticklabels Labels for the ticks of the colorbar
    @type cticklabels: list of str labels

    @param proj: Basemap projection parameter string, specifying the type of projections to be used ['robin','npstere']
    @type proj: str

    @param plot_method: specifies method how data shall be plotted (applies only when using Basemap)
                        allowed options are ['colormesh','scatter']
    @type plot_method: str


    """

    def _get_unstructured_collection(vlon,vlat,xm,vmin,vmax,basemap_object=None):
        """
        get collection of patches for plotting of unstructured grid

        vlon,vlat : vertex lon/lat [ncells,npoints]
        xm: mean field, generated by timmean function
        """

        #init
        Path = mpath.Path
        patches = []
        pdata = xm[0,:]*1. #full list of data
        vmsk  = np.ones_like(pdata).astype('bool') #mask to indicate which cells contain valid data

        for i in xrange(x.ncell):
            if np.any(vlon[i,:]) > 180.: #todo fix this properly !!!
                vmsk[i] = False
                continue
            if basemap_object is None:
                xv=vlon[i,:]
                yv=vlat[i,:]
            else:
                xv,yv = basemap_object(vlon[i,:], vlat[i,:])    #todo: how to properly deal with boundary problem ????
            if (vlon[i,:].min() < -100.) & (vlon[i,:].max() > 100.): #todo
                #... triangles across the boundaries of the projection are a problem
                # ... solution: generate two triangles ! TODO
                vmsk[i] = False
                continue

            verts=np.asarray([xv,yv]).T

            #--- specify how vertices are interconnected (here simple connection by lines)
            codes = [Path.MOVETO,Path.LINETO,Path.LINETO]

            #--- construct object and append to library of objects ---
            path = mpath.Path(verts, codes,closed=True)
            patches.append(mpatches.PathPatch(path))

        pdata = np.asarray(pdata)

        if vmin is None:
            vmin=pdata.min()
        if vmax is None:
            vmax=pdata.max()

        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        collection = PatchCollection(patches, cmap=cmap,norm=norm, alpha=1.,match_original=False,edgecolors='grey')  #construct library of all objects
        collection.set_array(pdata[vmsk]) #assign data values here

        return collection

    if savegraphicfile is not None:
        savegraphicfile = savegraphicfile.replace(' ', '_')


    if 'vmin' in kwargs.keys():
        vmin = kwargs['vmin']
    else:
        vmin = None
    if 'vmax' in kwargs.keys():
        vmax = kwargs['vmax']
    else:
        vmax = None

    if plot_method not in ['colormesh','scatter']:
        raise ValueError, 'Invalid plotting option ' + plot_method

    # checks
    if proj not in ['robin', 'npstere', 'spstere']:
        raise ValueError('ERROR: projection type not validated for map_plot so far: %s' % proj)
    if proj == 'npstere': #todo: for stereographic projection, scatter is used as method at the moment
        plot_method = 'scatter'
    if proj == 'spstere': #todo: for stereographic projection, scatter is used as method at the moment
        plot_method = 'scatter'

    if overlay is not None:
        if x.data.ndim == 2:
            if overlay.shape != x.data.shape:
                print overlay.shape, x.data.shape
                raise ValueError('Invalid geometry for overlay !')
        elif x.data.ndim == 3:
            if overlay.shape != x.data[0,:,:].shape:
                print overlay.shape, x.data.shape
                raise ValueError('Invalid geometry for overlay !')
        else:
            raise ValueError('Overlay for this geometry not supported!')

    #--- create new figure
    if ax is None:
        fig = plt.figure()

        # with timeseries plot?
        if show_timeseries:
            gs = gridspec.GridSpec(2, 1, wspace=0.05, hspace=0.05, bottom=0.2, height_ratios = [5,1])
            ax = fig.add_subplot(gs[0])
            ax2 = fig.add_subplot(gs[1])
        else:
            ax = fig.add_subplot(111)
    else:
        fig = ax.figure
        if show_timeseries:
            raise ValueError, 'Showing timeseries and providing some prior axis is currently not impelmented!' #todo

    # if cmap provided in kwargs, then remove it and set cmap_data
    kwargs1 = kwargs.copy()
    if 'cmap' in kwargs:
        cmap_data = kwargs1.pop('cmap')
    if ('levels' in kwargs) and (contours == False): #levels not needed
        dummy = kwargs1.pop('levels')

    #--- create colormap
    if hasattr(cmap_data, 'monochrome'):
        # colormap object was given
        cmap = cmap_data
    else:
        cmap = plt.cm.get_cmap(cmap_data, nclasses)

    # temporal mean fields as data to plot
    xm = x.timmean() #returns an array

    # logscale plot ?
    if logplot:
        if logoffset is None:
            if xm.min() < 0.: logoffset = abs(xm.min())*1.01
            else: logoffset = 0.
        else: logoffset = logoffset
        print '     logoffset: ', logoffset
        xm = np.log10( xm + logoffset )

    #--- save field that is plotted as file
    if savefile is not None:
        if savefile[:-3] != '.nc':
            savefile += '.nc'
        tmp = x.copy()
        tmp.data = xm*1.
        tmp.time = None
        tmp.save(savefile, varname='temporal_mean_field')
        del tmp

    #--- set projection parameters
    if proj == 'robin': #todo: more flexible handling of projection parameters (dictionary ??)
        lon_0=lon_0
        lat_0=lat_0
    elif proj == 'npstere':
        lon_0 = lon_0
        lat_0 = lat_0
        boundinglat = boundinglat
    elif proj == 'spstere':
        lon_0 = lon_0
        lat_0 = lat_0
        boundinglat = -boundinglat
    else:
        raise ValueError('Unsupported projection in map_plot (unsupported means, that it was not tested yet)')

    #--- plot using basemap
    if use_basemap:
        llcrnrlon=None
        llcrnrlat=None
        urcrnrlon=None
        urcrnrlat=None

        #if a region is specfied, then the plotting boundaries are set
        if region !=None:
            if not hasattr(region, 'lonmin'):
                print 'WARNING map boundaries can not be set, as region ' + region.label.upper() + ' has not lat/lon information'
            else:
                dlat = (region.latmax-region.latmin)*0.25
                dlon = (region.lonmax-region.lonmin)*0.25
                di = 0. # with 0 it works; for other values problems may occur for negative lon!
                llcrnrlon=region.lonmin - di; llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di; urcrnrlat=region.latmax + di
                proj='tmerc' #use mercator projection at regional scale as robinson does not work!

        ############################################
        # generate Basemap map
        ############################################
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax,
                   llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                   urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                   boundinglat=boundinglat)

        if bluemarble:
            m1.bluemarble()

        if x.gridtype == 'unstructured':
            #support of unstructured grid types, like e.g. provided by ICON model
            #it assumes that the data object has a list of center coordinates which correspond to the data
            #and vlon/vlat attributes with corresponding vertices corresponding to the center coordinates

            if not hasattr(x, 'vlon'):
                raise ValueError, 'Plotting for unstructured grid not possible, as VLON attribute missing!'
            if not hasattr(x, 'vlat'):
                raise ValueError, 'Plotting for unstructured grid not possible, as VLAT attribute missing!'

            #--- generate collection of patches for Basemap plot
            collection = _get_unstructured_collection(x.vlon, x.vlat, xm, vmin, vmax, basemap_object=m1)

        else: #unstructured gridtype


            #check if all longitudes are the same. If so, then a plotting with different options is possible.
            #otherwise, f_kdtree is activates as default option to ensure valid plotting
            #f_kdtree_act = f_kdtree
            #if x._lon360:
            #    if x._equal_lon():
            #        f_kdtree_act = f_kdtree
            #    else:
            #        print 'WARNING: longitudes are not all the same! f_kdtree=True is therefore used!'
            #        f_kdtree_act = True


            if f_kdtree:
                #use KDTRee nearest neighbor resampling to avoid stripes in plotting
                lons = np.unique(x.lon)
                lats = np.unique(x.lat)
                lons.sort(); lats.sort()
                TLON,TLAT = np.meshgrid(lons, lats)  #generate target coordinates
                XT,YT = m1(TLON,TLAT)
                X=XT.copy()
                Y=YT.copy()
                shape0 = np.shape(XT)
                XT.shape = (-1)
                YT.shape = (-1)  # ... vectorize them for inertpolation
                tree = KDTree(zip(XT,YT)) #generate tree from TARGET coordinates

                #prepare data and interpolate
                xmap,ymap = m1(x.lon,x.lat)
                xmap.shape = (-1)
                ymap.shape = (-1)
                pts  = zip(xmap,ymap) #generate points to interpolate from source data
                dist,idx = tree.query(pts,k=1)     #perform nearest neighbor interpolation (returns distance and indices)

                #- map data to output matrix for plotting
                Z = np.ones(shape0)*np.nan
                Z.shape = (-1) #generate target vector
                omask = np.ones(shape0).astype('bool')
                omask.shape = (-1)

                msk1 = xm.mask.copy()
                msk1.shape = (-1)
                omask[idx] = msk1

                #~ omask[dist != 0.] = True

                xm1 = xm.copy()
                xm1.shape = (-1)
                Z[idx]   = xm1 #assign data and reshape it and set generate masked array
                Z[omask] = np.nan
                Z = np.reshape(Z,shape0)
                Z = np.ma.array(Z,mask=np.isnan(Z))

            else: #f_kdtree --> not kdtree

                #/// in the following, we check if the longitudes are in the right order
                #    to allow for an appropirate plotting. Basemap assumes ascending order
                #    of longitudes. For global projections, problems occur if lon are given
                #    as 0...360 deg. It needs to be shifted then. Otherwise the data will produce
                #    strange stripes when plotting.
                #
                #    REFERENCES:
                #    * http://pl.digipedia.org/usenet/thread/15998/16891/
                if plot_method == 'colormesh':
                    print 'Projection: ', proj
                    if x._lon360:  # if lon 0 ... 360, then shift data
                        tmp_lon = x._get_unique_lon() #get unique longitudes
                        tmplon1 = tmp_lon.copy()
                        Z, tmp_lon = shiftgrid(180, xm, tmp_lon, start=False)
                        if overlay is not None:
                            overlay, nope = shiftgrid(180, overlay, tmplon1, start=False)
                        lon, lat = np.meshgrid(tmp_lon, np.arange(Z.shape[0]))
                        lat = x.lat
                    else:
                        print '*** WARNING: not lon360 not validated yet, try KDTREE option if stripes in plot ***'
                        lon = x.lon
                        lat=x.lat
                        Z = xm
                elif plot_method == 'scatter':
                        lon = x.lon
                        lat = x.lat
                        Z = xm
                else:
                    raise ValueError('Invalid option')

                X, Y = m1(lon, lat)

        #here is still a problem in the plotting over land; masking does not work properly,
        #while the data as such is o.k.!
        #~ im1=m1.pcolormesh(xmap,ymap,xm,cmap=cmap,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)

        if not bluemarble:

            if x.gridtype=='unstructured':
                #--- do actual plotting
                im1=m1.ax.add_collection(collection) #add unstructure grid plots (e.g. triangles)
            else:
                if contours:
                    if 'levels' in kwargs1.keys():
                        levels = kwargs1.pop('levels')
                    else:
                        raise ValueError('When plotting with contours, you need to specify the levels option (see contour documnetation)')
                    if contourf:
                        im1=m1.contourf(X, Y, Z,levels, cmap=cmap, **kwargs1)
                    else:
                        im1=m1.contour(X, Y, Z, levels, cmap=cmap, **kwargs1)
                        ax.clabel(im1, inline=1, fontsize=10)  # contour label

                else:
                    if plot_method == 'colormesh':
                        im1=m1.pcolormesh(X,Y,Z,cmap=cmap,**kwargs1) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
                    elif plot_method == 'scatter':
                        im1=m1.scatter(X,Y,c=Z,marker='8',edgecolor='None',cmap=cmap,**kwargs1)
                    else:
                        raise ValueError, 'Invalid plotting option! ' + plot_method

                if overlay != None:
                    #print 'Doing overlay plot! ...'
                    xcoordnew=X[overlay]
                    ycoordnew=Y[overlay]
                    m1.scatter(xcoordnew, ycoordnew, marker='x', s=50, c='white', edgecolors='white', linewidth=1)
                    #todo: there is still a problem that the coordinates are not properly aligned with the grid cells!!!

            #/// ANCILLARY PLOT FOR BASEMAP ///
            __basemap_ancillary(m1,latvalues=latvalues,lonvalues=lonvalues,drawparallels=drawparallels,drawcountries=drawcountries,land_color=land_color)

    else: #use_basemap = False

        if x.gridtype=='unstructured':
            #--- generate collection of patches for Basemap plot
            collection = _get_unstructured_collection(x.vlon,x.vlat,xm,vmin,vmax,basemap_object=None)
            im1 = ax.add_collection(collection)
            ax.set_xlim(max(x.vlon.min(), -180.), min(x.vlon.max(), 180.))
            ax.set_ylim(max(x.vlat.min(), -90.), min(x.vlat.max(), 90.))
        else:

            #- normal plots
            if contours:
                if contourf:
                    im1 = ax.contourf(xm, cmap=cmap, **kwargs1)
                else:
                    im1 = ax.contour(xm, cmap=cmap, **kwargs1)
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

        ax.set_xticks([])
        ax.set_yticks([])



    # set legend aligned with plot (nice looking)
    divider = make_axes_locatable(ax)
    if show_zonal:
        caxv = divider.new_horizontal(size="3%", pad=0.1, axes_class=maxes.Axes)
        caxh = divider.new_vertical(size="5%", pad=0.1, axes_class=maxes.Axes,pack_start=True)
        caxzonaldummy = divider.new_horizontal(size="15%", pad=0.1, axes_class=maxes.Axes, pack_start=True)
        # this is still not working properly !
    else:
        caxv = divider.new_horizontal(size="3%", pad=0.1, axes_class=maxes.Axes)
        caxh = divider.new_vertical(size="5%", pad=0.1, axes_class=maxes.Axes,pack_start=True)
    if colorbar_orientation == 'vertical':
        cax = caxv
        caxdummy = caxh
    elif colorbar_orientation == 'horizontal':
        cax = caxh
        caxdummy = caxv
    else:
        raise ValueError, 'Invalid option for colorbar! ' + colorbar_orientation

    ax.figure.add_axes(cax)



    vmin = im1.get_clim()[0]
    vmax = im1.get_clim()[1]
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # dummy axis to ensure equal spacing in multiple plots
    caxdummy.set_xticks([])
    caxdummy.set_yticks([])
    caxdummy.set_frame_on(False)

    if show_colorbar:
        #plot actual colorbar
        cb   = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks=cticks, orientation=colorbar_orientation)
        if cticklabels is not None:
          cb.set_ticklabels(cticklabels)
    else:
        cax.set_xticks([])
        cax.set_yticks([])
        cax.set_frame_on(False)

    #--- add a histogram below the map plot
    if show_histogram:
        add_histogram(ax,x,bins=bins)

    # Zonal plot
    if show_zonal:
        if x._latitudecheckok:
            add_zonal_plot(ax, x, timmean=zonal_timmean, vmin=vmin_zonal, vmax=vmax_zonal) #,vmin=im1.get_clim()[0],vmax=im1.get_clim()[1])
        else:
            print('WARNING: zonal plot not possible due to invalid latitude configurations')

    def _add_region_basemap(m, r, color='red', linewidth=1):
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
        lons = corners[:,0]
        lats=corners[:,1]
        x,y = m(lons,lats)
        xy = list(zip(x, y))
        mapboundary = Polygon(xy, edgecolor=color, linewidth=linewidth, fill=False, linestyle='dashed')
        m.ax.add_patch(mapboundary)

    def _add_region_standard(ax, r, color='red', linewidth=1):
        """
        plot region r on top of a normal map plot

        @param m: map
        @type m: C{Basemap} object

        @param r: region to plot
        @type r: C{Region}

        @param color: color to plot region
        @type color: str
        """
        corners = r.get_corners() #get list of corner coordinates
        corners = np.asarray(corners)
        xcoords = corners[:, 0]
        ycoords=corners[:, 1]
        xy = list(zip(xcoords, ycoords))
        mapboundary = Polygon(xy, edgecolor=color, linewidth=linewidth, fill=False, linestyle='dashed')
        ax.add_patch(mapboundary)

    # plot regions in the map ---
    if regions_to_plot is not None:
        if use_basemap:
            for region in regions_to_plot:
                if region.type=='latlon':
                    _add_region_basemap(m1, region, linewidth=regionlinewidth)
        else:
            for region in regions_to_plot:
                if region.type=='index':
                    _add_region_standard(ax, region, linewidth=regionlinewidth)

    # set title
    if title is None:
        title = x._get_label()
    else:
        pass

    #--- show field statistics in title ?
    # calculates first the temporal mean and results shown then
    # are the means and std of the temporal mean fields which are weighted
    # appropriately according to the cell area
    if show_stat:
        tmp_xm = x.timmean(return_object=True)  # from temporal mean
        if stat_type == 'mean':
            me = tmp_xm.fldmean()
            st = tmp_xm.fldstd()
            assert(len(me) == 1)
            assert(len(st) == 1)
            me = me[0]
            st=st[0]
            atitle = 'mean: $' + str(round(me, 2))  + ' \pm ' + str(round(st,2)) + '$'
        elif stat_type == 'sum':  # area sum
            me = tmp_xm.areasum()
            assert(len(me) == 1)
            me = me[0]
            atitle = 'sum: $' + str(round(me, 2))  + '$'
        else:
            me = np.ma.median(tmp_xm.data)
            atitle = 'median: $' + str(round(me, 2)) + '$'
        ax.set_title(atitle, size=titlefontsize-2, loc='left')

    ax.set_title(title + '\n', size=titlefontsize, loc='center')
    ax.set_title(x._get_unit(), size=titlefontsize-2, loc='right')

    #/// show timeseries? ///
    if show_timeseries:
        ax2.plot(x.num2date(x.time),x.fldmean())
        ax2.grid()
        ax2.set_ylim(im1.get_clim()[0]*scal_timeseries,im1.get_clim()[1]*scal_timeseries)
        ti = ax2.get_yticks(); n=len(ti) / 2
        ax2.set_yticks([ti[0],ti[n],ti[-1]])
        ax2.set_ylabel(x._get_unit())

    if savegraphicfile is not None:
        # save graphics to file
        if os.path.exists(savegraphicfile):
            os.remove(savegraphicfile)
        fig.savefig(savegraphicfile, bbox_inches='tight', dpi=200)
    return fig

#-----------------------------------------------------------------------

def add_histogram(ax, x, bins=10):
    """
    add a histogram to an axis

    @param ax: axis to add the histogram plot to
    @type ax: axis

    @param x: Data that should be vizualized as histogram
    @type x: Data
    """

    divider = make_axes_locatable(ax)
    zax = divider.append_axes("bottom", "30%", pad=0.1)
    ax.figure.add_axes(zax, axisbg=ax.figure.get_facecolor())

    H = HistogrammPlot(ax=zax, bins=bins)
    H.plot(x)

#-----------------------------------------------------------------------

def add_zonal_plot(ax, x, timmean=True, vmin=None, vmax=None):
    """
    add a zonal plot to the axis.
    An area weigting is automaticall performed

    Parameters
    ----------
    ax : axis
        axis where zonal plot should be added to
    x : Data
        data to plot
    timmean : bool
        temporal mean for zonal plot [default=True]
    vmin : float
        minimum value for zonal plot
    vmax : float
        maximum value for zonal plot
    """

    divider = make_axes_locatable(ax)
    zax = divider.new_horizontal("15%", pad=0.1,
                                    axes_class=maxes.Axes,
                                    pack_start=True)
    ax.figure.add_axes(zax, axisbg=ax.figure.get_facecolor())

    ZP = ZonalPlot(ax=zax, dir='y')

    if x.data.ndim == 2:
        pass
    elif x.data.ndim == 3:
        nt,ny,nx = x.data.shape
    ZP.plot(x, timmean=timmean, show_ylabel=False)

    # set limits
    if ((vmin is None) & (vmax is None)):
        vmin = zax.get_xlim()[0]
        vmax = zax.get_xlim()[1]
        # symmetry if neg. and posiitve limits
        if (vmin < 0.) & (vmax>0.):
            val = max(abs(vmin),abs(vmax))
            vmin = -val; vmax = val

    if vmin is None:
        vmin = zax.get_xlim()[0]
    if vmax is None:
        vmax = zax.get_xlim()[1]
    zax.set_xlim(vmin, vmax)

    #set only first and last label
    zax.set_xticks([vmin, vmax])
    zax.plot([0,0], zax.get_ylim(), linestyle='-', color='grey')

    for tick in zax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)

    return zax

#-----------------------------------------------------------------------

def add_nice_legend(ax, im,cmap, cticks=None, dummy=False, fontsize=8, label=None):
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

    if label != None:
        cax.set_ylabel(label)


#-----------------------------------------------------------------------

def hov_difference(x,y,climits=None,dlimits=None,data_cmap='jet',nclasses=15,cticks=None,cticks_dif=None,ax1=None,ax2=None,ax3=None,rescaley=6,grid=True,rescalex=1,clabel=None,**kwargs):
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

    hov1 = hovmoeller(x.num2date(x.time),xdata,rescaley=rescaley,lat=x.lat,rescalex=rescalex)
    hov2 = hovmoeller(y.num2date(y.time),ydata,rescaley=rescaley,lat=y.lat,rescalex=rescalex)

    hov1.time_to_lat(**kwargs); hov2.time_to_lat(**kwargs)

    cmap = plt.cm.get_cmap(data_cmap, nclasses)

    hov1.plot(title=x._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax1,showcolorbar=False,climits=climits,grid=grid)
    hov2.plot(title=y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax2,showcolorbar=False,climits=climits,grid=grid)

    add_nice_legend(ax1,hov1.im,cmap,cticks=cticks,label=clabel); add_nice_legend(ax2,hov2.im,cmap,cticks=cticks,label=clabel)

    if x.data.shape == y.data.shape:
        hov3 = hovmoeller(y.num2date(y.time),x.data - y.data,rescaley=rescaley,lat=y.lat,rescalex=rescalex)
        hov3.time_to_lat(**kwargs)
        cmap_diff = plt.cm.get_cmap('RdBu', nclasses)
        hov3.plot(title=x._get_label() + ' - ' + y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap_diff,ax=ax3,showcolorbar=False,climits=dlimits,grid=grid)
        add_nice_legend(ax3,hov3.im,cmap_diff,cticks=cticks_dif,label=clabel)
    else:
        msg = 'Difference plot not possible as data has different shape'
        ax3.text(0.5, 0.5,msg,
             horizontalalignment='center',
             verticalalignment='center') #,
             #transform = ax.transAxes)
        ax3.set_xticks([]); ax3.set_yticks([])

    return fig,hov1,hov2


#-----------------------------------------------------------------------


def map_difference(x, y, dmin=None, dmax=None, use_basemap=False,
                   ax=None, title=None, cticks=None,
                   region=None, nclasses=10, cmap_data='jet',
                   cmap_difference = 'RdBu_r', rmin=-1.,
                   rmax=1., absthres=None, show_stat=True,
                   show_zonal=True, zonal_timmean=False,
                   proj='robin', stat_type='mean', savefile=None,
                   savegraphics_prefix = None,
                   **kwargs):
    """
    Given two datasets, this routine generates a map plot of each dataset as
    well as of the difference of the two datasets

    Parameters
    ----------
    x : Data
        first dataset
    y : Data
        second dataset
    dmin : float
        minimum value of difference map
    dmax : float
        maximum value of difference map
    use_basemap : float
        flag if Basemap should be used for plotting
    ax : axis
        axis to plot to; if None [default], then new figure is generated
    title : str
        title of the plot
    nclasses : int
        number of classes for colormap
    cmap_data : str, colormap
        colormap for data to be plotted
    cmap_difference : str, colormap
        colormap for difference map to be plotted
    rmin : float
        minimum value for data plot
    rmax : float
        maximum value for data plot
    absthres : float
        threshold that will be estimated based on
        absolute difference and will be applied to
        relative difference maps
    show_stat : bool
        show statistics for each map
    show_zonal : bool
        show zonal means in each map
    zonal_timmean : bool
        use temporal mean for plotting zonal mean; applied only if
        show_zonal=True
    proj : str
        projection for map plotting
    stat_type : str
        type of statistic to be calculated (see map_plot documentation)
    savefile : str
        filename for netcdf file output
    savegraphics_prefix : str
        path and file prefix to save individual plots to graphic files.
        If given, the individual plots are saved to files. The
        filename needs to include already the file extension, as it is
        used to recognize the file extension.

        Example: savegraphics_prefix = '/my_path/outputfilename.png'
        will result the in output files:
           /my_path/outputfilename_X.png
           /my_path/outputfilename_Y.png
           /my_path/outputfilename_ADIFF.png
           /my_path/outputfilename_RDIFF.png
    """

    if savefile is not None:
        if '.nc' in savefile:
            savefile = savefile[:-3]

    if savegraphics_prefix is not None:
        graphic_rootname, extension = os.path.splitext(savegraphics_prefix)
    else:
        graphic_rootname = None
        extension = None

    if 'cticks_diff' in kwargs:
        cticks_diff = kwargs.pop('cticks_diff')
    else:
        cticks_diff = None

    if 'cticks_rdiff' in kwargs:
        cticks_rdiff = kwargs.pop('cticks_rdiff')
    else:
        cticks_rdiff = [-1.,-0.75,-0.5,-0.25,0.,0.25,0.5,0.75,1.]

    if 'colorbar_orientation' in kwargs:
        colorbar_orientation = kwargs.pop('colorbar_orientation')
    else:
        colorbar_orientation='horizontal'

    if 'drawparallels' in kwargs:
        drawparallels = kwargs.pop('drawparallels')
    else:
        drawparallels=False

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    #--- get colormap
    cmap = plt.cm.get_cmap(cmap_data, nclasses)

    #- temporal mean fields
    xm = x.timmean()
    ym = y.timmean()

    #proj='robin'; lon_0=0.; lat_0=0.

    #- plot first dataset
    if savefile is None:
        tmpoutname=None
    else:
        tmpoutname = savefile + '_xvar'

    if graphic_rootname is None:
        graphic_name = None
    else:
        graphic_name = graphic_rootname + '_X' + extension

    # plot dataset in entire figure
    map_plot(x, use_basemap=use_basemap,ax=ax1,cticks=cticks,region=region,nclasses=nclasses,
             cmap_data=cmap_data, title=title, show_stat=show_stat,show_zonal=show_zonal,
             zonal_timmean=zonal_timmean,proj=proj,stat_type=stat_type,savefile=tmpoutname,
             colorbar_orientation=colorbar_orientation,drawparallels=drawparallels, **kwargs)

    # ... do the same plot again but in an own figure!
    if graphic_rootname is not None:
        map_plot(x, use_basemap=use_basemap,cticks=cticks,region=region,nclasses=nclasses,
                 cmap_data=cmap_data, title=title, show_stat=show_stat,show_zonal=show_zonal,
                 zonal_timmean=zonal_timmean,proj=proj,stat_type=stat_type,savefile=tmpoutname,
                 colorbar_orientation=colorbar_orientation,drawparallels=drawparallels, savegraphicfile=graphic_name, **kwargs)

    #- plot second dataset
    if savefile is None:
        tmpoutname = None
    else:
        tmpoutname = savefile + '_yvar'

    if graphic_rootname is None:
        graphic_name = None
    else:
        graphic_name = graphic_rootname + '_Y' + extension

    # plot dataset in entire figure
    map_plot(y,use_basemap=use_basemap,ax=ax2,cticks=cticks,region=region,nclasses=nclasses,
             cmap_data=cmap_data, title=title,show_stat=show_stat,show_zonal=show_zonal,
             zonal_timmean=zonal_timmean,proj=proj,stat_type=stat_type,savefile=tmpoutname,
             colorbar_orientation=colorbar_orientation,drawparallels=drawparallels, **kwargs)

    # ... do the same plot again but in an own figure!
    if graphic_rootname is not None:
        map_plot(y,use_basemap=use_basemap,cticks=cticks,region=region,nclasses=nclasses,
                 cmap_data=cmap_data, title=title,show_stat=show_stat,show_zonal=show_zonal,
                 zonal_timmean=zonal_timmean,proj=proj,stat_type=stat_type,savefile=tmpoutname,
                 colorbar_orientation=colorbar_orientation,drawparallels=drawparallels, savegraphicfile=graphic_name, **kwargs)

    # first minus second dataset (absolute difference)
    adif = x.sub(y) #absolute difference #todo where to get std of seasonal means !!!! needs to be realized before beeing able to use significance ????

    if savefile is None:
        tmpoutname = None
    else:
        tmpoutname = savefile + '_absdif'

    if graphic_rootname is None:
        graphic_name = None
    else:
        graphic_name = graphic_rootname + '_ADIFF' + extension
    # entire plot
    map_plot(adif, use_basemap=use_basemap, ax=ax3, vmin=dmin, vmax=dmax,cticks=cticks_diff, region=region,
             nclasses=nclasses, cmap_data=cmap_difference, title='absolute difference',
             show_stat=show_stat, show_zonal=show_zonal, zonal_timmean=zonal_timmean, proj=proj, stat_type=stat_type, savefile=tmpoutname,
             colorbar_orientation=colorbar_orientation, drawparallels=drawparallels)

    # ... do the same plot again but in an own figure!
    if graphic_rootname is not None:
        map_plot(adif, use_basemap=use_basemap, vmin=dmin, vmax=dmax,cticks=cticks_diff, region=region,
                 nclasses=nclasses, cmap_data=cmap_difference, title='absolute difference',
                 show_stat=show_stat, show_zonal=show_zonal, zonal_timmean=zonal_timmean, proj=proj, stat_type=stat_type, savefile=tmpoutname,
                 colorbar_orientation=colorbar_orientation,
                 drawparallels=drawparallels, savegraphicfile=graphic_name)

    # relative error
    rdat = adif.div(x)
    if absthres is not None:
        mask = abs(x.timmean()) < absthres
        rdat._apply_mask(~mask)

    if savefile is None:
        tmpoutname = None
    else:
        tmpoutname = savefile + '_reldif'

    if graphic_rootname is None:
        graphic_name = None
    else:
        graphic_name = graphic_rootname + '_RDIFF' + extension

    map_plot(rdat, use_basemap=use_basemap, ax=ax4, vmin=rmin, vmax=rmax, title='relative difference',
             cticks=cticks_rdiff, region=region, nclasses=nclasses,
             cmap_data=cmap_difference,show_stat=show_stat, show_zonal=show_zonal,
             zonal_timmean=zonal_timmean, stat_type='median', proj=proj, savefile=tmpoutname,
             colorbar_orientation=colorbar_orientation, drawparallels=drawparallels)

    # ... do the same plot again but in an own figure!
    if graphic_rootname is not None:
        map_plot(rdat, use_basemap=use_basemap, vmin=rmin, vmax=rmax, title='relative difference',
                 cticks=cticks_rdiff, region=region, nclasses=nclasses,
                 cmap_data=cmap_difference,show_stat=show_stat, show_zonal=show_zonal,
                 zonal_timmean=zonal_timmean, stat_type='median', proj=proj, savefile=tmpoutname,
                 colorbar_orientation=colorbar_orientation, drawparallels=drawparallels, savegraphicfile=graphic_name)


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



    h = hovmoeller(x.num2date(x.time),x.data,rescaley=rescaley,lat=x.lat,rescalex=rescalex)
    h.time_to_lat(dlat=dlat,monthly = True, yearonly = True,monthsamp=monthsamp)
    h.plot(title=tit,xlabel=xlab,ylabel=ylabel,origin='lower',xtickrotation=xtickrotation,cmap=cmap,ax=ax,showcolorbar=False,climits=climits,grid=False,showxticks=showxticks)

    add_nice_legend(ax,h.im,cmap,cticks=cticks)

    return h



def calc_rms_error(x,y):
    """
    calculate RMS error

    x,y: masked arrays
    """

    xdat = x.flatten(); ydat=y.flatten()

    return np.sqrt(np.ma.mean((xdat-ydat)**2.))

def calc_centered_rms_error(x,y):
    """
    calculate centererd RMS error

    REFERENCES:
     * Taylor et al. (2001), eq. 2
    """
    xdat = x.flatten(); ydat=y.flatten()
    xm = np.ma.mean(xdat); ym = np.ma.mean(ydat)

    anom = np.sqrt(np.ma.mean(((xdat-xm) - (ydat-ym))**2.))

    return xm-ym,anom



