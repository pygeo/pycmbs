#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

'''
Module that contains relevant classes for diagnostic plots

@todo: implement writing of statistics to an ASCII file as export
@todo: implement taylor plots
@todo: faster implementation of Basemap plots. For large number of grid cells, the current KTree implementation is by far too slow!

'''

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
    '''
    thin xticks of axis

    If there are too manx xticks in a plot or the labels
    are overlapping, it makes sense to thin the mńumber of labels

    @param ax: axis that will be treated
    @type ax: matplotlib axis

    @param n: number of ticks to plot
    @type n: int
    '''
    ax.xaxis.set_major_locator(MaxNLocator(n+1))

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class CorrelationAnalysis():
        '''
        perform correlation analysis between two datasets
        and plot results
        '''

        def __init__(self,X,Y,mask=None,ax=None):
            '''
            constructor of class

            @param X: x dataset (either [time,sample] or [time,sample,sample]
            @type X: numpy array

            @param Y: y dataset (either [time,sample] or [time,sample,sample]
            @type Y: numpy array

            @param mask: mask to be applied to the data
            @type mask: numpy array(:,:) or (:)

            @param ax: axis to plot results to; new figure will be generated if ax==None
            @type ax: matplotlib axis
            '''

            self.x = X; self.y = Y
            self.mask = mask

            if ax == None:
                f = plt.figure()
                self.ax = f.add_subplot(111)
            else:
                self.ax = ax

#-----------------------------------------------------------------------

        def do_analysis(self):
            '''
            perform correlation analysis

            @todo: implement area weighting
            @todo: implement regional (condition) statisitcs based on a mask
            @todo: return value
            '''

            #--- calculate diagnostics
            D = Diagnostic(self.x,y=self.y)
            D._mat2vec(mask = self.mask) #here is the point fo rregional statistics
            rmse = D.get_rmse_value()
            r    = D.get_correlation_value()
            n    = D. get_n()

            print 'RMSE: ', rmse
            print 'R   : ', r
            print 'N   : ', n

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
class HovmoellerPlot():
    def __init__(self,D,rescaley=10,rescalex=10,dlat=1,yticksampling=1,monthly=False):
        '''

        D : C{Data} object

        if the argument lat is provided it is assumed that lat/lon are 2D matrices
        In this case the value is expected to be a 3D variables as
        value(time,ny,nx)
        '''
        #~ from python.hov import *
        self.hov = hovmoeller(pl.num2date(D.time),D.data,lat=D.lat,rescaley=rescaley,rescalex=rescalex)
        self.hov.time_to_lat(dlat=dlat,yticksampling=yticksampling,monthly=monthly)

    def plot(self,title=None,climits=None):
        if climits == None:
            raise ValueError, 'CLIMITS needs to be specified!'
        self.hov.plot(title=title,ylabel='lat',xlabel='days',origin='lower',xtickrotation=30,climits=climits)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class ReichlerPlot():
    '''
    class for Reichler plot generation

    @todo: add example how to use Reichler plotting
    @todo: provide references Glecher and Reichler + Kim
    '''
    def __init__(self,ax=None):
        '''
        constructor for Reichler plot

        @param ax: axis to plot data to; if None, new figure will be generated
        @type ax: matplotlib axis
        '''
        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

        self.e2 = [] #list to store RMS error results
        self.labels = []; self.colors=[]

#-----------------------------------------------------------------------

    def add(self,e2,label,color=None):
        '''
        register data to be plotted

        @param e2: reichler index that was already calculated
        @type e2: list

        @param label: label to be used for plotting
        @type label: str

        @param color: color to be used for plotting
        @type color: str
        '''
        self.e2.append(e2)
        self.labels.append(label)
        self.colors.append(color)

#-----------------------------------------------------------------------

    def bar(self,vmin=None,vmax=None,title='',**kwargs):
        '''
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
        '''
        print 'Doing Reichler plot as barplot ...'
        self._normalize()
        x = np.arange(len(self.e2_norm))
        self.ax.bar(x,self.e2_norm*100.,**kwargs)
        self.ax.set_xticks(x+0.5)
        self.ax.set_xticklabels(self.labels)
        if (vmin !=None) & (vmax != None):
            self.ax.set_ylim(vmin,vmax)
        self.ax.set_ylabel('$\\epsilon / \\bar{\\epsilon}$ [%]')
        self.ax.grid(); self.ax.set_title(title)


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

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class ScatterPlot():
    '''
    Class for generation of scatterplots
    '''
    def __init__(self,x,ax=None,ticksize=10):
        '''
        constructor of class C{ScatterPlot}

        @param x: Variable that will be used as the x-variable
        @type x: C{Data} object
        '''

        if ax == None:
            f = plt.figure()
            self.ax = f.add_subplot(111)
        else:
            self.ax = ax

        self.figure = self.ax.figure
        self.x = x
        self.lines = []; self.labels = []
        self.ticksize=ticksize

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

        #- calculate linear regression
        if regress:
            slope, intercept, r_value, p_value, std_err = stats.linregress(xdat,ydat)
            label = label + ' (r=' + str(round(r_value,2)) + ', p=' + str(round(p_value,2)) + ')'

        l = self.ax.plot(xdat,ydat,'.',label=label,**kwargs)[0]
        if regress:
            self.ax.plot(xdat,xdat*slope+intercept,'--',color=l.get_color())
        self.lines.append(l); self.labels.append(label)

        self.ax.set_xlabel(self.x._get_label(),size=self.ticksize )
        self.ax.set_ylabel(y._get_unit(),size=self.ticksize)

        self._change_ticklabels()

    def _change_ticklabels(self):
        for tick in self.ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(self.ticksize)
        for tick in self.ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)

#-----------------------------------------------------------------------

    def legend(self):
        '''
        plot legend
        '''
        self.ax.legend(self.lines,self.labels,prop={'size':8})

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class LinePlot():
    '''
    class for a pyCMBS Line Plot

    This class is usefull for plotting timeseries
    '''
    def __init__(self,ax=None,regress=False,title=None,show_xlabel=True,show_ylabel=True,ticksize=10,normx=1.):
        '''
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
        '''

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

        self.lines = []; self.labels = []

        self.ticksize = ticksize

        self.normx=normx

#-----------------------------------------------------------------------

    def legend(self):
        '''
        plot legend
        '''
        self.ax.legend(self.lines,self.labels,prop={'size':10})

#-----------------------------------------------------------------------

    def _change_ticklabels(self,ax=None):
        if ax == None:
            ax = self.ax

        #~ for tick in ax.xaxis.get_major_ticks():
            #~ tick.label.set_fontsize(self.ticksize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.ticksize)

#-----------------------------------------------------------------------

    def plot(self,x,ax=None,vmin=None,vmax=None,label = None, **kwargs):
        '''
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
                slope_print, intercept_print, r_value, p_value, std_err = stats.linregress(x.time/self.normx,y) #@todo: is it correct to use here time instead of .data?
                slope, intercept, r_value, p_value, std_err = stats.linregress(x.time,y) #@todo: is it correct to use here time instead of .data?
                label = label + ' (y=' + "%.1e" % slope_print + 'x+' + "%.1e" % intercept_print  + ', r=' + str(round(r_value,2)) + ', p=' + str(round(p_value,2)) + ')'

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
                ax.set_title(self.title)

            if vmin != None:
                if vmax != None:
                    ax.set_ylim(vmin,vmax)


            self._change_ticklabels(ax)

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

    def plot(self,x,areaweights,xlim=None,timmean = False,show_ylabel=True):
        '''
        plot zonal plot

        @param x: data to be plotted
        @type x: C{Data} object

        @param areaweights: cell weights for the area
        @type areaweights: numpy array

        @param xlim: limits for the x-axis (e.g. values)
        @type xlim: tuple

        @param timmean: temporal mean calculation
        @type timmean: bool

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

        if timmean:
            thex = x.timmean(return_object=True)
        else:
            thex = x

        if self.dir == 'y':
            dat = thex.get_zonal_statistics(areaweights) #no area weighting performed
        else:
            raise ValueError, 'Invalid option'

        if timmean:
            pass #todo
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
                #~ print 'Time in zonal: ', i
                #~ print dat[i,:]
                #~ self.ax.plot(dat[i,:],label='time='+str(i))
                self.ax.plot(dat[i,:],x.lat[:,0],label='time='+str(i))



        self.ax.set_ylim(-90.,90.)

        if show_ylabel:
            self.ax.set_ylabel('latitude [deg]')
        else:
            self.ax.set_yticks([])

        if xlim != None:
            self.ax.set_xlim(xlim)

        self.ax.grid()

#-----------------------------------------------------------------------

class GleckerPlot():
    '''
    Class to generate a plot that to illustrate multi-model, multi-variable scores

    It was introdcued by Glecker et al (2008)

    REFERENCES:
    * ﻿Gleckler, P.J., Taylor, K.E. & Doutriaux, C., 2008. Performance metrics for climate models. Journal of Geophysical Research, 113(D6). Available at: http://www.agu.org/pubs/crossref/2008/2007JD008972.shtml [Accessed February 29, 2012].

    EXAMPLE:
    G = GleckerPlot()
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
    '''

    def __init__(self,fig=None):
        '''
        constructor of C{GleckerPlot}

        @param fig: figure to which to plot to. If None, then a new figure will be generated
        @type fig: matplotlib figure
        '''
        if fig == None:
            color='grey'
            fig = plt.figure(facecolor=color,edgecolor=color)
        self.fig = fig

        self.models = []
        self.variables   = []
        self.data = {} #store data for plot
        self.pos = {} #store position of plot

    def add_model(self,label):
        '''
        register a model in the class
        @param label: string of the model
        @type label: str
        '''
        s = label.replace(' ','_')
        if s not in self.models:
            self.models.append(s)

    def add_variable(self,label):
        '''
        register a variable in the class
        @param label: string of variable
        @type label: str
        '''
        self.variables.append(label)

    def __set_ax_prop(self,ax):
        '''
        set axis properties of a subplot
        @param ax: subplot axis
        @type ax: matplotlib axis
        '''
        ax.set_xticks([]); ax.set_yticks([])

    def __value2color(self,v):
        '''
        return a color based on the value given
        the information on the colormap and its
        normalization is used for that purpose

        @param v: value of data
        @type v: float
        '''
        return self.cmap(self.norm(v))

    def __plot_triangle(self,ax,value,pos='top'):
        '''
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
        '''
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
        '''
        calculate for each observational data set
        the relative deviation from the average
        '''
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
        '''
        calculate mean value for a given observational dataset

        @param pos: position marker
        @type pos: int

        @param var: name of variable to analyze
        @type var: str
        '''
        x = []
        for k in self.pos:
            if (self.pos[k] == pos) & ('_' + var + '_' in k):
                x.append(self.data[k])
        x = np.asarray(x)
        return x.mean()



#-----------------------------------------------------------------------

    def plot(self,cmap_name='RdBu_r',vmin=-1.,vmax=1.,nclasses=15,normalize=True,size=10):
        '''
        plot Glecker diagram

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
        '''

        #~ todo implement normalization here per position

        if normalize:
            self._normalize_data()

        nm = len(self.models); nv = len(self.variables)

        #- colormap
        self.cmap = plt.cm.get_cmap(cmap_name, nclasses)
        self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        gs = gridspec.GridSpec(nm, nv, wspace=0.05,hspace=0.05,bottom=0.2) #generate grid for subplots

        cnt = 0; cnt_m = 0

        model_list = self.models.sort()

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


                cnt += 1
                cnt_v += 1
        #--- legend

        #- get positions of subplots to determine optimal position for legend
        def get_subplot_boundaries(g,f):
            x = g.get_grid_positions(f)
            b = x[0]; t = x[1]
            l = x[2]; r = x[3]
            return l[0], r[-1], b[-1], t[0]

        left,right,bottom,top = get_subplot_boundaries(gs,self.fig)
        #draw legend
        self._draw_legend(left,right-left)

#-----------------------------------------------------------------------

    def get_data(self,v,m,p):
        '''
        return data for a particular model and variable

        @param v: name of variable
        @type v: str

        @param m: model name
        @type m: str

        @param p: position
        @type p: int
        '''
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
        '''
        generate a unique key for dictionaries
        comprised of model name, variable name and position

        @param v: name of variable
        @type v: str

        @param m: model name
        @type m: str

        @param p: position
        @type p: int
        '''
        if m == None:
            return None
        if v == None:
            return None
        return m.replace(' ','_')+'_'+v.replace(' ','_')+'_'+str(p)

#-----------------------------------------------------------------------

    def add_data(self,v,m,x,pos=1):
        '''
        add a data for plotting

        @param v: name of variable
        @type v: str

        @param m: model name
        @type m: str

        @param x: value
        @type x: float

        @param pos: position where to plot data 1=top triangle, 2=lower triangle
        @type pos: int
        '''

        if v in self.variables:
            if m in self.models:
                self.data.update({ self.__gen_key(m,v,pos) :x})
                self.pos.update({ self.__gen_key(m,v,pos) : pos})

#-----------------------------------------------------------------------

    def calc_index(self,x,y,model,variable,weights=None):
        '''
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

        @return: returns performance index aggregated over time
        @rtype float
        '''

        if weights == None:
            print 'WARNING: no weights when calculating performance index'
            weights = np.ones(x.data[0,:].shape)

        from diagnostic import Diagnostic
        D = Diagnostic(x,y=y)
        e2 = D.calc_reichler_index(weights) #reichler performance index (might return a list if multiple times analyzed)

        r = np.nansum(e2) #temporal aggregation

        return r

#-----------------------------------------------------------------------

    def _draw_legend(self,left,width):
        '''
        draw legend for Glecker plot. Requires information on
        the positioning of the colormap axis which can be obtained from

        left,right,bottom,top = get_subplot_boundaries(gs,self.fig)

        @param left: left position of axis
        @type left: float

        @param width: width of the axis to plot colorbar
        @type width: float
        '''
        cax = self.fig.add_axes([left,0.05,width,0.05]) #left, bottom, width, height
        cb = mpl.colorbar.ColorbarBase(cax, cmap=self.cmap,
                                   norm=self.norm,
                                   orientation='horizontal')

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def __basemap_ancillary(m,latvalues = None, lonvalues = None):
    '''
    routine to plot ancillary data like coastlines
    or meridians on a basemap plot

    @param m: map to add features to
    @type m: C{Basemap} object

    @param latvalues: latitude values for drawing grid (optional)
    @type latvalues: list or numpy array

    @param lonvalues: longitude values for drawing grid (optional)
    @type lonvalues: list or numpy array

    '''

    if latvalues == None:
        latvalues=np.arange(-90.,120.,30.)
    if lonvalues == None:
        lonvalues= np.arange(-180.,180.,90.)
    m.drawcountries(); m.drawcoastlines()
    m.drawlsmask(lakes=True)
    m.drawmapboundary() # draw a line around the map region
    m.drawparallels(latvalues,labels=[1, 0, 0, 0])
    m.drawmeridians(lonvalues,labels=[0, 0, 0, 1]) # draw meridians

#-----------------------------------------------------------------------

def map_season(x,**kwargs):
    '''
    generate a seasonal plot
    all arguments are parsed directly to map_plot function

    if kwargs contain a 'figure' argument, then this figure fill be used
    for plotting. Otherwise a new figure will be generated

    @param x: C{Data} object
    @type x : C{Data}

    @return: returns the figure where plot was done
    @rtype: C{figure}

    '''

    #/// checks ///
    if x.data.ndim != 3:
        raise ValueError, 'only 3D data supported'

    if len(x.data) != 4:
        raise ValueError, 'only data with four seasons supported'

    #/// figure and axes
    if 'figure' in kwargs:
        f = kwargs['figure']
    else:
        f = pl.figure()

    #/// plot
    labels=['JFM','AMJ','JAS','OND']

    for i in range(4):
        ax = f.add_subplot(2,2,i+1)
        d = x.copy(); d.data = x.data[i,:,:]
        d.label = labels[i]
        map_plot(d,ax=ax,**kwargs); del d

    f.suptitle(x.label,size=12)

    return f

#-----------------------------------------------------------------------

def map_plot(x,use_basemap=False,ax=None,cticks=None,region=None,nclasses=10,cmap_data='jet', title=None,regions_to_plot = None,logplot=False,logoffset=None,show_stat=False, f_kdtree=False,show_colorbar=True,latvalues=None,lonvalues=None,show_zonal=False,zonal_timmean=True, **kwargs):
    '''
    produce a nice looking map plot

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

    @param show_stat: show statistic of field in figure title
    @type show_stat: bool

    @param f_kdtree: use kdTree for interpolation of data to grid (might be slow, but might solve problem of stripes in plots)
    @type f_kdtree: bool

    @param latvalues: latitude values for drawing grid (optional)
    @type latvalues: list or numpy array

    @param lonvalues: longitude values for drawing grid (optional)
    @type lonvalues: list or numpy array

    '''

    #--- create new figure
    if ax == None:
        fig = plt.figure(); ax = fig.add_subplot(111)
    else:
        fig = ax.figure

    #if cmap provided in kwargs, then remove it and set cmap_data
    kwargs1 = kwargs.copy()
    if 'cmap' in kwargs:
        cmap_data = kwargs1.pop('cmap')

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
                di = 0. #with 0 it works; for other values problems may occur for negative lon!
                llcrnrlon=region.lonmin - di; llcrnrlat=region.latmin - di
                urcrnrlon=region.lonmax + di; urcrnrlat=region.latmax + di
                proj='tmerc' #use mercator projection at regional scale as robinson does not work!
                #~ proj = 'cyl'

        #generate map
        m1=Basemap(projection=proj,lon_0=lon_0,lat_0=lat_0,ax=ax,llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)


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

        else: #f_kdtree
            X,Y = m1(x.lon,x.lat)
            Z = xm

        #here is still a problem in the plotting over land; masking does not work properly,
        #while the data as such is o.k.!
        #~ im1=m1.pcolormesh(xmap,ymap,xm,cmap=cmap,**kwargs) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        im1=m1.pcolormesh(X,Y,Z,cmap=cmap,**kwargs1) #,vmin=vmin,vmax=vmax,cmap=ccmap,norm=norm)
        __basemap_ancillary(m1,latvalues=latvalues,lonvalues=lonvalues)

    else: #use_basemap = False
        #- normal plots
        im1=ax.imshow(xm,cmap=cmap,interpolation='nearest', **kwargs1)
        ax.set_xticks([]); ax.set_yticks([])




    #set legend aligned with plot (nice looking)
    if show_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal("5%", pad=0.05, axes_class=maxes.Axes)
        ax.figure.add_axes(cax)
        norm = mpl.colors.Normalize(vmin=im1.get_clim()[0], vmax=im1.get_clim()[1])
        cb   = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,ticks=cticks)

    #Zonal plot
    if show_zonal:
        add_zonal_plot(ax,x,timmean=zonal_timmean) #,vmin=im1.get_clim()[0],vmax=im1.get_clim()[1])





    def _add_region(m,r,color='red'):
        '''
        plot region r on top of basemap map m

        @param m: map
        @type m: C{Basemap} object

        @param r: region to plot
        @type r: C{Region}

        @param color: color to plot region
        @type color: str
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
        title = x._get_label()
    else:
        pass

    #--- show field statistics in title ?
    if show_stat:
        me = xm.mean(); st=xm.std()
        title = title + '\n ($' + str(round(me,2))  + ' \pm ' + str(round(st,2)) + '$' + x._get_unit() + ')'


    ax.set_title(title,size=12)


    return fig

#-----------------------------------------------------------------------

def add_zonal_plot(ax,x,timmean=True,vmin=None,vmax=None):
    '''
    add a zonal plot to the axis

    @param ax: axis where zonal plot should be added to
    @type ax: axis

    @param x: data to plot
    @type x: C{Data} object

    @param timmean: temporal mean for zonal plot?
    @type timmean: bool
    '''

    divider = make_axes_locatable(ax)
    zax = divider.new_horizontal("15%", pad=0.05, axes_class=maxes.Axes,pack_start=True)
    ax.figure.add_axes(zax,axisbg=ax.figure.get_facecolor())

    ZP = ZonalPlot(ax=zax,dir='y')

    if x.data.ndim == 2:
        weights = np.ones(x.data.shape)
    elif x.data.ndim == 3:
        nt,ny,nx = x.data.shape
        weights = np.ones((ny,nx))
    weights = np.ma.array(weights,mask = weights != weights)
    ZP.plot(x,weights,timmean=timmean,show_ylabel=False) #@todo: area weighting??

    #set only first and last label
    #~ lab = zax.get_xticklabels()
    ti  = zax.get_xticks()
    zax.set_xticks([ti[0],ti[-1]])

    return zax

    #~ zax.set_xlim(vmin,vmax)

#-----------------------------------------------------------------------

def add_nice_legend(ax,im,cmap,cticks=None,dummy=False):
    '''
    add a nice looking legend

    @param ax: major axis with plot
    @type ax: matpltlib axis object

    @param im: result from command like imshow
    @param im: matplotlib im object (???)

    @param dummy: add colorbar axis as a dummy axis which is not visible
                  this is useful if you have multiple subplots which should
                  have the same size. Adding a colorbar will slightly change the size
    @type dummy: bool
    '''

    #set legend aligned with plot (nice looking)
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal("5%", pad=0.05, axes_class=maxes.Axes)
    ax.figure.add_axes(cax,axisbg=ax.figure.get_facecolor())
    if dummy:
        cax.set_xticks([])
        cax.set_yticks([])
        cax.set_frame_on(False)
    else:
        norm = mpl.colors.Normalize(vmin=im.get_clim()[0], vmax=im.get_clim()[1])
        cb   = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm,ticks=cticks)

#-----------------------------------------------------------------------

def hov_difference(x,y,climits=None,dlimits=None,data_cmap='jet',nclasses=15,cticks=None,cticks_dif=None,ax1=None,ax2=None,ax3=None,rescaley=6,grid=True,rescalex=1,**kwargs):
    '''

    class to plot hovmoeller diagrams of two datasets
    and their difference

    x,y two Data structures

    axextra: plot difference on separate axis
    '''

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
    xdata = x.data
    ydata = y.data

    #~ xdata = x.data.data
    #~ xdata[x.data.mask] = np.nan
#~
    #~ ydata = y.data.data
    #~ ydata[y.data.mask] = np.nan


    hov1 = hovmoeller(num2date(x.time),xdata,rescaley=rescaley,lat=x.lat,rescalex=rescalex)
    hov2 = hovmoeller(num2date(y.time),ydata,rescaley=rescaley,lat=y.lat,rescalex=rescalex)


    hov1.time_to_lat(**kwargs)
    hov2.time_to_lat(**kwargs)


    cmap = plt.cm.get_cmap(data_cmap, nclasses)


    hov1.plot(title=x._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax1,showcolorbar=False,climits=climits,grid=grid)
    hov2.plot(title=y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax2,showcolorbar=False,climits=climits,grid=grid)

    add_nice_legend(ax1,hov1.im,cmap,cticks=cticks)
    add_nice_legend(ax2,hov2.im,cmap,cticks=cticks)

    #~ plt.colorbar(hov1.im,ax=ax1,shrink = 0.5,orientation='vertical')
    #~ plt.colorbar(hov2.im,ax=ax2,shrink = 0.5,orientation='vertical')


    if x.data.shape == y.data.shape:
        hov3 = hovmoeller(num2date(y.time),x.data - y.data,rescaley=rescaley,lat=y.lat,rescalex=rescalex)
        hov3.time_to_lat(**kwargs)
        cmap_diff = plt.cm.get_cmap('RdBu', nclasses)
        hov3.plot(title=x._get_label() + ' - ' + y._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap_diff,ax=ax3,showcolorbar=False,climits=dlimits,grid=grid)
        #~ plt.colorbar(hov3.im,ax=ax3,shrink = 0.5,orientation='vertical')
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
    '''
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
    map_plot(x,use_basemap=use_basemap,ax=ax1,cticks=cticks,region=region,nclasses=nclasses,cmap_data=cmap_data, title=title,show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean, **kwargs)

    #- plot second dataset
    map_plot(y,use_basemap=use_basemap,ax=ax2,cticks=cticks,region=region,nclasses=nclasses,cmap_data=cmap_data, title=title,show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean,  **kwargs)

    #-first minus second dataset
    adif = x.sub(y) #absolute difference
    map_plot(adif,use_basemap=use_basemap,ax=ax3,vmin=dmin,vmax=dmax,cticks=None,region=region,nclasses=nclasses,cmap_data=cmap_difference, title='absolute difference [' + x.unit + ']',show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean)

    #~ map_plot(adif,use_basemap=use_basemap,vmin=dmin,vmax=dmax,cticks=None,region=region,nclasses=nclasses,cmap_data=cmap_difference, title='absolute difference [' + x.unit + ']')


    #- relative error
    rdat = adif.div(x) #y.div(x).subc(1.) #relative data
    if absthres != None:
        mask = abs(x.timmean()) < absthres
        rdat._apply_mask(~mask)

    map_plot(rdat,use_basemap=use_basemap,ax=ax4,vmin=rmin,vmax=rmax,title='relative difference',cticks=None,region=region ,nclasses=nclasses,cmap_data=cmap_difference,show_stat=show_stat,show_zonal=show_zonal,zonal_timmean=zonal_timmean)


    return fig

#-----------------------------------------------------------------------

def plot_hovmoeller(x,rescaley=10,rescalex=1,monthsamp=24,dlat=1.,cmap=None,ax=None,climits=None):
    '''
    plot hovmoeller plots given a C{Data} object

    @param x: C{Data} object
    @type x: C{Data} object
    '''

    if climits == None:
        raise ValueError, 'climits need to be given!'

    if cmap == None:
        cmap = plt.cm.get_cmap('RdBu', 15)
    if ax == None:
        f = plt.figure()
        ax = f.add_subplot(111)

    h = hovmoeller(pl.num2date(x.time),x.data,rescaley=rescaley,lat=x.lat,rescalex=rescalex)
    h.time_to_lat(dlat=dlat,monthly = True, yearonly = True,monthsamp=monthsamp)
    h.plot(title=x._get_label(),ylabel='lat',xlabel='time',origin='lower',xtickrotation=30,cmap=cmap,ax=ax,showcolorbar=False,climits=climits,grid=False)

    add_nice_legend(ax,h.im,cmap,cticks=None)

    return h



