#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1.4"
__date__ = "2012/10/29"
__email__ = "alexander.loew@mpimet.mpg.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@mpimet.mpg.de
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

from matplotlib import pylab as plt

import numpy as np


class Taylor(object):
    def __init__(self, stdmax=2., plot_r_mesh=True, plot_std_mesh=True,
                 ref_std=1., plot_reference=True, r_meshstep=0.1, std_meshstep=0.25,
                 edgecolor='black', r_meshcolor='red', std_meshcolor='blue', title='',
                 r_equidistant=False, normalize=False, figsize=(8, 6), rms_meshstep=0.1, maxrms = 20., show_negative_r=True):
        """
        initialize Taylor diagramm class

        Options
        * edgecolor: specifies the color of the edges
        * stdmax   : specifies the limits of the x/y axis (needs to be defined from data in advance)
        * ref_std  : reference standard deviation
        * plot_r_mesh : plot meshgrid for correlation
        * plot_std_mesh : plot meshgrid for standard deviation
        * r_meshcolor: color of correlation mesh
        * std_meshcolor : color of standard deviation mesh
        * plot_reference : plot the reference circle
        * r_meshstep : stepsize for correlation mesh
        * std_meshstep : stepsize for standard deviation mesh
        * r_equidistant: plot correlations equidistant, if FALSE, then standard approach of Taylor is used with NOT equidistant plots
        * normalize : normalize standard deviation by value provided in ref_std


        EXAMPLE:
        #generate some sample data
        corr1 = asarray([0.,.5,1.]); corr2=asarray([0.1,0.25,0.5])
        s1 = asarray([0.1,1.,0.6]); s2 = asarray([0.9,1.5,1.7])

        #Taylor plot
        tay=taylor() #initialize the Taylor diagram
        tay.plot(corr1,s1,markerfacecolor='green',marker='^',label='test1') #plot some data
        tay.plot(corr2,s2,markerfacecolor='magenta',marker='*',label='test2') #plot some other data
        tay.ax.legend(loc='lower center',ncol=2,fancybox=True,shadow=True) #if you like, access the axis object and specify legend properties

        #plotting with legeneds for groups
        tay.plot(res_station['LE_daily']['corr'],res_station['LE_daily']['rms'],marker='d',markerfacecolor=color) #generate some plot
        tay.plot(res_station['LE_daily']['corr'],res_station['LE_daily']['rms'],marker='d',markerfacecolor=color) #and another plot for the same group
        tay.set_legend_plot('station',color=color,marker='o') #define the legend
        ...
        tay.ax.legend(tay.plots,tay.labelnames,loc=[0.8,0.8],ncol=1,fancybox=True,shadow=True, scatterpoints=1)  #plot the legend



        TODO
        * implement RMSE plot
        * implement normalization option
        * negative correlations




        """

        self.normalize = normalize
        self.figure = plt.figure(facecolor='w', edgecolor='w', figsize=figsize)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_aspect('equal')
        self.ax.set_title(title)
        self.ax.set_frame_on(False)
        self.stdmax = stdmax
        self.ax.plot([0, stdmax], [0, 0.00001], color=edgecolor)
        self.ax.plot([0, 0], [0, stdmax], color=edgecolor)
        self.ref_std = ref_std

        self.r_equidistant = r_equidistant  # plot correlations equidistant?

        self.rms_meshstep = rms_meshstep
        self.maxrms = maxrms

        #set ticks only for bottom and left
        self.ax.xaxis.set_ticks_position('none')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('none')
        self.ax.yaxis.set_ticks_position('left')

        #plot list for legend
        self.plots = []
        self.labelnames = []

        self.r_color = r_meshcolor
        self.rms_color = 'grey'
        self.std_color = std_meshcolor

        self.negative = show_negative_r

        #/// meshlines ///
        if plot_r_mesh:
            self.plot_r_meshlines(step=r_meshstep)
        if plot_std_mesh:
            self.plot_std_meshlines(step=std_meshstep)

        #/// plot reference circle and point ///
        if plot_reference:
            self.plot_reference_circle()
            self.plot_reference_point()

    def plot_r_meshlines(self, step=0.1):
        """
        plots mesh lines for correlation
        """

        color = self.r_color

        nstdmax = self.stdmax
        if self.negative:
            axmin = -1.
            addlist = [-0.99, -0.95, 0.95, 0.99]
        else:
            axmin = 0.
            addlist = [0.95, 0.99]
        corrs = np.arange(axmin, 1.+step, step)
        if not self.r_equidistant:
            corrs = np.asarray(list(corrs) + addlist)

        for corr1 in corrs:  # arange(0.,1.+step,step):
            c_theta = np.deg2rad((1.-corr1)*90.)
            c_r = nstdmax
            if self.r_equidistant:
                self.ax.plot([0, c_r*np.cos(c_theta)], [0, c_r*np.sin(c_theta)], ':', color=color)
                xa, ya = self.map2xy(corr1, 1.0*nstdmax)
                #self.ax.annotate(str(corr1), [c_theta, nstdmax*1.11], xycoords='polar',color=color)
                self.ax.annotate(str(corr1), [xa, ya], color=color, horizontalalignment='center')

            else:
                self.ax.plot([0, c_r*corr1], [0, c_r*np.sin(np.arccos(corr1))], ':', color=color)
                xa, ya = self.map2xy(corr1, 1.0*nstdmax)
                #self.ax.annotate(str(corr1), [np.arccos(corr1), nstdmax*1.11], xycoords='polar',color=color)
                self.ax.annotate(str(corr1), [xa, ya], color=color, horizontalalignment='center')

        if self.r_equidistant:
            pass
        else:
            self.ax.annotate('Correlation coefficient $R$', [np.deg2rad(50.), nstdmax*1.1], xycoords='polar', rotation=-45., color=color)

    def plot_circle(self, x, y, r, color='grey', label=None, size=8):
        """
        plot a circle at point x,y with radius r
        """

        X = np.linspace(-r, +r, 1000.)
        Y = np.sqrt(r*r - X*X)
        XN = X+x
        YN = Y+y
        m = XN < 0.
        XN[m] = np.nan
        YN[m] = np.nan
        m = (XN**2 + YN**2) >= self.stdmax**2
        XN[m] = np.nan
        YN[m] = np.nan

        #calculate label position
        alpha = np.arctan(self.stdmax / x)
        xlab = x - r * np.cos(alpha)
        ylab = y + r * np.sin(alpha)

        self.ax.plot(XN, YN, '--', color=color)
        if label is not None:
            if ylab < 0.9*self.stdmax:  # avoid labels at boundaries
                self.ax.annotate(label, [xlab, ylab], xycoords='data', rotation=30., color=color, backgroundcolor='white', verticalalignment='center', horizontalalignment='center', size=size)

    def plot_rms_meshlines(self):
        """
        plot RMS meshlines

        todo labesl still not work properly !
        """

        color = self.rms_color

        rms = np.arange(0., self.maxrms, self.rms_meshstep)  # rms values to plot

        n = len(rms)
        for i in range(n):
            r = rms[i]
            self.plot_circle(self.ref_std, 0., r, color=color, label=str(r))

    def plot_std_meshlines(self, step=0.1):
        '''
        plot mesh circles for stdv
        '''

        color = self.std_color

        nstdmax = self.stdmax
        if self.negative:
            axmin = -np.pi/2.
        else:
            axmin = 0.

        th = np.arange(axmin, np.pi/2, 0.01)

        for ra in np.arange(0, nstdmax+0.1*step, step):
            self.ax.plot(ra*np.sin(th), ra*np.cos(th), ':', color=color)

        if self.normalize:
            self.ax.set_ylabel('$\sigma / \sigma_{obs}$', color=color)
            self.ax.set_xlabel('$\sigma / \sigma_{obs}$', color=color)
        else:
            self.ax.set_ylabel('Standard Deviation', color=color)
            self.ax.set_xlabel('Standard Deviation', color=color)

        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, color=color)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, color=color)

    def plot_reference_circle(self, linestyle='-', linewidth=2, color='orange', marker='None'):
        """
        plot reference circle of equal std which corresponds to ref_std
        """
        R = np.linspace(-1., 1., 1000)
        S = np.ones(len(R))*self.ref_std
        self.plot(R, S, linestyle=linestyle, linewidth=linewidth, color=color, marker=marker)

    def plot_reference_point(self, color='orange', marker='o'):
        """
        plot reference point at location with r=1 and std=ref_std
        """
        R = np.asarray([1])
        S = np.asarray([self.ref_std])
        self.plot(R, S, color=color, marker=marker)

    def plot_taylor_legend(self):
        """
        plots a legend for the relationship between
        R,RMSE and std
        """

        ax = self.figure.add_axes([0.8, 0.8, 0.15, 0.15], frameon=False)
        ax.set_xticks([])
        ax.set_yticks([])

        lw = 2
        P1 = [0., 0.]
        P2 = [1., 1.]
        P3 = [0.7, 0.]
        ax.plot([P1[0], P2[0]], [P1[1], P2[1]], '-', color=self.r_color, linewidth=lw)
        ax.plot([P1[0], P3[0]], [P1[1], P3[1]], '-', color=self.std_color, linewidth=lw)
        ax.plot([P3[0], P2[0]], [P3[1], P2[1]], '-', color=self.rms_color, linewidth=lw)

    def plot(self, R, S, marker='o', linestyle='None', markerfacecolor='black', markeredgecolor=None,
             color='black', linewidth=0, label=None, markersize=8, markeredgewidth=1, R1=None, S1=None,
             labels=None, shiftcolor='grey', plot_mean=False, normvalue=np.nan):
        """
        R: correlation coefficient
        S: standard deviation

        normvalue: value to use for normalization (e.g. stdv of observations)
        """

        if markersize is None:
            markersize = self.size * 0.67

        def nanmean(x):
            return np.mean(x[~np.isnan(x)])

        if self.normalize:
            S = S / normvalue

        show_shift = False
        if (R1 is not None) & (S1 is not None):
            if (len(R1) == len(S1)) & (len(R1) == len(R)):
                show_shift = True
            else:
                print len(R1), len(S1), len(R)
                sys.exit('Not possible to show shift in taylor as different sizes!')

        if markeredgecolor is None:
            markeredgecolor = markerfacecolor

        #/// calculate mean values ///
        try:
            self.R_mean = nanmean(R)
        except:
            self.R_mean = np.nan
        try:
            self.S_mean = nanmean(S)
        except:
            self.S_mean = np.nan
        if R1 is not None:
            self.R1_mean = nanmean(R1)
        else:
            self.R1_mean = None
        if S1 is not None:
            self.S1_mean = nanmean(S1)
        else:
            self.S1_mean = None

        #/// get coordinates
        x, y = self.map2xy(R, S)

        #/// generate plots
        if labels is None:
            self.ax.plot(x, y, marker=marker, linestyle=linestyle, color=color, markerfacecolor=markerfacecolor, linewidth=linewidth,
                         label=label, markeredgecolor=markeredgecolor, markersize=markersize, markeredgewidth=markeredgewidth)
        else:
            for i in range(len(x)):
                self.ax.text(x[i], y[i], labels[i], color=color, fontsize=markersize, horizontalalignment='center')

        if plot_mean:
            #mean of R/S data
            xmean, ymean = self.map2xy(self.R_mean, self.S_mean)
            self.ax.plot(xmean, ymean, marker=marker, linestyle=linestyle, color=color, markerfacecolor=markerfacecolor, linewidth=linewidth, label=label, markeredgecolor='blue', markersize=markersize*2, markeredgewidth=markeredgewidth*2.)

            if (self.R1_mean is not None) & (self.S1_mean is not None):
                xmean, ymean = self.map2xy(self.R1_mean, self.S1_mean)
                self.ax.plot(xmean, ymean, marker=marker, linestyle=linestyle, color=shiftcolor, markerfacecolor=markerfacecolor, linewidth=linewidth, label=label, markeredgecolor=shiftcolor, markersize=markersize*2, markeredgewidth=markeredgewidth*2)

        #/// shift plot
        if show_shift:
            x1, y1 = self.map2xy(R1, S1)
            #self.ax.plot(x1,y1,'ro')
            for i in np.arange(len(x1)):
                dx = x1[i]-x[i]
                dy = y1[i]-y[i]
                self.ax.arrow(x[i], y[i], dx, dy, edgecolor=shiftcolor, alpha=0.5, linewidth=2.)

                #this is a nice looking more flexible arrow
                #self.ax.annotate('', xy=(x[i]+dx, y[i]+dy),  xycoords='data',xytext=(x[i], y[i]), textcoords='data', \
                #arrowprops=dict(facecolor='grey', shrink=0.0,edgecolor='None',width=2), \
                #horizontalalignment='right', verticalalignment='top')

        #/// rescale axes
        if self.negative:
            axmin = -self.stdmax
        else:
            axmin = 0.
        self.ax.set_xlim(axmin, self.stdmax)
        self.ax.set_ylim(0., self.stdmax)

    def map2xy(self, R, S):
        """
        map R and S coordinates into x and y
        """

        theta = np.deg2rad((1.-R)*90.)
        r = S
        if self.r_equidistant:
            x = r*np.cos(theta)
            y = r*np.sin(theta)
        else:
            x = S*R
            y = S*np.sin(np.arccos(R))
        return x, y

    def set_legend_plot(self, label, color='black', marker='o'):
        """
        set a plot to be registered for the legend
        """
        self.labelnames.append(label)
        self.plots.append(self.ax.scatter([-1], [-1], color=color, marker=marker, s=8))  # fake plot for legend
        self.ax.set_xlim(0., self.stdmax)
        self.ax.set_ylim(0., self.stdmax)


def test():
    '''
    test taylor diagram
    this routine is for developing purposes
    '''
    #/// testing ///
    plt.close('all')
    #generate some sample data
    corr1 = np.asarray([1., .5, 0.9, 1., -0.8])
    corr2 = np.asarray([0.1, 0.25, 0.5, 0.6])
    s1 = np.asarray([1., 1., 1.5, 0.5])
    s2 = np.asarray([0.9, 1.5, 1.7])

    s1 = [0.5, .5, 1.5, 1.2]
    corr1 = [0.6, 1., 0.95, -0.8]
    S1 = [0.4, .2, 1.7, 0.7]
    R1 = [0.4, 0.9, 0.5, 0.6]

    s1 = np.asarray(s1)
    corr1 = np.asarray(corr1)
    S1 = np.asarray(S1)
    R1 = np.asarray(R1)

    #Taylor plot with equidistant correlations
    tay1 = taylor(r_equidistant=True)  # initialize the Taylor diagram
    tay1.plot(corr1, s1, markerfacecolor='green', marker='^', label='test1')  # plot some data
    tay1.ax.legend(loc='lower center', ncol=2, fancybox=True, shadow=True)  # if you like, access the axis object and specify legend properties
    tay1.ax.set_title('equidistant')

    plt.close('all')

    print 'corr1: ', corr1
    #Taylor plot with no equidistant correlations (as original in Taylor 2001)
    tay2 = taylor()  # initialize the Taylor diagram
    tay2.plot(corr1, s1, markerfacecolor='green', marker='^', label='test1', R1=R1, S1=S1)  # plot some data
    #tay2.plot_rms_meshlines()
    tay2.ax.legend(loc='lower center', ncol=2, fancybox=True, shadow=True)  # if you like, access the axis object and specify legend properties
    tay2.ax.set_title('original taylor')

    tay2.plot_taylor_legend()

    return tay2





#test()
#show()
