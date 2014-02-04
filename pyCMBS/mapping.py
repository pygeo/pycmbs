# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""

"""
This module implements generic map plotting capabilities
"""

import matplotlib as mpl
mpl.use('agg')
import os
from mpl_toolkits.basemap import Basemap,shiftgrid

from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as grd

from . data import Data
from . plots import ZonalPlot

class MapPlotGeneric(object):
    """
    Generic class to produce map plots
    """

    def __init__(self, backend='imshow', format='png', savefile=None,
                    show_statistic = True, stat_type='mean', figure=None):

        """
        Parameters
        ----------
        backend : str
            specifies plotting backend for map plotting. Currently the
            following backends are supported
            'imshow' : =default, simply uses plt.imshow for vizualization
                       this is the fastest approach
            'basemap' : uses Basemap as backend for plotting
        """
        self.backend = backend
        self.format = format
        self.savefile = savefile
        self.show_statistic = show_statistic
        self.stat_type = stat_type
        self._dummy_axes=[]

        if figure is None:
            self.figure = plt.figure()
        else:
            self.figure = figure

        # consistency checks
        self._check()

        # set plotting backend routine
        if self.backend == 'imshow':
            self._draw = self._draw_imshow
        elif self.backend == 'basemap':
            self._draw = self._draw_basemap
        else:
            raise ValueError('Unknown backend!')

    def _save_data_to_file(self, timmean=True):
        """
        save data to file. The filename is automatically generated
        depednend on self.savefile and the option timmean

        Parameters
        ----------
        timmean : bool
            if True, then the temporal mean field (like in the plot) is
            stored. Otherwise the entire dataset which was available for
            plotting will be stored
        """
        if timmean:
            tok = '_timmean'
        else:
            tok = '_all'
        if self.savefile is None:
            return
        if self.savefile[:-3] != '.nc':
            self.savefile +=  tok + '.nc'
        else:
            self.savefile = self.savefile[:-3] + tok +  '.nc'
        self.x.save(self.savefile, timmean=timmean, delete=True, mean=False)

    def save(self, save_mean=True, save_all=False):
        """
        save data to file

        Parameters
        ----------
        save_mean : bool
            save temporal mean field to file
        save_all : bool
            save entire field which was available for plotting
            to file
        """
        if save_mean:
            self._save_data_to_file(timmean=True)
        if save_all:
            self._save_data_to_file(timmean=False)

    def _check(self):
        if self.stat_type not in ['mean', 'median', 'sum']:
            raise ValueError('Invalid statistic type: %s' % self.stat_type)
        if self.backend not in ['imshow','basemap']:
            raise ValueError('Invalid plotting backend: %s' % self.backend)

    def _draw_basemap(self, proj_prop=None, **kwargs):
        if proj_prop is None:
            raise ValueError('No projection properties are given! Please modify or choose a different backend!')

        the_map = Basemap(ax=self.pax, **proj_prop)
        xm = self.x.timmean()

        Z = xm
        lon = self.x.lon
        lat = self.x.lat

        X, Y = the_map(lon, lat)
        self.im = the_map.pcolormesh(X, Y, Z, **kwargs)

        self.__basemap_ancillary(the_map)


    def __basemap_ancillary(self, m, latvalues=None, lonvalues=None, drawparallels=True, drawcountries=True, land_color=0.8):
        """
        routine to plot ancillary data like coastlines
        or meridians on a basemap plot

        Parameters
        ----------
        m : Basemap
            map to add features to

        latvalues : ndarray
            latitude values for drawing grid (optional)

        lonvalues : ndarray
            longitude values for drawing grid (optional)

        """

        if latvalues is None:
            latvalues=np.arange(-90., 120., 30.)
        if lonvalues is None:
            lonvalues= np.arange(-180., 180., 90.)
        if drawcountries:
            m.drawcountries()
        m.drawcoastlines()
        #~ m.drawlsmask(lakes=True,land_color=land_color)
        m.drawmapboundary() # draw a line around the map region
        if drawparallels:
            m.drawparallels(latvalues,labels=[0, 0, 0, 0])
            m.drawmeridians(lonvalues,labels=[0, 0, 0, 0])

    def _draw_imshow(self, **kwargs):
        """
        draw data using imshow command
        """
        if self.pax is None:
            raise ValueError('Fatal Error: no axis for plotting specified')
        # set dummy axes invisible
        for ax in self._dummy_axes:
            self._set_axis_invisible(ax, frame=False)
        if self.pax is not None:
            self._set_axis_invisible(self.pax, frame=True)
        if self.cax is not None:
            self._set_axis_invisible(self.cax, frame=True)
        if self.zax is not None:
            self._set_axis_invisible(self.zax, frame=True)

        # do plotting
        self.im = self.pax.imshow(self.x.timmean(), interpolation='nearest', **kwargs)

    def _get_cticks(self):
        """ get cticks from dictionary if available """
        if self.ctick_prop is None:
            return None
        if 'ticks' in self.ctick_prop.keys():
            return self.ctick_prop['ticks']

    def _get_cticklabels(self):
        if self.ctick_prop is None:
            return None
        if 'labels' in self.ctick_prop.keys():
            l = self.ctick_prop['labels']
            if len(l) != len(self._get_cticks()):
                raise ValueError('CTICKS and CTICKLABELS need to have the same length!')
            return l
        else:
            return None




    def _set_colorbar(self, im):
        """
        create colorbar and return colorbar object

        Parameters
        ----------
        im : plot
            results from e.g. an imshow command
        """
        if not self.show_colorbar:
            raise ValueError('Colorbar can not be generated when not requested')

        vmin = im.get_clim()[0]
        vmax = im.get_clim()[1]
        self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        cb   = mpl.colorbar.ColorbarBase(self.cax, cmap=self.cmap, norm=self.norm, ticks=self._get_cticks(), orientation=self.colorbar_orientation)
        if self.ctick_prop is not None:
            cb.set_ticklabels(self._get_cticklabels())

    def _set_axis_invisible(self, ax, frame=True):
        """
        set axis parameteres in a way that it is invisible
        """
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_frame_on(frame)





class SingleMap(MapPlotGeneric):
    """
    A class to generate a plot with a single figure
    """
    def __init__(self, x, **kwargs):
        """
        Parameters
        ----------
        X : Data
            Data object with data to plot
        """
        assert(isinstance(x,Data))
        super(SingleMap, self).__init__(**kwargs)
        self.x = x

        self.pax = None  # axis for plot
        self.cax = None  # axis for colorbar
        self.tax = None  # axis for timeseries
        self.hax = None  # axis for histogram
        self.zax = None  # axis for zonal plot

        self.cmap = 'jet'

    def _plot_zonal(self):
        if self.show_zonal:
            if self.x._latitudecheckok:
                self._draw_zonal_plot(self)
            else:
                print('WARNING: zonal plot not possible due to invalid latitude configurations')


    def _set_cmap(self, nclasses):
        """
        generate a colormap. If self.cmap is already a colormap
        object, then nothing happens. Otherwise a new colormap object
        is created which has nclasses

        Parameters
        ----------
        nclasses : int
            number of classes for colormap
        """
        if hasattr(self.cmap,'monochrome'):
            # colormap object was given
            self.cmap = cmap_data
        else:
            self.cmap = plt.cm.get_cmap(self.cmap, nclasses)


    def _draw_title(self, title=None, fontsize=14):
        """
        draw title, units and statistics

        Parameters
        ----------
        title : str
            title for figure. If not specified, then the label
            of the data will be used
        fontsize : int
            fontsize for the title

        """
        stat = self._get_statistics_str()
        if title is None:
            title = self.x._get_label()
        unit = self.x._get_unit()

        self.pax.set_title(title + '\n', size=fontsize)
        self.pax.set_title(unit, loc='right', size=fontsize-2)
        self.pax.set_title(stat, loc='left', size=fontsize-2)

    def _get_statistics_str(self):
        tmp_xm = self.x.timmean(return_object=True)  # from temporal mean
        s = ''
        if self.show_statistic:
            if self.stat_type == 'mean':
                me = tmp_xm.fldmean()
                st = tmp_xm.fldstd()
                assert(len(me) == 1)
                assert(len(st) == 1)
                me = me[0]
                st=st[0]
                s ='mean: $' + str(round(me,2))  + ' \pm ' + str(round(st,2)) + '$'
            elif stat_type == 'sum': #area sum
                me = tmp_xm.areasum()
                assert(len(me) == 1)
                me = me[0]
                s = 'sum: $' + str(round(me,2))  + '$'
            else:
                me = np.ma.median(tmp_xm.data)
                s = 'median: $' + str(round(me,2)) + '$'
        return s


    def _draw_zonal_plot(self, timmean=True, vmin=None, vmax=None, fontsize=8):
        """
        calculate zonal statistics and add to zonal axis

        Parameters
        ----------
        timmean : bool
            temporal mean for zonal plot [default=True]
        vmin : float
            minimum value for zonal plot
        vmax : float
            maximum value for zonal plot
        """

        ZP = ZonalPlot(ax=self.zax, dir='y')

        if self.x.ndim == 2:
            pass
        elif self.x.ndim == 3:
            nt,ny,nx = self.x.shape

        ZP.plot(self.x, timmean=timmean, show_ylabel=False)

        # set limits
        if ((vmin is None) & (vmax is None)):
            vmin = self.zax.get_xlim()[0]
            vmax = self.zax.get_xlim()[1]
            # symmetry if neg. and positive limits
            if (vmin < 0.) & (vmax > 0.):
                val = max(abs(vmin) ,abs(vmax))
                vmin = -val
                vmax = val

        if vmin is None:
            vmin = self.zax.get_xlim()[0]
        if vmax is None:
            vmax = self.zax.get_xlim()[1]
        self.zax.set_xlim(vmin, vmax)

        # set only first and last label
        self.zax.set_xticks([vmin, vmax])
        self.zax.plot([0,0], self.zax.get_ylim(), linestyle='-', color='grey')

        for tick in self.zax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)



    def plot(self, show_zonal=False, show_histogram=False,
            show_timeseries=False, show_colorbar=True,
            colorbar_orientation='vertical', cmap='jet', ctick_prop=None,
            vmin=None, vmax=None, nclasses=10,
            title=None, proj_prop=None):
        """
        routine to plot a single map

        Parameters
        ----------

        proj_prop : dict
            dictionary for projection properties. Is needed if projected
            maps shall be drawn. The projection parameters given need
            to be compliant with the backend used for plotting


        ctick_prop : dict
            dictionary that specifies the properties for the colorbar
            ticks. Currently the following keys are supported:

            'ticks' : float list : specifies locations of ticks
            'labels' : str list : user defined label list; needs to
                                  have same length as 'ticks'

             Example:
                ctick_prop={'ticks':[-15, 0., 3.], 'labels':['A','B','C']
        """

        if colorbar_orientation not in ['vertical','horizontal']:
            raise ValueError('Invalid colorbar orientation')

        self.show_zonal=show_zonal
        self.colorbar_orientation=colorbar_orientation
        self.show_histogram=show_histogram
        self.show_timeseries=show_timeseries
        self.show_colorbar=show_colorbar
        self.ctick_prop = ctick_prop  # dictionary
        self.vmin = vmin
        self.vmax = vmax


        # set colormap and ensure to have a colormap object
        self.cmap = cmap
        self._set_cmap(nclasses)

        # set axes layout
        self._set_layout()

        # do plot using current backend
        if self.backend == 'basemap':
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, proj_prop=proj_prop)
        else:
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap)

        # colorbar
        if self.show_colorbar:
            self._set_colorbar(self.im)
        self._plot_zonal()
        self._draw_title(title=title)

        # adjust plots to minimize spaces between subplots
        self._adjust_figure()

        # save data if required
        self.save()

    def _adjust_figure(self):
        """
        adjust subplot sizes
        """
        # ensure that full space is covered by data
        self.pax.set_aspect('auto', adjustable='datalim')

        #~ if self.colorbar_orientation == 'vertical':
            #~ # pos = [left, bottom, width, height]
#~
            #~ cleft, cbottom, cright, ctop
            #~ res = self.cax.get_position().get_points()
            #~ cleft = res[0,0]
            #~ cbottom = res[0,1]
            #~ cright = res[1,0]
            #~ ctop = res[1,1]
            #~ cax_width = cright-cleft
#~
            #~ print cleft, cbottom, cright, ctop
#~
            #~ res = self.pax.get_position().get_points()
            #~ pleft = res[0,0]
            #~ pbottom = res[0,1]
            #~ pright = res[1,0]
            #~ ptop = res[1,1]
            #~ pheight = ptop-pbottom
            #~ pos = [cleft, pbottom, cax_width, pheight]
            #~ print 'pos: ', pos
            #~ self.cax.set_position(pos)
        #~ self.figure.tight_layout(w_pad=0., h_pad=0.)


    def _set_layout(self):
        """
        routine specifies layout of different axes
        """
        # check if option combinations are possible
        if self.show_timeseries and self.show_histogram:
            raise ValueError('Combination of histogram and timeseries not supported')

        if self.show_colorbar:
            # timeseries or histogram require an additional lower axis
            if (self.show_timeseries or self.show_histogram):
                raise ValueError('Combination with timeseries not supported yet!')
            else:
                if self.show_zonal:
                    self._set_layout2()  # layout with colorbar and zonal plot
                else:
                    self._set_layout1()  # layout with only colorbar

        else:
            raise ValueError('Layout without colorbar not supported yet')




    def _set_layout1(self):
        """
        setlayout with only colorbar. This might be oriented either
        vertically or horizontally

        vertical
        ----------- +-+
        |  pax      |c|
        ----------- +-+

        horizontal
        -----------
        |  pax    |
        |         |
        -----------
        -----------
        |  cax    |
        -----------

        """

        wspace=0.05

        if not self.show_colorbar:
            raise ValueError('This routine was called by fault!')
        if self.colorbar_orientation == 'horizontal':
            self.gs = grd.GridSpec(2, 1, height_ratios=[95,5], wspace=wspace)
        elif self.colorbar_orientation == 'vertical':
            self.gs = grd.GridSpec(1, 2, width_ratios=[95,5], wspace=wspace)
        else:
            raise ValueError('Invalid option')
        self.pax = self.figure.add_subplot(self.gs[0])
        self.cax = self.figure.add_subplot(self.gs[1])


    def _set_layout2(self):
        """
        layout with zonal mean and colorbar
        """
        if not self.show_zonal:
            raise ValueError('Only WITH zonal mean supported here!')
        if not self.show_colorbar:
            raise ValueError('Only WITH colorbar supported here!')

        if self.colorbar_orientation == 'horizontal':
            self.gs = grd.GridSpec(2, 2, height_ratios=[95,5], width_ratios = [15, 85], wspace=0.05)
            self.zax = self.figure.add_subplot(self.gs[0])
            self.pax = self.figure.add_subplot(self.gs[1])
            self.cax = self.figure.add_subplot(self.gs[3])
            self._dummy_axes.append(self.figure.add_subplot(self.gs[2]))
        elif self.colorbar_orientation == 'vertical':
            self.gs = grd.GridSpec(1, 3, width_ratios = [15, 80, 5], wspace=0.05)
            self.zax = self.figure.add_subplot(self.gs[0])
            self.pax = self.figure.add_subplot(self.gs[1])
            self.cax = self.figure.add_subplot(self.gs[2])
        else:
            raise ValueError('Invalid colorbar option')





class MultipleMap(MapPlotGeneric):
    def __init__(self,geometry=None,**kwargs):
        """
        Parameters
        ----------
        geometry : tuple
            specifies number of subplots (nrows,mcols); e.g. geometry=(2,3)
            give a plot with 2 rows and 3 cols.
        """
        assert(isinstance(geometry,tuple))
        self.n = geometry[0]
        self.m = geometry[1]
        super(MultipleMap,self).__init__(**kwargs)
        self.fig = plt.figure()
        stop
        #todo: define axes here !
