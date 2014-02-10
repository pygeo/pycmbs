# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please refer to
the pyCMBS licensing details.
"""

"""
This module implements generic map plotting capabilities
"""

installed_backends=[]

import os
try:
    from mpl_toolkits.basemap import Basemap
    installed_backends.append('basemap')
except:
    print('WARNING: BASEMAP seems not to be installed and can therefore not be used as plotting backend')

try:
    import cartopy.crs as ccrs
    installed_backends.append('cartopy')
except:
    print('WARNING: CARTOPY seems not to be installed and can therefore not be used as plotting backend')

try:
    from matplotlib import pyplot as plt
    installed_backends.append('imshow')
except:
    raise ValueError('Fatal error: You need to have a valid matplotlib installation to be able to work with pycmbs')

import matplotlib as mpl
import numpy as np
import matplotlib.gridspec as grd

from pycmbs.data import Data
from pycmbs.plots import ZonalPlot

class MapPlotGeneric(object):
    """
    Generic class to produce map plots
    """

    def __init__(self, backend='imshow', format='png', savefile=None,
                 show_statistic=True, stat_type='mean', figsize=(10, 5),
                 **kwargs):
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

        if 'ax' in kwargs.keys():
            ax = kwargs['ax']
            if ax is not None:
                self.ax_main = ax
            else:
                self._set_default_axis(figsize)  # set ax_main
        else:
            self._set_default_axis(figsize)  # set ax_main

        self.figure = self.ax_main.figure
        self._set_axis_invisible(self.ax_main, frame=False)

        # set default layout parameters for map
        self._set_layout_parameters()  # TODO make this dependent on the ax_main size!

        # consistency checks
        self._check()

        # set plotting backend routine
        if self.backend == 'imshow':
            self._draw = self._draw_imshow
        elif self.backend == 'basemap':
            self._draw = self._draw_basemap
        elif self.backend == 'cartopy':
            self._draw = self._draw_cartopy

        else:
            raise ValueError('Unknown backend!')

    def _set_default_axis(self, figsize):
        fig = plt.figure(figsize=figsize)
        # full size axis in case of new figure
        self.ax_main = fig.add_axes([0., 0., 1., 1.],
                                    label='ax_main')

    def _set_layout_parameters(self, left=0.1, right=0.9, bottom=0.1,
                               top=0.9, wspace=0.05, hspace=0.05,
                               wcolorbar=0.03, hcolorbar=0.04,
                               wzonal=0.1):
        """
        set layout parameters for map plotting
        all parameters are given in relative figure units

        If an axis 'ax_main' is already defined, then a scaling
        factor is calculated that transforms the coordinates in such
        a way that they are related to the ax_main position and size

        Parameters
        ----------
        left : float
            position of leftmost axis
        right : float
            position of the rightmost edge of all axes
        bottom : float
            position of bottom position of lowes axis
        top : float
            position of upper edge of uppermost axis
        wspace : float
            width space between different axes
        hspace : float
            height space between different axes
        wcolorbar : float
            width of the colorbar axis (only used when vertical oriented)
        hcolorbar : float
            height of the colorbar axis (only used when horizontal oriented)
        wzonal : float
            width of zonal mean plot
        """

        if self.ax_main is not None:
            # rescale dimensions to current ax_main properties
            b = self.ax_main.get_position()
            left = b.x0 + left * b.width
            right = b.x1 - (1. - right) * b.width
            bottom = b.y0 + bottom * b.height
            top = b.y1 - (1. - top) * b.height
            wspace = wspace * b.width
            hspace = hspace * b.height
            wcolorbar = wcolorbar * b.width
            hcolorbar = hcolorbar * b.height
            wzonal = wzonal * b.width

        self._layout_parameters = {'left': left, 'right': right,
                                   'bottom': bottom, 'top': top,
                                   'wspace': wspace, 'hspace': hspace,
                                   'wcolorbar': wcolorbar,
                                   'hcolorbar': hcolorbar,
                                   'wzonal': wzonal}

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
            self.savefile += tok + '.nc'
        else:
            self.savefile = self.savefile[:-3] + tok + '.nc'
        self.x.save(self.savefile, timmean=timmean, delete=True,
                    mean=False)

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
        if self.backend not in installed_backends:
            raise ValueError('Invalid plotting backend: %s' % self.backend)
        if self.backend == 'basemap':
            print('INFO: You have chosen BASEMAP as plotting backend. It is recommended to use CARTOPY instead as it is faster and also provides higher quality plotting capabilities.')

    def _draw_basemap(self, proj_prop=None, drawparallels=True, **kwargs):
        if proj_prop is None:
            raise ValueError('No projection properties are given! Please modify or choose a different backend!')

        the_map = Basemap(ax=self.pax, **proj_prop)
        xm = self.x.timmean()

        Z = xm
        lon = self.x.lon
        lat = self.x.lat

        X, Y = the_map(lon, lat)
        self.im = the_map.pcolormesh(X, Y, Z, **kwargs)

        self.__basemap_ancillary(the_map, drawparallels=drawparallels)

    def _draw_cartopy(self, proj_prop=None, **kwargs):
        if proj_prop is None:
            raise ValueError('No projection properties are given! Please modify or choose a different backend!')

        if proj_prop['projection'] == 'robin':
           pass
        else:
            raise ValueError('Unsupported projection type')

        xm = self.x.timmean()

        Z = xm
        lon = self.x.lon
        lat = self.x.lat

        # convert normal axis to GeoAxis
        def _ax2geoax(ax, ccrs_obj):
            """
            This routine converts a given matplotlib axis to a GeoAxis.
            It is ensured that the axis has the same position and size.

            Parameters
            ----------
            ax : axis
                matplotlib axis to be modified
            ccrs_obj : cartopy.crs
                reference system object

            Example
            -------
            ax2 = _ax2geoax(ax2, ccrs.Robinson())

            """

            b = ax.get_position()
            rect = [b.x0, b.y0, b.width, b.height]
            ax.set_visible(False)
            return ax.figure.add_axes(rect, label="pax", projection=ccrs_obj)

        self.pax = _ax2geoax(self.pax, ccrs.Robinson())

        # plot and ancillary plots
        self.pax.set_global()  # ensure global plot
        self.pax.coastlines()
        self.im = self.pax.pcolormesh(lon, lat, Z, transform=ccrs.PlateCarree(), **kwargs)
        self.pax.gridlines()


    def __basemap_ancillary(self, m, latvalues=None, lonvalues=None,
                            drawparallels=True, drawcountries=True,
                            land_color=0.8):
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
            latvalues = np.arange(-90., 120., 30.)
        if lonvalues is None:
            lonvalues = np.arange(-180., 180., 90.)
        if drawcountries:
            m.drawcountries()
        m.drawcoastlines()
        #~ m.drawlsmask(lakes=True,land_color=land_color)
        m.drawmapboundary()  # draw a line around the map region
        if drawparallels:
            m.drawparallels(latvalues, labels=[0, 0, 0, 0])
            m.drawmeridians(lonvalues, labels=[0, 0, 0, 0])

    def _draw_imshow(self, **kwargs):
        """
        draw data using imshow command
        """
        if self.pax is None:
            raise ValueError('Fatal Error: no axis for plotting specified')
        if self.pax is not None:
            self._set_axis_invisible(self.pax, frame=True)
        if self.cax is not None:
            self._set_axis_invisible(self.cax, frame=True)
        if self.zax is not None:
            self._set_axis_invisible(self.zax, frame=True)

        # do plotting
        self.im = self.pax.imshow(self.x.timmean(),
                                  interpolation='nearest', **kwargs)

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
            if l is None:
                return None
            else:
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

        cb = mpl.colorbar.ColorbarBase(self.cax, cmap=self.cmap,
                                       norm=self.norm,
                                       ticks=self._get_cticks(),
                                       orientation=self.colorbar_orientation)
        if self.ctick_prop is not None:
            l = self._get_cticklabels()
            if l is not None:
                cb.set_ticklabels(l)

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
        assert(isinstance(x, Data))
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
        if hasattr(self.cmap, 'monochrome'):
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
        self.pax.set_title(unit, loc='right', size=fontsize - 2)
        self.pax.set_title(stat, loc='left', size=fontsize - 2)

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
                st = st[0]
                s = 'mean: $' + str(round(me, 2)) + ' \pm ' + str(round(st, 2)) + '$'
            elif self.stat_type == 'sum':  # area sum
                me = tmp_xm.areasum()
                assert(len(me) == 1)
                me = me[0]
                s = 'sum: $' + str(round(me, 2)) + '$'
            else:
                me = np.ma.median(tmp_xm.data)
                s = 'median: $' + str(round(me, 2)) + '$'
        return s

    def _draw_zonal_plot(self, timmean=True, vmin=None, vmax=None,
                         fontsize=8):
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
            nt, ny, nx = self.x.shape

        ZP.plot(self.x, timmean=timmean, show_ylabel=False)

        # set limits
        if ((vmin is None) & (vmax is None)):
            vmin = self.zax.get_xlim()[0]
            vmax = self.zax.get_xlim()[1]
            # symmetry if neg. and positive limits
            if (vmin < 0.) & (vmax > 0.):
                val = max(abs(vmin), abs(vmax))
                vmin = -val
                vmax = val

        if vmin is None:
            vmin = self.zax.get_xlim()[0]
        if vmax is None:
            vmax = self.zax.get_xlim()[1]
        self.zax.set_xlim(vmin, vmax)

        # set only first and last label
        self.zax.set_xticks([vmin, vmax])
        self.zax.plot([0, 0], self.zax.get_ylim(), linestyle='-',
                      color='grey')

        for tick in self.zax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)

    def plot(self, show_zonal=False, show_histogram=False,
             show_timeseries=False, show_colorbar=True,
             colorbar_orientation='vertical', cmap='jet',
             ctick_prop=None,
             vmin=None, vmax=None, nclasses=10,
             title=None, proj_prop=None, drawparallels=True, titlefontsize=14):
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

        if colorbar_orientation not in ['vertical', 'horizontal']:
            raise ValueError('Invalid colorbar orientation')

        self.show_zonal = show_zonal
        self.colorbar_orientation = colorbar_orientation
        self.show_histogram = show_histogram
        self.show_timeseries = show_timeseries
        self.show_colorbar = show_colorbar
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
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, proj_prop=proj_prop, drawparallels=drawparallels)
        elif self.backend == 'cartopy':
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, proj_prop=proj_prop)
        else:
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap)

        # colorbar
        if self.show_colorbar:
            self._set_colorbar(self.im)
        self._plot_zonal()
        self._draw_title(title=title, fontsize=titlefontsize)

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

        # case 1: only the standard plot is desired
        if (not self.show_timeseries) and (not self.show_histogram) and (not self.show_colorbar):
            self._set_layout0()
            return

        if self.show_colorbar:
            # timeseries or histogram require an additional lower axis
            if (self.show_timeseries or self.show_histogram):
                raise ValueError('Combination with timeseries not supported yet!')
            else:
                if self.show_zonal:
                    self._set_layout2()  # layout with colorbar and zonal plot
                    return
                else:
                    self._set_layout1()  # layout with only colorbar
                    return
        else:
            raise ValueError('Layout without colorbar not supported yet')

    def _set_layout0(self):
        """
        set layout with only a plotting axis
        +----------+
        |  pax     |
        +----------+
        """
        left = self._layout_parameters['left']
        bottom = self._layout_parameters['bottom']
        right = self._layout_parameters['right']
        width = right - left
        height = self._layout_parameters['top'] - bottom
        self.pax = self.figure.add_axes([left, bottom, width, height], label="pax")

    def _set_layout1(self):
        """
        setlayout with only colorbar. This might be oriented either
        vertically or horizontally

        vertical
        +----------+ +-+
        |  pax     | |c|
        +----------+ +-+

        horizontal
        +---------+
        |  pax    |
        |         |
        +---------+
        +---------+
        |  cax    |
        +---------+
        """

        if not self.show_colorbar:
            raise ValueError('This routine was called by fault!')

        if self.colorbar_orientation == 'horizontal':
            # main plot axis
            left = self._layout_parameters['left']
            bottom = self._layout_parameters['bottom'] + self._layout_parameters['hspace']
            right = self._layout_parameters['right']
            width = right - left
            height = self._layout_parameters['top'] - bottom
            self.pax = self.figure.add_axes([left, bottom, width, height], label="pax")

            # colorbar axis
            bottom = self._layout_parameters['bottom']
            height = self._layout_parameters['hcolorbar']
            self.cax = self.figure.add_axes([left, bottom, width, height], label="cax")

        elif self.colorbar_orientation == 'vertical':
            # main plot axis
            left = self._layout_parameters['left']
            right = self._layout_parameters['right']
            bottom = self._layout_parameters['bottom']
            top = self._layout_parameters['top']
            wcolorbar = self._layout_parameters['wcolorbar']
            wspace = self._layout_parameters['wspace']
            width = right - left - wspace - wcolorbar
            height = top - bottom
            self.pax = self.figure.add_axes([left, bottom, width, height], label="pax")

            # colorbar axis
            left = right - wcolorbar
            rect = [left, bottom, wcolorbar, height]
            self.cax = self.figure.add_axes(rect, label="cax")

        else:
            raise ValueError('Invalid option')

    def _set_layout2(self):
        """
        layout with zonal mean and colorbar

        colorbar_orientation = 'vertical'
        +--+ +----------------+ +-+
        |  | |                | | |
        |  | |                | | |
        +--+ +----------------+ +-+

        colorbar_orientation = 'horizontal'
        +--+ +----------------+
        |  | |                |
        |  | |                |
        +--+ +----------------+
             +----------------+
             |                |
             +----------------+

        """
        if not self.show_zonal:
            raise ValueError('Only WITH zonal mean supported here!')
        if not self.show_colorbar:
            raise ValueError('Only WITH colorbar supported here!')

        if self.colorbar_orientation == 'vertical':
            # zonal axis
            left = self._layout_parameters['left']
            right = left + self._layout_parameters['wzonal']
            zwidth = right - left

            height = self._layout_parameters['top'] - self._layout_parameters['bottom']
            bottom = self._layout_parameters['bottom']

            rect = [left, bottom, zwidth, height]
            self.zax = self.figure.add_axes(rect, label="zax")

            # main plotting axis
            left = self._layout_parameters['left'] + zwidth + self._layout_parameters['wspace']
            right = self._layout_parameters['right'] - self._layout_parameters['wspace'] - self._layout_parameters['wcolorbar']
            rect = [left, bottom, right - left, height]
            self.pax = self.figure.add_axes(rect, label="pax")

            # colorbar axis
            left = self._layout_parameters['right'] - self._layout_parameters['wcolorbar']
            right = self._layout_parameters['right']
            rect = [left, bottom, right - left, height]
            self.cax = self.figure.add_axes(rect, label="cax")

        elif self.colorbar_orientation == 'horizontal':
            # colorbar axis
            left = self._layout_parameters['left'] + self._layout_parameters['wzonal'] + self._layout_parameters['wspace']
            right = self._layout_parameters['right']
            bottom = self._layout_parameters['bottom']
            height = self._layout_parameters['hcolorbar']
            rect = [left, bottom, right - left, height]
            self.cax = self.figure.add_axes(rect, label="cax")

            # main plotting axis
            left = self._layout_parameters['left'] + self._layout_parameters['wzonal'] + self._layout_parameters['wspace']
            right = self._layout_parameters['right']
            bottom = self._layout_parameters['bottom'] + self._layout_parameters['hspace'] + self._layout_parameters['hcolorbar']
            top = self._layout_parameters['top']
            rect = [left, bottom, right - left, top - bottom]
            self.pax = self.figure.add_axes(rect, label="pax")

            # zonal axis
            left = self._layout_parameters['left']
            right = self._layout_parameters['left'] + self._layout_parameters['wzonal']
            rect = [left, bottom, right - left, top - bottom]
            self.zax = self.figure.add_axes(rect, label="zax")

        else:
            raise ValueError('Invalid colorbar option')


class MultipleMap(MapPlotGeneric):
    def __init__(self, geometry=None, **kwargs):
        """
        Parameters
        ----------
        geometry : tuple
            specifies number of subplots (nrows,mcols); e.g. geometry=(2,3)
            give a plot with 2 rows and 3 cols.
        """
        if geometry is None:
            raise ValueError('You need to specify a valid geometry.')
        raise ValueError('MultipleMap not further implemented yet.')
        assert(isinstance(geometry, tuple))
        self.n = geometry[0]
        self.m = geometry[1]
        super(MultipleMap, self).__init__(**kwargs)
        self.fig = plt.figure()
        stop
        # TODO : define axes here !


def map_plot(x, use_basemap=False, show_zonal=False,
             show_histogram=False, show_timeseries=False,
             show_stat=False, stat_type='mean', savefile=None,
             nclasses=10, colorbar_orientation='vertical',
             show_colorbar=True, cmap_data='jet', title=None,
             vmin=None, vmax=None, proj='robin', lon_0=0., lat_0=0.,
             cticks=None, cticklabels=None, ax=None,
             drawparallels=True, overlay=None, titlefontsize=14,
             zonal_timmean=None, region=None, savegraphicfile=None):
    """
    This is a wrapper function to replace the old map_plot routine
    It provides a similar interface, but makes usage of the new
    SingleMap object for plotting

    Parameters
    ----------
    x : Data
        data object to be plotted
    """

    print('WARNING: usage of map_plot is depreciated. This routine will be removed in future versions. Please use SingleMap instead.')

    # unused keyword arguments! These are not further handled and are
    # just for compatability reasons with old map_plot routine
    if overlay is not None:
        print('WARNING: overlay function not supported yet!')  # TODO
    if zonal_timmean is not None:
        if not zonal_timmean:
            print('WARNING: zonal_timmean option not supported yet!')  # TODO
    if region is not None:
        print('WARNING: region option not supported yet!')  # TODO

    if use_basemap:
        backend = 'basemap'
    else:
        backend = 'imshow'

    proj_prop = {'projection': proj, 'lon_0': lon_0, 'lat_0': lat_0}
    ctick_prop = {'ticks': cticks, 'labels': cticklabels}

    M = SingleMap(x, backend=backend, show_statistic=show_stat,
                  stat_type=stat_type, savefile=savefile, ax=ax)
    M.plot(title=title, show_zonal=show_zonal, show_histogram=False,
           show_timeseries=False, nclasses=nclasses,
           colorbar_orientation=colorbar_orientation,
           show_colorbar=show_colorbar, cmap=cmap_data, vmin=vmin,
           vmax=vmax, proj_prop=proj_prop, ctick_prop=ctick_prop,
           drawparallels=drawparallels, titlefontsize=titlefontsize)

    if savegraphicfile is not None:  # save to graphics file
        if os.path.exists(savegraphicfile):
            os.remove(savegraphicfile)
        M.figure.savefig(savegraphicfile, dpi=200)


    # TODO
    #arguments from original map plot which are not covered yet
    #           ~  ax=None, , ,
     #~ ,
     #~ , regions_to_plot=None, logplot=False,
     #~ logoffset=None,
     #~ f_kdtree=False, latvalues=None,
     #~ lonvalues=None
     #~
     #~ scal_timeseries=1., vmin_zonal=None, vmax_zonal=None,
     #~ bluemarble=False, contours=False,
     #~, drawcountries=True,
     #~
     #~ contourf=False, land_color=(0.8, 0.8, 0.8),
     #~ regionlinewidth=1, bins=10,
     #~
     #~ cax_rotation=0.,
     #~ plot_method='colormesh', boundinglat=60.,
     #~  , , savegraphicfile=None,
     #~ **kwargs):
