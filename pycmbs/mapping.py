# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
This module implements generic map plotting capabilities
"""

installed_backends = []

import os
try:  # note that this import statement needs to come BEFOREbasemap, as otherwise some cartopy libraries are not found for some strange reasons (at least on some machines)
    import cartopy.crs as ccrs
    installed_backends.append('cartopy')
except:
    print('WARNING: CARTOPY seems not to be installed and can therefore not be used as plotting backend')

try:
    from mpl_toolkits.basemap import Basemap
    installed_backends.append('basemap')
except:
    print('WARNING: BASEMAP seems not to be installed and can therefore not be used as plotting backend')


try:
    from matplotlib import pyplot as plt
    installed_backends.append('imshow')
except:
    raise ValueError('Fatal error: You need to have a valid matplotlib installation to be able to work with pycmbs')

import matplotlib as mpl
import numpy as np

import matplotlib.gridspec as grd
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

from pycmbs.data import Data
from pycmbs.plots import ZonalPlot
from pycmbs.polygon_utils import Polygon


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
            if 'cartopy' in installed_backends:
                # use cartopy instead of Basemap, when possible
                print('INFO: The backend has been automatically switched to CARTOPY as this provides higher quality and faster plotting on your machine')
                self.backend = 'cartopy'
            else:
                print('INFO: You have chosen BASEMAP as plotting backend. It is recommended to use CARTOPY instead as it is faster and also provides higher quality plotting capabilities.')

    def _draw_basemap(self, proj_prop=None, drawparallels=True, vmin_polygons=None, vmax_polygons=None, **kwargs):
        """
        """
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

        # add polygons to map
        if self.polygons is not None:
            if False:  # individual polygons
                for p in self.polygons:
                    self._add_single_polygon_basemap(the_map, p)
            else:  # plot all polygons at once
                self._add_polygons_as_collection_basemap(the_map, vmin=vmin_polygons, vmax=vmax_polygons)

    def _add_polygons_as_collection_basemap(self, plot_handler, **kwargs):
        collection = self._polygons2collection(plot_handler=plot_handler, **kwargs)
        self._add_collection(collection)

    def _add_single_polygon_basemap(self, m, p, color='red', linewidth=1):
        """
        plot region r on top of basemap map m

        Parameters
        ----------
        m : map
            Basemap object
        p : Polygon
            Polygon object as defined in polygon.py. Note that
            this is different from the matpltlib.Polygon object
        color : str
            color to plot region
        linewidth : int
            width of outline border for polygon to plot
        """
        from matplotlib.patches import Polygon as mplPolygon

        lons = p._xcoords()
        lats = p._ycoords()

        x, y = m(lons, lats)
        xy = list(zip(x, y))
        mapboundary = mplPolygon(xy, edgecolor=color, linewidth=linewidth, fill=False)
        self.pax.add_patch(mapboundary)

    # convert normal axis to GeoAxis
    def _ax2geoax(self, ax, ccrs_obj):
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

    def _draw_cartopy(self, proj_prop=None, vmin_polygons=None, vmax_polygons=None, **kwargs):
        if proj_prop is None:
            raise ValueError('No projection properties are given! Please modify or choose a different backend!')

        if proj_prop['projection'] in ['robin','TransverseMercator', 'mercator']:
            pass
        else:
            raise ValueError('Unsupported projection type')

        xm = self.x.timmean()
        Z = xm
        lon = self.x.lon
        lat = self.x.lat

        if proj_prop['projection'] == 'robin':
            act_ccrs = ccrs.Robinson()
        elif proj_prop['projection'] == 'TransverseMercator':
            act_ccrs = ccrs.TransverseMercator(central_longitude=proj_prop.pop('central_longitude', 0.), central_latitude=proj_prop.pop('central_latitude', 0.))
        elif proj_prop['projection'] == 'mercator':
            if 'extent' in proj_prop.keys():
                ymin = proj_prop['extent']['ymin']
                ymax = proj_prop['extent']['ymax']
            else:
                raise ValueError('Need to specify extent!')
            act_ccrs = ccrs.Mercator(central_longitude=proj_prop.pop('central_longitude', 0.), min_latitude=ymin, max_latitude=ymax)
        else:
            raise ValueError('Unsupported projection')

        self.pax = self._ax2geoax(self.pax, act_ccrs)

        # add cyclic coordinates if possible
        if self.x._equal_lon():
            try:
                lon1, lat1, Z1 = self._add_cyclic_to_field(self.x._get_unique_lon(), lat, Z)
            except:
                lon1 = None
            if lon1 is not None:
                lon = lon1
                lat = lat1
                Z = Z1

        # plot and ancillary plots
        if 'extent' in proj_prop.keys():
            if proj_prop['projection'] == 'mercator':
                pass
            else:
                xmin = proj_prop['extent']['xmin']
                xmax = proj_prop['extent']['xmax']
                ymin = proj_prop['extent']['ymin']
                ymax = proj_prop['extent']['ymax']

                self.pax.set_extent([xmin, xmax, ymin, ymax])
        else:
            self.pax.set_global()  # ensure global plot
        self.pax.coastlines()
        self.im = self.pax.pcolormesh(lon, lat, Z, transform=ccrs.PlateCarree(), **kwargs)
        self.pax.gridlines()  #draw_labels=kwargs.pop('draw_labels', True))

        # plot polygons
        if self.polygons is not None:
            if False:  # plot all polygons individually
                for p in self.polygons:
                    self._add_single_polygon_cartopy(p)
            else:  # all polygons as collection
                self._add_polygons_as_collection_cartopy(act_ccrs, vmin=vmin_polygons, vmax=vmax_polygons)

    def _add_collection(self, collection):
        if self.backend == 'imshow':
            raise ValueError('Collections not tested yet with backend IMSHOW')
        elif self.backend == 'cartopy':
            self.pax.add_collection(collection)
        elif self.backend == 'basemap':
            self.pax.add_collection(collection)
        else:
            raise ValueError('INVALID backend!')

    def _add_polygons_as_collection_cartopy(self, plot_handler, **kwargs):
        """
        add polygons as collection
        """
        collection = self._polygons2collection(plot_handler=plot_handler, **kwargs)
        self._add_collection(collection)

    def _polygons2collection(self, vmin=None, vmax=None, color='red', cmap='jet', plot_handler=None):
        """
        generate collection from list of polygons

        Parameters
        ----------
        color : str
            color for edges of polygons
        vmin : float
            minimum for scaling
        vmax : float
            maximum for scaling
        cmap : str, colormap object
            colormap specification
        """

        if plot_handler is None:
            raise ValueError('No plotting handler provided!')

        Path = mpath.Path
        patches = []
        pdata = np.ones(len(self.polygons)) * np.nan

        cnt = 0

        def _get_codes(Path, n):
            """
            specify how vertices are interconnected (here simple connection by lines)
            """
            codes = [Path.MOVETO]
            for i in xrange(n-1):
                codes.append(Path.LINETO)
            return codes


        for p in self.polygons:

            # convert lat/lon to map coordinates
            x, y = self._get_map_coordinates(p._xcoords(), p._ycoords(), plot_handler=plot_handler)
            x.append(x[0])
            y.append(y[0])
            verts = np.asarray([x, y]).T

            codes = _get_codes(Path, len(verts))

            # construct object and append to library of objects
            path = mpath.Path(verts, codes, closed=True)
            patches.append(mpatches.PathPatch(path))

            # store data information
            if p.value is not None:
                pdata[cnt] = p.value

            cnt += 1

        pdata = np.asarray(pdata)
        pdata = np.ma.array(pdata, mask=np.isnan(pdata))

        # generate collection
        if vmin is None:
            vmin = pdata.min()
        if vmax is None:
            vmax = pdata.max()

        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        collection = PatchCollection(patches, cmap=cmap, norm=norm, alpha=1., match_original=False, edgecolors=color)
        collection.set_array(pdata)

        return collection

    def _get_map_coordinates(self, lons, lats, plot_handler=None):
        if self.backend == 'imshow':
            print ValueError('Not implemented for backend IMSHOW')
        elif self.backend == 'basemap':
            if plot_handler is None:
                raise ValueError('ERROR: no basemap handler provided!')
            x, y = plot_handler(lons, lats)
            return list(x), list(y)
        elif self.backend == 'cartopy':
            X = plot_handler.transform_points(ccrs.PlateCarree(), lons, lats)  # gives an array (N, 3) with x,y,z as columns
            x = X[:,0]
            y = X[:,1]
            return list(x), list(y)
        else:
            raise ValueError('Invalid backend!')

    def _add_single_polygon_cartopy(self, p0, color='red', linewidth=1):
        """
        add a polygon to a map

        Parameters
        ----------
        p : Polygon, dict
            this is a Polygon object from polygon_utils or a dictionary
            if a dictionary is provided, this needs to have the following structure
            {'id' : int, 'polygon' : Polygon, 'color' : str}
            This allows to specifiy properties for each polygon individually
        """
        from matplotlib.patches import Polygon as mplPolygon

        if isinstance(p0, dict):
            p = p0['polygon']
            color = p0['color']
        else:
            p = p0

        lons = list(p._xcoords())
        lats = list(p._ycoords())

        lons.append(lons[0])
        lats.append(lats[0])
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        self.pax.plot(lons, lats, transform=ccrs.PlateCarree(), color=color, linewidth=linewidth)  # TODO create a polygon that might be also filled

    def _add_cyclic_to_field(self, lon, lat, z):
        """
        add an additional column to a dataset to avoid plotting problems
        around the 0degree longitude

        Parameters
        ----------
        lon : ndarray
            VECTOR of unique longitudes. Note that this needs to be a vector!
        lat : ndarray
            2D array of latitudes with same geometry than the data field Z
        z : ndarray
            2D array of values

        Returns
        -------
        lon : ndarray
            2D
        lat : ndarray
            2D
        Z : ndarray
            2D

        References
        ----------
        [1] https://github.com/SciTools/cartopy/issues/393
        [2] http://stackoverflow.com/questions/21864512/cartopy-behavior-when-plotting-projected-data
        [3] https://github.com/SciTools/cartopy/pull/394
        """
        try:
            from cartopy import util as ut
        except:
            print('Longitude shift can not be performed as most recent CARTOPY version seems not to be installed')
            return None, None, None
        assert lon.ndim == 1
        assert lat.shape == z.shape

        lat_out, lon1 = ut.add_cyclic_point(lat, coord=lon)
        z_out, lon1 = ut.add_cyclic_point(z, coord=lon)
        lon_out = np.ones_like(lat_out) * lon1
        return lon_out, lat_out, z_out

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
            self.cmap = self.cmap
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
             title=None, proj_prop=None, drawparallels=True,
             titlefontsize=14, polygons=None, vmin_polygons=None, vmax_polygons=None):
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
s
             Example:
                ctick_prop={'ticks':[-15, 0., 3.], 'labels':['A','B','C']
        polygons : list
            list of Polygon object of e.g. a regions to draw
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
        self.polygons = polygons

        # set colormap and ensure to have a colormap object
        self.cmap = cmap
        self._set_cmap(nclasses)

        # set axes layout
        self._set_layout()

        # do plot using current backend
        if self.backend == 'basemap':
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, proj_prop=proj_prop, drawparallels=drawparallels, vmin_polygons=vmin_polygons, vmax_polygons=vmax_polygons)
        elif self.backend == 'cartopy':
            self._draw(vmin=self.vmin, vmax=self.vmax, cmap=self.cmap, proj_prop=proj_prop, vmin_polygons=vmin_polygons, vmax_polygons=vmax_polygons)
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

    def plot_around_coordinate(self, clon, clat, radius, show_center=True, **kwargs):
        """
        generate map plot around a specific location
        The projection properties are ignored if they are provided
        in kwargs

        Parameters
        ----------
        clon : float
            center longitude coordinate [deg]
        clat : float
            center latitude coordinate [deg]
        show_center : bool
            show center coordinates by a marker
        """

        # update projection parameters
        proj_prop = kwargs.pop('proj_prop', None)
        proj_prop = {}

        xmin = clon - radius  # todo how to handle across dateline or near the poles
        xmax = clon + radius
        ymin = clat - radius
        ymax = clat + radius

        proj_prop.update({'projection' : 'TransverseMercator'})
        proj_prop.update({'central_longitude' : clon})
        proj_prop.update({'central_latitude' : clat})
        proj_prop.update({'extent' : {'xmin' : xmin, 'xmax' : xmax, 'ymin' : ymin, 'ymax' : ymax}})

        # generate Polygon to plot center location that gives a circle
        if show_center:
            theta = np.linspace(0., 2.*np.pi, 360)
            r = 0.02 * radius  # size as function of overall radius
            x = clon + r*np.cos(theta)
            y = clat + r*np.sin(theta)
            P = Polygon(1, zip(x,y))

        if 'polygons' in kwargs.keys():
            polygons = kwargs.pop('polygons')
            if polygons is None:
                polygons = [P]
            else:
                polygons.append(P)
        else:
            polygons = [P]

        self.plot(proj_prop=proj_prop, polygons=polygons, **kwargs)


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
             zonal_timmean=None, region=None, savegraphicfile=None, return_plot_handler=False, logplot=False):
    """
    This is a wrapper function to replace the old map_plot routine
    It provides a similar interface, but makes usage of the new
    SingleMap object for plotting

    Parameters
    ----------
    x : Data
        data object to be plotted
    """

    # print('WARNING: usage of map_plot is depreciated. This routine will be removed in future versions. Please use SingleMap instead.')

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

    if logplot:
        raise ValueError('Logplot option not supported yet!')

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

    if return_plot_handler:
        return M
    else:
        return M.figure


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
