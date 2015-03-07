# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

"""
MODULE for spatial grid manipulations
"""

import numpy as np
import matplotlib.delaunay as triang
from matplotlib import pylab as pl
import matplotlib as mpl

import matplotlib.patches as mpatches
import matplotlib.path as mpath
from matplotlib.collections import PatchCollection


class Grid(object):

    """
    general class for spatial grids
    """

    def __init__(self, lat_rad, lon_rad, sphere_radius=None,
                 gridtype='regular'):
        """
        Parameters
        ----------

        lat_rad : float
            latitude [rad]
        lon_rad : float
            longitude [rad]
        sphere_radius : float
            Radius of sphere [m]
        gridtype : str
            identifier to specify the type of grid to be used ['regular']
        """
        self._lat0 = lat_rad
        self._lon0 = lon_rad
        self.lat = self._lat0.flatten()
        self.lon = self._lon0.flatten()
        self.gridtype = gridtype

        if sphere_radius is None:
            raise ValueError('Grid needs to be initialized \
                               with sphere radius!')
        self.radius = sphere_radius
        self._interpolated = False

    def orthodrome(self, lon1, lat1, lon2, lat2):
        """
        calculate the orthodrome between two points with coordinates
        given in radians
        http://en.wikipedia.org/wiki/Great-circle_distance
        see also CDO code in file Gridcell.c

        @todo: how to deal with latitudes across the dateline ?
        """
        return np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1)
                         * np.cos(lat2) * np.cos(np.abs(lon2 - lon1))) * self.radius

    def calc_cell_area(self):
        """
        calculate cell area of the grid cell using spherical triangles
        within the CDOs, the area caluculation is implemented in grid.c
        """
        raise ValueError('Calculation of cell area not implemented so far')

    def _delaunay_triangulation(self):
        self.cens, self.edg, self.tri, self.neig = triang.delaunay(
            self.lon, self.lat)
        self._interpolated = True

    def plot(self, ax=None):
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111)
        ax.plot(self.lon, self.lat, 'x')

    def draw_edges(self, ax=None):
        if not self._interpolated:
            self._delaunay_triangulation()
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111)

        for e in self.edg:
            ax.plot(self.lon[e], self.lat[e], '--', color='grey')

    def _calc_triangles(self, ax=None):
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111)
        if self.gridtype != 'regular':
            raise ValueError('only supported for regular grids so far!')
        if self._lon0.ndim != 2:
            raise ValueError('2D grid required!')

        ny, nx = self._lon0.shape
        lons = []
        lats = []
        for i in range(1, ny - 1):   # todo revise this calculation method!!!
            for j in range(1, nx - 1):
                # half way to corner neighbors?

                # upper left
                lons.append(self._lon0[i - 1, j - 1] + 0.5
                            * (self._lon0[i, j] - self._lon0[i - 1, j - 1]))
                lats.append(self._lat0[i - 1, j - 1] + 0.5
                            * (self._lat0[i, j] - self._lat0[i - 1, j - 1]))

                # lower right
                lons.append(self._lon0[i, j] + 0.5
                            * (self._lon0[i + 1, j + 1] - self._lon0[i, j]))
                lats.append(self._lat0[i, j] + 0.5
                            * (self._lat0[i + 1, j + 1] - self._lat0[i, j]))

                # lower left
                # lons.append(self._lon0[i-1,j]+0.5
                #*(self._lon0[i+1,j+1]-self._lon0[i,j]))
                # add also current center coordinate for
                # plotting vizualization
                lons.append(self._lon0[i, j])
                lats.append(self._lat0[i, j])
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        cens, edg, tri, neig = triang.delaunay(lons, lats)
        self._do_delaunay_plot(tri, lons, lats, ax)

    def plot_voronoi(self, ax=None, show_latlon=True):
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111)

        # calculate half size position of each edge and
        # do delaunay interpolation again. This
        # gives then the Voronoi diagram
        if not self._interpolated:
            self._delaunay_triangulation()

        lons = []
        lats = []
        for e in self.edg:
            # half way position
            plon = self.lon[e[0]] + 0.5 * (self.lon[e[1]] - self.lon[e[0]])
            plat = self.lat[e[0]] + 0.5 * (self.lat[e[1]] - self.lat[e[0]])
            lons.append(plon)
            lats.append(plat)
        lons = np.asarray(lons)
        lats = np.asarray(lats)

        cens, edg, tris, neig = triang.delaunay(lons, lats)
        self._do_delaunay_plot(tris, lons, lats, ax)
        ax.plot(lons, lats, marker='x', linestyle='None', color='green')

        if show_latlon:
            ax.plot(self.lon, self.lat, 'o', color='black')

    def plot_delaunay_grid(self, ax=None, show_latlon=True):
        if ax is None:
            fig = pl.figure()
            ax = fig.add_subplot(111)
        if not self._interpolated:
            self._delaunay_triangulation()

        #--- do actual plotting ---
        self._do_delaunay_plot(self.tri, self.lon, self.lat, ax)

        if show_latlon:
            ax.plot(self.lon, self.lat, 'o', color='black')

    def _do_delaunay_plot(self, tri, lon, lat, ax):
        patches = []
        Path = mpath.Path
        for t in tri:
            # t[0], t[1], t[2] are the points indexes of the triangle
            t_i = [t[0], t[1], t[2]]
            #... vertices for a single triangle
            verts = np.asarray([lon[t_i], lat[t_i]]).T

            #--- specify how vertices are interconnected
            #    (here simple connection by lines)
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO]

            #--- construct object and append to library of objects ---
            path = mpath.Path(verts, codes, closed=True)
            patches.append(mpatches.PathPatch(path))

        cmap = pl.cm.get_cmap('jet', 50)
        norm = mpl.colors.Normalize(vmin=None, vmax=None)  # colorbar mapping
        # construct library of all objects
        self._collection = PatchCollection(patches, cmap=cmap, norm=norm,
                                           alpha=1., match_original=False)
        colors = 100 * np.random.rand(len(patches))
        # assign data values here
        self._collection.set_array(np.array(colors))

        im = ax.add_collection(self._collection)
        ax.set_xlim(lon.min(), lon.max())
        ax.set_ylim(lat.min(), lat.max())
