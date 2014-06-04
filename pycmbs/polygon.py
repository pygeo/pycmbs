# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import numpy as np


class Polygon(object):
    """
    define a polygon
    """
    def __init__(self, id, coordinates):
        """
        Parameters
        ----------
        id : int
            unique identifier
        coordinates : list of tuples
            list of tuples (x,y) defining the polygon
        """
        self.poly = coordinates
        self.id = id

    def _xcoords(self):
        return np.asarray([t[0] for t in self.poly])

    def _ycoords(self):
        return np.asarray([t[1] for t in self.poly])

    def _xmin(self):
        return self._xcoords().min()

    def _xmax(self):
        return self._xcoords().max()

    def _ymin(self):
        return self._ycoords().min()

    def _ymax(self):
        return self._ycoords().max()

    def bbox(self):
        """
        returns bbox for rasterization of data

        Todo
        ----
        how to handle bbox across coordinate borders ???
        """
        return [self._xmin(), self._xmax(), self._ymin(), self._ymax()]

    def point_in_poly(self, x, y):
        """
        Parameters
        ----------
        x : float
            x-coordinate of the point to be investigated
        y : float
            y-coordinate of the point to be investigated

        TODO
        ----
        does that work also across the deadline and datum line?
        """

        n = len(self.poly)
        inside = False

        p1x, p1y = self.poly[0]
        for i in xrange(n + 1):
            p2x, p2y = self.poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside


class Raster(object):

    def __init__(self, lon, lat):
        """
        Parameters
        ----------
        lat : ndarray
            2D latitude array
        lon : ndarray
            2D longitude array

        Todo
        ----
        TODO support different geodetic systems
        """

        self.lon = lon
        self.lat = lat

        self._check()

    def _check(self):
        if self.lon.ndim != 2:
            raise ValueError('Longitude need to be 2D')
        if self.lat.ndim != 2:
            raise ValueError('Latitude need to be 2D')
        if self.lat.shape != self.lon.shape:
            raise ValueError('ERROR: Inconsistent shapes')

    def rasterize_polygons(self, polygons, method='full'):
        """
        rasterize polygons and return a rasterdataset
        with the same geometry

        Parameters
        ----------
        polygons : list
            list of Polygon objects that need to be rasterized
        method : str
            full: itterate over entire spatial domain (most precise)
            fast: calculate first bbox and then itterate (faster)
        """
        self.mask = np.zeros(self.lon.shape) * np.nan
        for P in polygons:
            self._rasterize_single_polygon(P, method=method)
        self.mask = np.ma.array(self.mask, mask=np.isnan(self.mask))

    def _rasterize_single_polygon(self, P, method='full'):
        """
        rasterize a single polygon

        requires that the output array has already been initialized

        Parameters
        ----------
        P : Polygon
            Polygon object
        id : int
            identifier
        method : str
            full: itterate over entire spatial domain (most precise)
            fast: calculate first bbox and then itterate (faster)
        """

        if not hasattr(self, 'mask'):
            raise ValueError('Output array not existing yet!')
        if not isinstance(P, Polygon):
            raise ValueError('No Polygon object provided!')
        if P.id is None:
            raise ValueError('ERROR: ID value must not be None')

        id = float(P.id)

        # check if id already existing
        msk = self.mask == float(id)
        if msk.sum() > 0:
            raise ValueError('The ID value is already existing!')

        if method == 'full':
            ny, nx = self.lon.shape
            for i in xrange(ny):
                for j in xrange(nx):
                    if P.point_in_poly(self.lon[i, j], self.lat[i, j]):
                        if np.isnan(self.mask[i, j]):
                            self.mask[i, j] = id
                        else:
                            print i, j, self.lon[i, j], self.lat[i, j]
                            raise ValueError('Overlapping polygons not supported yet!')
                    else:
                        pass

        elif method == 'fast':
            xmin, xmax, ymin, ymax = P.bbox()
            valid_points = (self.lon >= xmin) & (self.lon <= xmax) & (self.lat >= ymin) & (self.lat <= ymax)
            plon = self.lon[valid_points]  # preselect points that are likely to fall within polygon
            plat = self.lat[valid_points]

            resmsk = np.zeros_like(plon) * np.nan
            for i in xrange(len(plon)):
                if P.point_in_poly(plon[i], plat[i]):
                    if np.isnan(resmsk[i]):
                        resmsk[i] = id
                    else:
                        raise ValueError('Overlapping polygons not supported yet!')
            # reassign mask
            newmsk = np.ones_like(self.lon) * np.nan
            newmsk[valid_points] = resmsk
            newvalues = ~np.isnan(newmsk)
            if np.any(~np.isnan(self.mask[newvalues])):
                print newmsk[newvalues]
                print sum(~np.isnan(self.mask[newvalues]))
                raise ValueError('Some values had been set already before!')
            self.mask[newvalues] = newmsk[newvalues]

        else:
            raise ValueError('ERROR: Invalid method')
