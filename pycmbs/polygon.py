# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import numpy as np
from pycmbs.polygon_utils import Polygon
from polygon_utils import fast_point_in_poly


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
            faster: same as full, but fully implemented in cython. However
            this is not necessarily faster than 'full'
        """

        if not hasattr(self, 'mask'):
            raise ValueError('Output array not existing yet!')
        if not isinstance(P, Polygon):
            raise ValueError('No Polygon object provided!')
        if P.id < 0:
            raise ValueError('ERROR: ID value must not be negative!')

        if method != 'full':
            raise ValueError('Only FULL method currently thoroughly validated!')

        # check if data across the dateline. This is currently not supported yet!
        if method != 'full':
            if P.across_dateline():  # check routine not existing yet!
                print 'WARNING: seems that polygon is across the dateline. This is currently not supported yet! SKIPPING polygon with ID: ', P.id
                return

        id = float(P.id)

        # check if id already existing
        msk = self.mask == float(id)
        if msk.sum() > 0:
            raise ValueError('The ID value is already existing!')

        if method == 'full':
            ny, nx = self.lon.shape
            for i in xrange(ny):
                #~ if i % 10 == 0:
                print 'Rasterization ... ', 100. * float(i) / float(ny), '%'
                for j in xrange(nx):
                    if P.point_in_poly(self.lon[i, j], self.lat[i, j]):
                        if np.isnan(self.mask[i, j]):
                            self.mask[i, j] = id
                        else:
                            print i, j, self.lon[i, j], self.lat[i, j]
                            raise ValueError('Overlapping polygons not supported yet!')
                    else:
                        pass
        elif method == 'faster':  # an alternative implementation. This is however not necessarily faster than 'full'
            print 'Using CYTHON method for rasterization!'
            self.mask = fast_point_in_poly(self.lon, self.lat, P)
        elif method == 'fast':
            xmin, xmax, ymin, ymax = P.bbox()
            # determine bounding box very roughly NOT THAT THIS MIGHT BE NOT CORRECT
            valid_points = (self.lon >= xmin) & (self.lon <= xmax) & (self.lat >= ymin) & (self.lat <= ymax)
            plon = self.lon[valid_points]  # preselect points that are likely to fall within polygon
            plat = self.lat[valid_points]

            resmsk = np.zeros_like(plon) * np.nan
            for i in xrange(len(plon)):
                if i % 500 == 0:
                    print i, len(plon)
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
