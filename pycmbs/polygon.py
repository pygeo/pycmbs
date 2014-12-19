# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import numpy as np
from pycmbs.polygon_utils import Polygon
from pycmbs.polygon_utils import get_point_in_poly_mask

import multiprocessing
import itertools

THE_POLYGONS=None
THE_MASK = None


def argument_mapper_rasterize(x):

    _rasterize_polygon(x[0], x[1], x[2], x[3])

def _rasterize_polygon(lon, lat, i, method):
    """
    routine to rasterize polygon
    as separate function to allow for parallel operations
    """
    global THE_POLYGONS
    global THE_MASK

    P = THE_POLYGONS[i]

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
    msk = THE_MASK == float(id)
    if msk.sum() > 0:
        raise ValueError('The ID value is already existing!')

    if method == 'full':
        #~ hlp = np.ones_like(THE_MASK)*np.nan

        #~ uvals = np.unique(THE_MASK[~np.isnan(THE_MASK)])
        #~ print 'THE_MASK before: ', i, P.id, len(uvals)
        get_point_in_poly_mask(THE_MASK, lon, lat, P)  # use fast cython based method
        #~ print 'THE_MASK after: ', hlp

        #~ mm = ~np.isnan(hlp)
        #~ print 'N-valid: ', np.sum(mm)
        #~ THE_MASK[mm] = hlp[mm]*1.
        #~ uvals = np.unique(THE_MASK[~np.isnan(THE_MASK)])
        #~ print 'THEMASK-after: ', i, len(uvals)
        #~ del hlpe

    elif method == 'fullold':
        assert False
        ny, nx = lon.shape
        for i in xrange(ny):
            #~ if i % 1 == 0:
            print 'Rasterization complete by ', np.round(100. * float(i) / float(ny),0), '%               \r',
            for j in xrange(nx):
                if P.point_in_poly(lon[i, j], lat[i, j]):
                    if np.isnan(mask[i, j]):
                        mask[i, j] = id
                    else:
                        print i, j, lon[i, j], lat[i, j]
                        raise ValueError('Overlapping polygons not supported yet!')
                else:
                    pass

    elif method == 'faster':  # an alternative implementation. This is however not necessarily faster than 'full'
        assert False, 'Option currently not supported!'
        print 'Using CYTHON method for rasterization!'
        mask = fast_point_in_poly(lon, lat, P)
    elif method == 'fast':
        assert False
        xmin, xmax, ymin, ymax = P.bbox()
        # determine bounding box very roughly NOT THAT THIS MIGHT BE NOT CORRECT
        valid_points = (lon >= xmin) & (lon <= xmax) & (lat >= ymin) & (lat <= ymax)
        plon = lon[valid_points]  # preselect points that are likely to fall within polygon
        plat = lat[valid_points]

        resmsk = np.zeros_like(plon) * np.nan
        for i in xrange(len(plon)):
            #~ if i % 500 == 0:
                #~ print i, len(plon)
            if P.point_in_poly(plon[i], plat[i]):
                if np.isnan(resmsk[i]):
                    resmsk[i] = id
                else:
                    raise ValueError('Overlapping polygons not supported yet!')
        # reassign mask
        newmsk = np.ones_like(lon) * np.nan
        newmsk[valid_points] = resmsk
        newvalues = ~np.isnan(newmsk)
        if np.any(~np.isnan(mask[newvalues])):
            print newmsk[newvalues]
            print sum(~np.isnan(mask[newvalues]))
            raise ValueError('Some values had been set already before!')
        mask[newvalues] = newmsk[newvalues]

    else:
        raise ValueError('ERROR: Invalid method')




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

    def rasterize_polygons(self, polygons, method='full', nproc=1):
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
        nproc : int
            number of processors to be used in parallel
        """

        global THE_POLYGONS
        global THE_MASK

        THE_POLYGONS = polygons

        nproc = 1  #no parallization as otherwise error occurs

        #~ self.mask = np.zeros(self.lon.shape) * np.nan
        THE_MASK = np.zeros(self.lon.shape) * np.nan

        if nproc == 1:
            #~ for P in polygons:
            for i in xrange(len(polygons)):
                self._rasterize_single_polygon(i, method=method)
        else:
            assert False, 'Parallelization does not work properly here yet'
            # problem is that the final array concatenation resuls in an invalid array
            # no idea why!
            # caution: this might not work in case of overlapping polygons

            N = len(polygons)
            pool = multiprocessing.Pool(processes=nproc)
            the_args = itertools.izip(itertools.repeat(self.lon,N), itertools.repeat(self.lat,N), xrange(len(polygons)), itertools.repeat(method,N))
            r = pool.map(argument_mapper_rasterize, the_args)
            pool.close()


        self.mask = np.ma.array(THE_MASK, mask=np.isnan(THE_MASK))

        #~ print self.mask
        #~ print  self.mask[~np.isnan(self.mask)]
        #~ print THE_MASK
        #~ print  THE_MASK[~np.isnan(THE_MASK)]

        #~ stop


    def _rasterize_single_polygon(self, i, method='full'):
        """
        rasterize a single polygon; this function should in general not
        be used alone, but only be called from rasterize_polygons !

        requires that the output array has already been initialized

        Parameters
        ----------
        i : int, Polygon
            if an integer is provided it is assumed that it points to
            a global list THE_POLYGONS

            if a Polygon is provided, THE_POLYGONS are set

        id : int
            identifier

        method : str
            full: itterate over entire spatial domain (most precise)
            fast: calculate first bbox and then itterate (faster)
            faster: same as full, but fully implemented in cython. However
            this is not necessarily faster than 'full'
        """

        if isinstance(i, Polygon):
            # a single polygon is provided
            global THE_POLYGONS
            global THE_MASK

            THE_POLYGONS = [i]
            THE_MASK = np.zeros(self.lon.shape) * np.nan
            _rasterize_polygon(self.lon, self.lat, 0, method)
            self.mask = np.ma.array(THE_MASK, mask=np.isnan(THE_MASK))
        elif isinstance(i, int):
            # an index is provided
            _rasterize_polygon(self.lon, self.lat, i, method)
        else:
            raise ValueError('Only integers or Polygon objects are allowed as arguments')


