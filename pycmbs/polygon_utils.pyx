# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

from __future__ import division

STUFF = "Hello world"  # this is done to avoid this problem: http://stackoverflow.com/questions/8024805/cython-compiled-c-extension-importerror-dynamic-module-does-not-define-init-fu



#http://docs.cython.org/src/tutorial/numpy.html#efficient-indexing
import numpy as np
import copy
try:
    from osgeo import ogr
except:
    print 'WARNING: import of OSGEO did not work. Could cause trouble in usage of polygon funcitonalities!'

cimport numpy as np

ctypedef np.double_t DTYPE_t  # double type for numpy arrays







cdef class Polygon(object):
    """
    define a polygon
    """
    # in cython classes the class attributes need to be specified explicitely!
    #http://docs.cython.org/src/userguide/extension_types.html

    # Note however that the attributes can obviously not be directly accessed.
    # functions to return the attributes are therefore explicitely needed!

    cdef int id
    cdef float value
    cdef list poly
    cdef list _ogr_poly

    cdef bint ensure_positive

    def __init__(self, int id, list coordinates, ensure_positive=False):
        """
        Parameters
        ----------
        id : int
            unique identifier
        coordinates : list of tuples
            list of tuples (x,y) defining the polygon
        ensure_positive : bool
            if True, coordinates are shifted such that the longitude
            ranges from 0...360 degrees
        """
        self.poly = coordinates
        self.id = id
        self.value = np.nan
        self._ogr_poly = None

        if ensure_positive:
            self._shift_coordinates()

    property value:
        def __get__(self):
          return self.value
        def __set__(self, float value):
          self.value = value

    property id:
        def __get__(self):
          return self.id
        def __set__(self, int value):
          self.id = value

    property poly:
        def __get__(self):
          return self.poly
        def __set__(self, int value):
          self.poly = value

    property _ogr_poly:
        def __get__(self):
          return self._ogr_poly
        def __set__(self, int value):
          self._ogr_poly = value


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



    def point_in_poly_latlon(self, double lon, double lat, slow_warning=True):

        """
        This routine enables the calculation of the point-in-polygon
        problem for geographical coordinates. In particular it ensures
        that it is also working for polygons across the dateline.

        This is done, by shifting the polygon first in a way
        that no number across the dateline occur.

        It is assumed that longitudes go from -180. ... 180.

        Parameters
        ----------
        lon : float
            longitude [deg]
        lat : float
            latitude [deg]
        """

        if slow_warning:
            raise ValueError('This is far too slow for pixel wise processing!!!')

        if self._xmin() <= -180.:
            raise ValueError('minimum longitudes need to be >= -180.: ' + str(self._xmin()))
        if self._xmax() >= 180.:
            raise ValueError('maximum longitudes need to be <= 180.: ' + str(self._xmax()))

        # shift coordinates to ensure that only positive coordinates occur
        if lon < 0.:
            lon1 = lon + 360.
        else:
            lon1 = lon*1.

        tmp = Polygon(self.id, copy.deepcopy(self.poly), ensure_positive=True)  # ensure positive numbers for coordinates!
        res = tmp.point_in_poly(lon1, lat)
        del tmp
        return res

    def _get_point_count(self):
        return len(self.poly)


    def convertToOGRPolygon(self, ensure_positive=False):
        """
        convert polygon to an OGR polygon and return this
        http://pcjericks.github.io/py-gdalogr-cookbook/geometry.html#calculate-intersection-between-two-geometries

        Parameters
        ----------
        ensure_positive : bool
            if True, then the longitudes are shifted to positive numbers
            and wrapped around the dateline if necessary, resulting in coordinates
            from 0 ... 360 instead of -180 ... 180
        """

        if self._get_point_count() == 0:
            return None

        ring = ogr.Geometry(ogr.wkbLinearRing)
        for p in self.poly:
            x = p[0]*1.
            if ensure_positive:
                if x < 0.:
                    x += 360.
            ring.AddPoint(x, p[1])

        if not self.is_closed():
            # close polygon if needed
            x = self.poly[0][0]*1.
            if ensure_positive:
                if x < 0.:
                    x += 360.
            ring.AddPoint(x, self.poly[0][1])
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        self._ogr_poly = [poly]

        return poly

    def is_closed(self):
        """
        check if polygon is closed
        """
        p1 = self.poly[0]
        p2 = self.poly[self._get_point_count()-1]
        if (p1[0] == p2[0]) and (p1[1] == p2[1]):
            return True
        else:
            return False

    def _shift_coordinates(self):
        """
        shift coordinates to ensure values for longitude
        between 0...360 degrees
        """
        opoly = []
        for p in self.poly:
            if p[0] < 0.:
                x = p[0] + 360.
            else:
                x = p[0]
            opoly.append((x, p[1]))
        self.poly = opoly

    def xxx_not_so_fast_point_in_poly(self, double x, double y):
        """
        solve the point in area problem using OGR.
        in case that the point is within the polygon
        a valid geometry is returned. Otherwise the
        geometry object is empty

        Parameters
        ----------
        x : float
            x-coordinate = longitude
        y : float
            y-coordinate = latitude
        """
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)

        if self._ogr_poly is None:
            poly = self.convertToOGRPolygon()  # todo this could be done more efficiently by storing the converted polygon once and then simply use it.
        else:
            poly = self._ogr_poly[0]

        # now check for intersection
#~         print poly.ExportToWkt()
#~         print point.ExportToWkt()
        intersection = poly.Intersection(point)
#~         print intersection.ExportToWkt()
        N = intersection.GetPointCount()
        if N == 1:
            return True
        elif N == 0:
            return False
        else:
            print N
            raise ValueError('Some invalid number of points! Should be one or zero!')


    def point_in_poly(self, double x, double y):
        """
        Parameters
        ----------
        x : float
            x-coordinate of the point to be investigated
        y : float
            y-coordinate of the point to be investigated
        """

        cdef int i
        cdef double p1x, p1y, p2x, p2y

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


def fast_point_in_poly(np.ndarray[DTYPE_t, ndim=2] lon, np.ndarray[DTYPE_t, ndim=2] lat, Polygon P):

    cdef int i, j
    cdef int ny = lon.shape[0]
    cdef int nx = lon.shape[1]

    cdef id

    cdef np.ndarray[DTYPE_t, ndim=2] mask = np.zeros([ny, nx]) * np.nan
    DTYPE = np.double

    assert lon.dtype == DTYPE
    assert lat.dtype == DTYPE

    id = float(P.id)

    for i in range(ny):
        if i % 10 == 0:
            print 'Processing line: ', i, ny
        for j in xrange(nx):
            if P.point_in_poly(lon[i, j], lat[i, j]):
                if np.isnan(mask[i, j]):
                    mask[i, j] = id
                else:
                    print i, j, lon[i, j], lat[i, j]
                    raise ValueError('Overlapping polygons not supported yet!')
            else:
                pass

    return mask



def get_point_in_poly_mask(np.ndarray[DTYPE_t, ndim=2] mask, np.ndarray[DTYPE_t, ndim=2] lon, np.ndarray[DTYPE_t, ndim=2] lat, Polygon P):
    """
    routine to calculte point in polygon problem in a fast way using cython

    note that no checks for geometry consistenc<y are performed

    Parameters
    ----------
    mask : ndarray
    lon : ndarray
    lat : ndarray
    P : Polygon
    id : valu to assign
    """

    cdef int i, j
    cdef double id

    id = float(P.id)

    ny, nx = np.shape(lon)
    for i in xrange(ny):
        if i % 1000 == 0:
            print 'Rasterization complete by ', np.round(100. * float(i) / float(ny),0), '%               \r',
        for j in xrange(nx):
            if P.point_in_poly(lon[i, j], lat[i, j]):
                if np.isnan(mask[i, j]):
                    mask[i, j] = id
                else:
                    print ''
                    print ''
                    print i, j, lon[i, j], lat[i, j], mask[i,j]
                    print ''
                    print ''
                    raise ValueError('Overlapping polygons not supported yet!')
            else:
                pass
