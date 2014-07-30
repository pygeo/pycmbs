from __future__ import division
#http://docs.cython.org/src/tutorial/numpy.html#efficient-indexing
import numpy as np
cimport numpy as np

#~ ctypedef np.double_t DTYPE_t



cdef class Polygon(object):
    """
    define a polygon
    """
    # in cython classes the class attributes need to be specified explicitely!
    #http://docs.cython.org/src/userguide/extension_types.html

    # Note however that the attributes can obviously not be directly accessed.
    # functions to return the attributes are therefore explicitely needed!

    cdef int id
    cdef list poly

    def __init__(self, int id, list coordinates):
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

    def point_in_poly(self, double x, double y):
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



def fast_point_in_poly(np.ndarray lon, np.ndarray lat, Polygon P):

    cdef int i, j
    cdef int ny = lon.shape[0]
    cdef int nx = lon.shape[1]

    cdef np.ndarray mask = np.zeros([ny, nx]) * np.nan
    DTYPE = np.double

    assert lon.dtype == DTYPE
    assert lat.dtype == DTYPE

    for i in range(ny):
        if i % 10 == 0:
            print 'Processing line: ', i, ny
        for j in range(nx):
            if P.point_in_poly(lon[i, j], lat[i, j]):
                if np.isnan(mask[i, j]):
                    mask[i, j] = float(P.id)
                else:
                    print i, j, lon[i, j], lat[i, j]
                    raise ValueError('Overlapping polygons not supported yet!')
            else:
                pass


#~ # look again more into efficient indexing in CYTHON documentation!!!


    return mask
