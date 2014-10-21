# -*- coding: utf-8 -*-
from __future__ import division
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Variogram modelling
"""

import numpy as np
cimport numpy as np

from scipy.optimize import minimize

ctypedef np.double_t DTYPE_t  # double type for numpy arrays

cdef class Variogram(object):

    def __init__(self, **kwargs):
        pass

    def _orthodrome(self, double lon1_deg, double lat1_deg, double lon2_deg, double lat2_deg, double radius=6371000.):
        """
        calculate the orthodrome between two points with coordinates
        given in radians
        http://en.wikipedia.org/wiki/Great-circle_distance
        see also CDO code in file Gridcell.c

        @todo: how to deal with latitudes across the dateline ?

        Note that the same routine is also implemented for Data object at the moment!!

        Parameters
        ----------
        lon/lat : float
            coordinates of two points [degree]
        radius : float
            Earth radius (sphere) in [m]
        """

        cdef double lat1, lon1, lat2, lon2

        lat1 = np.deg2rad(lat1_deg)
        lon1 = np.deg2rad(lon1_deg)
        lat2 = np.deg2rad(lat2_deg)
        lon2 = np.deg2rad(lon2_deg)

        return np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1)
                         * np.cos(lat2) * np.cos(np.abs(lon2 - lon1))) * radius

    def _orthodrome_arr(self, double lon1_deg, double lat1_deg, np.ndarray[DTYPE_t, ndim=1] lon2_deg, np.ndarray[DTYPE_t, ndim=1] lat2_deg, double radius=6371000.):
        """
        ARRAY VERSION

        calculate the orthodrome between two points with coordinates
        given in radians
        http://en.wikipedia.org/wiki/Great-circle_distance
        see also CDO code in file Gridcell.c

        @todo: how to deal with latitudes across the dateline ?

        Note that the same routine is also implemented for Data object at the moment!!

        Parameters
        ----------
        lon/lat : float
            coordinates of two points [degree]
        radius : float
            Earth radius (sphere) in [m]
        """

        cdef double lat1
        cdef double lon1
        cdef np.ndarray[DTYPE_t, ndim=1] lat2
        cdef np.ndarray[DTYPE_t, ndim=1] lon2
        cdef deg2rad = np.pi / 180.

        lat1 = lat1_deg * deg2rad
        lon1 = lon1_deg * deg2rad
        lat2 = lat2_deg * deg2rad
        lon2 = lon2_deg * deg2rad

        return np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1)
                          * np.cos(lat2) * np.cos(np.abs(lon2 - lon1))) * radius


    #~ def _paired_distance(self, lon, lat, radius=6371.):
        #~ """
        #~ calculate paired distance between data points [km]
        #~ returns a matrix with squared distances [km]
        #~ """
#~
        #~ assert len(lon) == len(lat)
        #~ N = len(lon)
#~
        #~ pd = np.zeros((N,N))*np.nan
#~
        #~ for i in xrange(N):
            #~ for j in xrange(i+1,N):
                #~ d = self._orthodrome(lon[i], lat[i], lon[j], lat[j], radius=radius)
                #~ pd[i,j] = d
                #~ pd[j,i] = d
        #~ for i in xrange(N):
            #~ pd[i,i] = 0.
#~
        #~ return pd


    def _semivariance(self, np.ndarray[DTYPE_t, ndim=1] x, np.ndarray[DTYPE_t, ndim=1] lon, np.ndarray[DTYPE_t, ndim=1] lat, double h_km, double dh_km, double radius=6371.):
        """
        calculate semivariogram for a single lag

        Parameters
        ----------
        h_km : float
            distance lag [km]
        dh_km : float
            buffer zone for distance lag h [km]
        radius : float
            earth radius [km]
        """

        cdef int N
        cdef int i
        cdef int j
        cdef np.ndarray[DTYPE_t, ndim=1] d
        cdef double zval
        cdef int zcnt

        assert (x.ndim == 1)

        N = len(x)

        # calculate semivariance
        zcnt = 0
        zval = 0.
        for i in xrange(N):
            if i % 1000 == 0:
                print '    Variogramm calculation:', i, N

            # calculate distances for whole array to be fast
            d = self._orthodrome_arr(lon[i], lat[i], lon[i+1:], lat[i+1:], radius=radius)
            md = (d>=h_km-dh_km) & (d <= h_km+dh_km)
            Z = (x[i+1:] - x[i])**2.
            zval += np.sum(Z[md])
            zcnt += np.sum(md)

        if zcnt > 0:
            return 0.5 * zval / float(zcnt)
        else:
            return np.nan


    def semivariogram(self, x, lon, lat, lags, dlag):
        """
        calculate semivariogram for different lags

        Returns
        -------
        lags : ndarray
            vector fo lags
        gamma : ndarray
            and corresponding semivariance

        Parameters
        ----------
        x : ndarray
            array with data values
        lags : ndarray, list
            array with lag values
        lon : ndarray
            longitude coordinates
        lat : ndarray
            latitude coordinates
        lags : ndarray
            lags [km]
        dlag : float
            distance of lags [km] to specifiy valid buffer zone
        """

        cdef int i

        assert (lon.shape == lat.shape)
        assert(lon.shape == x.shape)

        gamma = np.ones(len(lags)) * np.nan
        for i in xrange(len(lags)):
            gamma[i] = self._semivariance(x, lon, lat, float(lags[i]), dlag)
        return lags, gamma

