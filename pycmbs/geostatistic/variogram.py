# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Variogram modelling
"""

import numpy as np
from scipy.optimize import minimize

class Variogram(object):

    def __init__(self, **kwargs):
        pass

    def _orthodrome(self, lon1_deg, lat1_deg, lon2_deg, lat2_deg, radius=6371000.):
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

        lat1 = np.deg2rad(lat1_deg)
        lon1 = np.deg2rad(lon1_deg)
        lat2 = np.deg2rad(lat2_deg)
        lon2 = np.deg2rad(lon2_deg)

        return np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1)
                         * np.cos(lat2) * np.cos(np.abs(lon2 - lon1))) * radius


    def _paired_distance(self, lon, lat, radius=6371.):
        """
        calculate paired distance between data points [km]
        returns a matrix with squared distances [km]
        """

        assert len(lon) == len(lat)
        N = len(lon)

        pd = np.zeros((N,N))*np.nan

        for i in xrange(N):
            for j in xrange(i+1,N):
                d = self._orthodrome(lon[i], lat[i], lon[j], lat[j], radius=radius)
                pd[i,j] = d
                pd[j,i] = d
        for i in xrange(N):
            pd[i,i] = 0.

        return pd


    def _semivariance(self, x, lon, lat, h_km, dh_km, radius=6371.):
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

        assert (x.ndim == 1)

        N = len(x)

        # calculate pairwise distance
        # TODO: this is calculating only Eucledian distance at the moment!
        # TODO replace this by proper calculation of orthodrome!

        #~ pd = self._paired_distance(lon, lat)
        #~ assert pd.shape[0] == N

        # calculate semivariance
        Z = list()
        for i in xrange(N):  # TODO: do this more efficient (e.g. only looking for points which are within distance anyway)
            print 'Variogramm calculation:', i, N
            for j in xrange(i+1,N):
                # calculate distance between points
                d = self._orthodrome(lon[i], lat[i], lon[j], lat[j], radius=radius)
                if (d >= h_km-dh_km) and (d <= h_km+dh_km):
                    Z.append((x[i]-x[j])**2.)
        if len(Z) > 0:
            return np.sum(Z) / (2. * len(Z))
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
        assert (lon.shape == lat.shape)
        assert(lon.shape == x.shape)

        gamma = np.ones(len(lags)) * np.nan
        for i in xrange(len(lags)):
            gamma[i] = self._semivariance(x, lon, lat, float(lags[i]), dlag)
        return lags, gamma

class SphericalVariogram(Variogram):

    def __init__(self, **kwargs):
        super(SphericalVariogram, self).__init__(**kwargs)

    def _get_initial_parameters(self, sill=5., nugget=0., range=2.):
        return [nugget, sill, range]

    def fit(self, h, gamma):
        """
        fit theoretical model to empirical data

        Returns
        -------
        returns a dictionary with model parameters
        """

        self._h = np.asarray(h)*1.
        self._gamma = gamma

        x0 = self._get_initial_parameters()

        res = minimize(self.cost, x0, method='nelder-mead',
                   options={'xtol': 1e-8, 'disp': False})
        self.model_parameters = {'sill' : res.x[1], 'range' : res.x[2], 'nugget' : res.x[0]}
        return self.model_parameters

    def cost(self, x):
        nugget = x[0]
        sill = x[1]
        range = x[2]

        y = self.model(self._h, sill, nugget, range)

        return np.sum((y - self._gamma)**2.)

    def model(self, h, sill, nugget, range):
        """
        References
        ----------
        Roman et al. (2009): doi:10.1016/j.rse.2009.07.009
        """
        c = sill
        c0 = nugget
        a = range

        if np.any(h < 0.):
            raise ValueError('Distances are not allowed to be smaller than zero!')

        gamma = c0 + c * (1.5*h/a - 0.5*((h/a)**3.))
        gamma[h>a] = c0+c

        return gamma

    def plot(self, h, gamma):
        """
        plot semivariogram
        """
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.plot(h, gamma, 'x')
        ax.set_ylabel('$\gamma$')
        ax.set_xlabel('$lag distance [km]$')

        gmodel = self.model(h, self.model_parameters['sill'], self.model_parameters['nugget'], self.model_parameters['range'])
        ax.plot(h, gmodel, '-', color='red')




