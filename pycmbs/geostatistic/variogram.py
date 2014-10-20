"""
Variogram modelling
"""

import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform



class Variogram(object):

    def __init__(self, **kwargs):
        pass

    def _semivariance(self, x, lon, lat, h, dh):
        """
        calculate semivariogram for a single lag

        Parameters
        ----------
        h : float
            distance lag
        dh : float
            buffer zone for distance lag h
        """

        assert (x.ndim == 1)


        N = len(x)

        # calculate pairwise distance
        # TODO: this is calculating only Eucledian distance at the moment!
        # TODO replace this by proper calculation of orthodrome!

        pd = squareform( pdist( np.vstack([lon, lat]).T, 'eucledian' ) )
        assert pd.shape[0] == N

        # calculate semivariance
        Z = list()
        for i in xrange(N):  # TODO: do this more efficient (e.g. only looking for points which are within distance anyway)
            for j in xrange(i+1,N):
                if (pd[i,j] >= h-dh) and (pd[i,j] <= h+dh):
                    Z.append((x[i]-x[j])**2.)
        return np.sum(Z) / (2. * len(Z))


    def semivariogram(self, x, lon, lat, lags):
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
        """
        assert (lon.shape == lat.shape)
        assert(lon.shape == x.shape)

        gamma = np.ones(len(lags)) * np.nan
        for i in xrange(len(lags)):
            gamma[i] = self._semivariance(x, lon, lat, lags[i])
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




