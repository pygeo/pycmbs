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

from variogram_base import Variogram

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

        self.dlag = np.diff(h)[0]

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




