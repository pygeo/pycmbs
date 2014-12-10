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
from scipy import stats

from variogram_base import Variogram
import matplotlib.pyplot as plt
import pickle


class SphericalVariogram(Variogram):

    def __init__(self, **kwargs):
        super(SphericalVariogram, self).__init__(**kwargs)

    def _get_initial_parameters(self, sill=5., nugget=0., range=2.):
        return [nugget, sill, range]

    def fit(self, h, gamma):
        """
        fit theoretical model to empirical data

        Parameters
        ----------
        h : ndarray
            ranges
        gamma : ndarray
            experimental semivariance

        Returns
        -------
        returns a dictionary with model parameters
        """

        self.dlag = np.diff(h)[0]

        msk = ~np.isnan(gamma)

        self._h = np.asarray(h[msk]) * 1.
        self._gamma = gamma[msk]

        if len(self._gamma) < 1:
            self.model_parameters = {'sill': np.nan, 'range': np.nan, 'nugget': np.nan, 'fit_success': False, 'r_value' : np.nan}
            return self.model_parameters

        x0 = self._get_initial_parameters(sill=self._gamma.max(), range=self._h[np.argmax(self._gamma)])

        res = minimize(self.cost, x0, method='nelder-mead',
                       options={'xtol': 1e-8, 'disp': True})

        self.model_parameters = {'sill': res.x[1], 'range': res.x[2], 'nugget': res.x[0], 'fit_success': res['success']}

        # calculate correlation parameters between model fit and experimental data
        yfit = self.model(self._h, self.model_parameters['sill'] , self.model_parameters['nugget'], self.model_parameters['range'] )
        slope, intercept, r_value, p_value, std_err = stats.linregress(self._gamma, yfit)
        self.model_parameters.update({'r_value' : r_value, 'slope' : slope, 'intercept' : intercept, 'p_value' : p_value})

        return self.model_parameters

    def cost(self, x):
        nugget = x[0]
        sill = x[1]
        range = x[2]

        y = self.model(self._h, sill, nugget, range)
        r = np.sum((y - self._gamma) ** 2.)
        if range < 0.:  # no negative range values!
            r += 10. ** 5.
        return r

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

        gamma = c0 + c * (1.5 * h / a - 0.5 * ((h / a) ** 3.))
        gamma[h > a] = c0 + c

        return gamma

    def plot(self, h, gamma, ax=None, color='red'):
        """
        plot semivariogram
        """
        if ax is None:
            f = plt.figure()
            ax = f.add_subplot(111)
        else:
            ax = ax

        # plot original data
        ax.plot(h, gamma, 'x', color=color)
        ax.set_ylabel('$\gamma$')
        ax.set_xlabel('$lag distance [km]$')

        if self.model_parameters['fit_success']:  # if variogram fitting was not sucessfull it is not fitted!
            lab = 'r=' + str(np.round(self.model_parameters['r_value'], 2))
            gmodel = self.model(h, self.model_parameters['sill'], self.model_parameters['nugget'], self.model_parameters['range'])
            ax.plot(h, gmodel, '-', color=color, label=lab)
            ax.legend(prop={'size' : 8})
