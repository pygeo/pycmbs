# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from pycmbs.data import Data
import matplotlib.pyplot as plt
import numpy as np


class Geostatistic(object):
    def __init__(self, x, range_bins=None):
        assert isinstance(x, Data)
        if range_bins is None:
            raise ValueError('ERROR: you need to specifiy the range bins!')
        self.x = x
        self.range_bins = np.asarray(range_bins)
        self._check()

    def _check(self):
        if np.any(np.diff(self.range_bins)<0.):
            raise ValueError('Bins are not in ascending order!')
        if self.x.data.ndim != 2:
            raise ValueError('Currently only support for 2D data')

    def set_center_position(self, lon, lat):
        """
        set position of location to be used for analysis as center
        lon : float
            longitude of center position [deg]
        lat : float
            latitude of center position [deg]
        """
        self.lon_center = lon
        self.lat_center = lat

    def _get_center_pos(self):
        """
        get indices of position of center coordinate within grid
        """
        if not hasattr(self, 'lon_center'):
            raise ValueError('ERROR: You need to specify first the center position!')
        d = np.abs((self.x.lon - self.lon_center)**2. + (self.x.lat - self.lat_center)**2.)
        dmin = d.min()
        m = d == dmin

        idx = np.indices(d.shape)
        i = idx[0][m][0]
        j = idx[1][m][0]

        if (np.abs(1.-self.x.lon[i,j]/self.lon_center) > 0.05) or (np.abs(1.-self.x.lat[i,j]/self.lat_center) > 0.05):  # at least 5% acc.
            i = None
            j = None
        return i, j


    def plot_semivariogram(self, ax=None, color='red', logy=False):
        """
        plot semivariogram
        """
        if not hasattr(self, '_distance'):
            self._calculate_distance()

        f, ax = self._get_figure_ax(ax)
        r, sigma = self.calc_semivariance()

        if logy:
            ax.semilogy(r, sigma, 'x', color=color)
        else:
            ax.plot(r, sigma, 'x', color=color)
        ax.set_ylabel('$\sigma^2$ / 2 (isotropic)')
        ax.set_xlabel('distance [km]')
        ax.grid()
        return ax

    def plot_percentiles(self, p, ax=None, logy=False):
        """
        plot percentiles

        p : list
            list of percentiles [0 ... 1]
        """
        f, ax = self._get_figure_ax(ax)
        for e in p:
            r, v = self.calc_percentile(e)
            if logy:
                ax.semilogy(r, v, 'x-', label='p='+str(e))
            else:
                ax.plot(r, v, 'x-', label='p='+str(e))
        ax.set_ylabel(self.x._get_unit())
        ax.set_xlabel('distance [km]')
        ax.grid()
        ax.legend(loc='upper left', prop={'size':10}, ncol=2)
        return ax

    def calc_percentile(self, p):
        """
        p [0 ... 1]
        """
        bounds = self.range_bins
        r = []
        v = []
        for b in bounds:
            d = self._get_data_distance(0., b)
            if len(d) < 1:
                continue
            r.append(b)
            v.append(np.percentile(d, p*100.))  # percentile value
        return np.asarray(r), np.asarray(v)




    def _get_data_distance(self, lb, ub):
        """
        get data for distance where

        lb<= r < ub
        """
        if not hasattr(self, '_distance'):
            self._calculate_distance()
        m = (self._distance >= lb) & (self._distance < ub)
        return self.x.data[m].flatten()


    def calc_semivariance(self):
        """
        calculate semivariance for selected range bins
        """
        bounds = self.range_bins
        r = []
        v = []
        for b in bounds:
            d = self._get_data_distance(0., b)
            r.append(b)
            v.append(0.5*d.var())  #semivariance
        return np.asarray(r), np.asarray(v)

    def _get_figure_ax(self, ax):
        if ax is None:
            f = plt.figure()
            ax = f.add_subplot(111)
        else:
            f = ax.figure
        return f, ax


    def _calculate_distance(self):
        """
        calculate distance [km]
        """
        if not hasattr(self, 'lon_center'):
            raise ValueError('ERROR: You need to specify first the center position!')
        self._distance = self.x.distance(self.lon_center, self.lat_center) / 1000.
