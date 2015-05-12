# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

from pycmbs.data import Data
from variogram import SphericalVariogram

import matplotlib.pyplot as plt
import numpy as np


class Geostatistic(object):
    def __init__(self, x, lags=None, maxdist=None):
        """
        Geostatistical calculations on a data object

        Parameters
        ----------
        x : Data
            data to be analyzed
        lags : list
            list of bins to perform analysis
        maxdist : float
            maximum distance for analyis [km]
        """
        assert isinstance(x, Data)
        if lags is None:
            raise ValueError('ERROR: you need to specifiy the range bins!')
        self.x = x
        self.lags = np.asarray(lags)
        self._check()
        self.statistic = {}

        self.maxdist = maxdist
        self._distfiltered = False

    def _filter_data_maxdist(self):
        """
        filters the data obect to contain only values within a maximum distance
        This allows more efficient calculations afterwards
        """

        # calculate distance
        d = self.x.distance(self.lon_center, self.lat_center, earth_radius=6371.) / 1000.
        msk = d < self.maxdist
        x = self.x.copy()
        del self.x

        x._apply_mask(msk)
        self.x = x.cut_bounding_box(return_object=True)

        lon, lat, data = self.x.get_valid_data()
        del x

        self._distfiltered = True

    def _check(self):
        if np.any(np.diff(self.lags) < 0.):
            raise ValueError('Bins are not in ascending order!')

        # ensure qual binning  of lags
        di = np.diff(self.lags)
        if np.any(np.abs(di - di[0]) > 1.E-10):
            print di
            print di[0]
            raise ValueError('Only equal bins currently supported"')
        if np.any(np.diff(self.lags) < 0.):
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
        d = np.abs((self.x.lon - self.lon_center) ** 2. + (self.x.lat - self.lat_center) ** 2.)
        dmin = d.min()
        m = d == dmin

        idx = np.indices(d.shape)
        i = idx[0][m][0]
        j = idx[1][m][0]

        if (np.abs(1. - self.x.lon[i, j] / self.lon_center) > 0.05) or (np.abs(1. - self.x.lat[i, j] / self.lat_center) > 0.05):  # at least 5% acc.
            print 'lon: ', self.x.lon[i, j], self.lon_center
            print 'lat: ', self.x.lat[i, j], self.lat_center
            i = None
            j = None
        return i, j

    def plot_semivariogram(self, ax=None, color='red', logy=False, ref_lags=None, fit_variogram=False):
        """
        plot semivariogram

        Parameters
        ----------
        ref_lags : list
            list of reference lags. If given, then these lags are plotted in the figure as well
        """
        #~ if not hasattr(self, '_distance'):
            #~ self._calculate_distance()

        if self.maxdist is not None:
            if not self._distfiltered:
                self._filter_data_maxdist()

        f, ax = self._get_figure_ax(ax)
        V = self.calc_semivariance()
        r = self.statistic['semivariogram']['r']
        sigma = self.statistic['semivariogram']['sigma']

        # plot fitted semivariogram if desired
        if fit_variogram:
            param = V.fit(r, sigma)
            if False:
                print r
                print sigma
                print 'Parameters: ', param
                print param['fit_success']
            V.plot(V._h, V._gamma, ax=ax)  # plots experimental variogram and fitted model
        else:
            if logy:
                ax.semilogy(r, sigma, 'x', color=color)
            else:
                ax.plot(r, sigma, 'x', color=color)

        # set additional labels
        ax.set_ylabel('$\sigma^2$ / 2 (isotropic)')
        ax.set_xlabel('distance from center [km]')
        ax.grid()

        if ref_lags is not None:
            for d in ref_lags:
                ax.plot([d, d], ax.get_ylim(), color='grey')

        return ax

    def plot_percentiles(self, p, ax=None, logy=False, ref_lags=None):
        """
        plot percentiles

        p : list
            list of percentiles [0 ... 1]
        """

        if self.maxdist is not None:
            if not self._distfiltered:
                self._filter_data_maxdist()

        f, ax = self._get_figure_ax(ax)
        for e in p:
            self.calc_percentile(e)
            r = self.statistic['percentiles'][e]['r']
            v = self.statistic['percentiles'][e]['value']
            if logy:
                ax.semilogy(r, v, 'x-', label='p=' + str(e))
            else:
                ax.plot(r, v, 'x-', label='p=' + str(e))
        ax.set_ylabel(self.x._get_unit())
        ax.set_xlabel('distance [km]')
        ax.grid()
        ax.legend(loc='upper left', prop={'size': 10}, ncol=3)
        if ref_lags is not None:
            for d in ref_lags:
                ax.plot([d, d], ax.get_ylim(), color='grey')

        return ax

    def calc_percentile(self, p):
        """
        p [0 ... 1]
        """
        bounds = self.lags
        r = []
        v = []
        for b in bounds:
            d = self._get_data_distance(0., b)
            if d is None:
                continue
            if len(d) < 1:
                continue
            r.append(b)
            v.append(np.percentile(d, p * 100.))  # percentile value

        r = np.asarray(r)
        np.asarray(v)

        o = {'r': np.asarray(r), 'value': np.asarray(v)}
        if 'percentiles' not in self.statistic.keys():
            self.statistic.update({'percentiles': {}})

        self.statistic['percentiles'].update({p: o})

    def _get_data_distance(self, lb, ub):
        """
        get data for distance where

        lb<= r < ub
        """
        if not hasattr(self, '_distance'):
            self._calculate_distance()
        m = (self._distance >= lb) & (self._distance < ub)
        if m.sum() == 0:
            return None
        o = self.x.data[m].flatten()
        if isinstance(o, np.ma.core.MaskedArray):
            o = o.data[~o.mask]  # ensure that nparray is returned
        return o

    def calc_semivariance(self, model='spherical'):
        """
        calculate semivariance for selected range bins

        Parameters
        ----------
        maxdist : float
            maximum distance [km]

        """
        assert self.x.ndim == 2

        # determine first all data within a certain distance only

        # get flattened data
        lon, lat, data = self.x.get_valid_data()


        if type(data) is np.ma.core.MaskedArray:
            # convert to ndarray
            msk = ~data.mask
            if type(lon) is not np.ma.core.MaskedArray:
                lon = np.ma.array(lon, mask=lon!=lon)
            if type(lat) is not np.ma.core.MaskedArray:
                lat = np.ma.array(lat, mask=lat!=lat)

            lon = lon.data[msk]
            lat = lat.data[msk]
            data = data.data[msk]

        # estimate experimental variogram
        if model == 'spherical':
            V = SphericalVariogram()
        else:
            raise ValueError('Invalid variogram type')
        dlag = self.lags[1] - self.lags[0]  # assume equal lag binning

        r, v = V.semivariogram(data.astype('float'), lon.astype('float'), lat.astype('float'), self.lags.astype('float'), dlag)

        # store results
        o = {'r': np.asarray(r), 'sigma': np.asarray(v)}
        self.statistic.update({'semivariogram': o})

        return V

    def _get_figure_ax(self, ax):
        if ax is None:
            f = plt.figure()
            ax = f.add_subplot(111)
        else:
            f = ax.figure
        return f, ax

    def _calculate_distance(self, data=None):
        """
        calculate distance [km]

        Parameters
        ----------
        data : Data
            if given, then this object is taken for distance calculations instead of self.x
        """
        if not hasattr(self, 'lon_center'):
            raise ValueError('ERROR: You need to specify first the center position!')
        if data is None:
            ref = self.x
        else:
            ref = data

        self._distance = ref.distance(self.lon_center, self.lat_center) / 1000.

    def get_coordinates_at_distance(self, d, N=360, dist_threshold=None, oversampling_factor=None):
        """
        return an ordered list of coordinates that are closest to the
        specified radius.

        This can be used to generate e.g. a Polygon for plotting.
        The approach is that a circle is drawn around the center location
        and the all points with cloes radius in each direction is stored

        Parameters
        ----------
        d : float
            distance [km]
        N : int
            number of samples
        dist_threshold : float
            threshold [km] to identify closest points
        oversampling_factor : float
            if None, then nothing happens. If > 1. then the distance
            calculations will be done on an oversampled grid. This allows
            to plot e.g. nicer overlays (e.g. circles) with the data
        """

        if not hasattr(self, 'lon_center'):
            raise ValueError('Missing coordinate!')
        if not hasattr(self, 'lat_center'):
            raise ValueError('Missing coordinate!')
        if dist_threshold is None:
            raise ValueError('You need to provide a distance threshold!')

        if oversampling_factor is None:
            refobj = self.x
            if not hasattr(self, '_distance'):
                self._calculate_distance()
        else:
            if oversampling_factor < 1.:
                raise ValueError('Oversampling factor needs to be > 1!')
            refobj = self.x.copy()

            if not isinstance(refobj.lon, np.ma.core.MaskedArray):
                refobj.lon = np.ma.array(refobj.lon, mask = refobj.lon != refobj.lon)
            if not isinstance(refobj.lat, np.ma.core.MaskedArray):
                refobj.lat = np.ma.array(refobj.lat, mask = refobj.lat != refobj.lat)

            if not isinstance(refobj.lon, np.ma.core.MaskedArray):
                print 'lon', refobj.lon
                print type(refobj.lon)
                assert False, 'ERROR: longitudes are expected to be masked arrays!  get_coordinates'
            if not isinstance(refobj.lat, np.ma.core.MaskedArray):
                print 'lat', refobj.lat
                print type(refobj.lat)
                assert False, 'ERROR: latitudes are expected to be masked arrays! get_coordinates'

            ny = int(refobj.ny * oversampling_factor)
            nx = int(refobj.nx * oversampling_factor)

            # in case that no valid coordinates available, no analysis will be made
            if np.sum(~refobj.lon.mask) < 1:
                print 'No proper coordinates found! (get_coordinates_at_distance) A'
                return None, None

            lonn = np.linspace(refobj.lon.min(), refobj.lon.max(), nx)
            latn = np.linspace(refobj.lat.min(), refobj.lat.max(), ny)
            refobj.lon, refobj.lat = np.meshgrid(lonn, latn)

            self._calculate_distance(data=refobj)

        # get closest points as preselection
        di = np.abs(self._distance - d)
        msk = di < dist_threshold
        lons = refobj.lon[msk].flatten()
        lats = refobj.lat[msk].flatten()
        dist = self._distance[msk].flatten()

        if len(lons) == 0:
            print 'No proper coordinates found! (get_coordinates_at_distance) B'
            return None, None

        theta = np.linspace(0., 2. * np.pi, N)  # angle
        LON = []
        LAT = []
        for t in theta:
            x = self.lon_center + d * np.cos(t)
            y = self.lat_center + d * np.sin(t)

            # search for closest point
            dd = np.sqrt((lons - x) ** 2. + (lats - y) ** 2.)

            LON.append(lons[dd.argmin()])
            LAT.append(lats[dd.argmin()])

        return LON, LAT
