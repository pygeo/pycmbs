# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from unittest import TestCase
import unittest

from nose.tools import assert_raises


import numpy as np

from pycmbs.data import Data
from pycmbs.geostatistic import Geostatistic, Variogram


class TestData(unittest.TestCase):

    def setUp(self):
        D = Data(None, None)
        D.data = np.random.random((55, 20))
        lon = np.arange(-10.,10.)  # -10 ... 9
        lat = np.arange(-60., 50., 2.)  # -60 ... 48
        D.lon, D.lat = np.meshgrid(lon, lat)
        self.x = D

    def tearDown(self):
        pass

    def test_invalid_geometry(self):
        x = self.x
        x.data = np.random.random((5,4,3))
        with self.assertRaises(ValueError):
            G = Geostatistic(x)

    def test_init(self):
        with self.assertRaises(ValueError):  # missing range bins
            G = Geostatistic(self.x)

        # 3D is invalid
        y = self.x.copy()
        y.data = np.random.random((10,20,30))
        with self.assertRaises(ValueError):  # missing range bins
            G = Geostatistic(self.x, range_bins = np.random.random(10))

        y = self.x.copy()  # missing 3D
        y.data = np.random.random((10,20,30))
        with self.assertRaises(ValueError):
            G = Geostatistic(self.x)

        bins = np.random.random(10)
        with self.assertRaises(ValueError):
            G = Geostatistic(self.x, range_bins=bins)
        bins = np.arange(10)
        G = Geostatistic(self.x, range_bins=bins)

    def test_plot(self):
        bins = np.arange(3) / 6.
        G = Geostatistic(self.x, range_bins=bins)
        G.set_center_position(0., 0.)
        G.plot_semivariogram()

    def test_percentiles(self):
        bins = np.arange(10)
        G = Geostatistic(self.x, range_bins=bins)
        G.set_center_position(5., -20.)
        p = [0.05, 0.1, 0.5]
        G.plot_percentiles(p, ax=None, logy=False, ref_lags=None)

    def test_set_center(self):
        bins = np.arange(10)
        G = Geostatistic(self.x, range_bins=bins)
        G.set_center_position(5., -20.)
        self.assertEqual(G.lon_center, 5.)
        self.assertEqual(G.lat_center, -20.)

    def test_center_pos(self):
        bins = np.arange(10)
        G = Geostatistic(self.x, range_bins=bins)
        # forgot to set center
        with self.assertRaises(ValueError):
            G._get_center_pos()

        G.set_center_position(-10., -60.)
        i_lat, i_lon = G._get_center_pos()
        self.assertEqual(i_lat, 0)
        self.assertEqual(i_lon, 0)

        G.set_center_position(9., 48.)
        i_lat, i_lon = G._get_center_pos()
        self.assertEqual(i_lat, 54)
        self.assertEqual(i_lon, 19)

        G.set_center_position(-10.1, 48.)
        i_lat, i_lon = G._get_center_pos()
        self.assertEqual(i_lat, 54)
        self.assertEqual(i_lon, 0)

    def test_get_coordinate_at_distance(self):
        bins = np.arange(10)
        G = Geostatistic(self.x, range_bins=bins)

        with self.assertRaises(ValueError):  # missing center
            lon, lat = G.get_coordinates_at_distance(2.5, dist_threshold=1.)
        G.set_center_position(5., -20.)
        lon, lat = G.get_coordinates_at_distance(2.5, dist_threshold=1.)


    def test_variogram_semivariance(self):

        V = Variogram()
        #~ V._semivariance(self, x, lon, lat, h_km, dh_km)
        #~ _semivariance(self, x, lon, lat, h, dh):

        #~ assert False

    def test_variogram_paired_distance(self):
        V = Variogram()
        #~ lon =
        #~ lat =
        # orthodrome calculations
        #~ pd = V._paired_distance()
        #~ self.assertEqual(...)

    def test_semivariogram_calculation(self):
        V = Variogram()
        x = np.random.random(100)
        lon = np.random.random(100)
        lat = np.random.random(100)
        lags = [1,2,3,4]
        V.semivariogram(x, lon, lat, lags)

    def test_variogram_orthodrome(self):
        # example distance Berlin-Tokio
        # http://de.wikipedia.org/wiki/Orthodrome

        lat_berlin = 52.517
        lon_berlin = 13.4
        lat_tokio = 35.70
        lon_tokio = 139.767

        V = Variogram()
        r = V._orthodrome(lon_tokio, lat_tokio, lon_berlin, lat_berlin, radius=6370.*1000.)
        print r, r-8918000.
        self.assertTrue(abs(r-8918000.)<1000.)



