# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from unittest import TestCase
import unittest

from nose.tools import assert_raises
from pycmbs.geostatistic import Geostatistic

import numpy as np
from pycmbs.data import Data


class TestData(unittest.TestCase):

    def setUp(self):
        D = Data(None, None)
        D.data = np.random.random((10, 20))
        lon = np.arange(-10.,10.)  # -10 ... 9
        lat = np.arange(-60., 50., 2.)  # -60 ... 48
        D.lon, D.lat = np.meshgrid(lon, lat)
        self.x = D

    def tearDown(self):
        pass

    def test_init(self):
        with self.assertRaises(ValueError):
            G = Geostatistic(self.x)

        bins = np.random.random(10)
        with self.assertRaises(ValueError):
            G = Geostatistic(self.x, range_bins=bins)
        bins = np.arange(10)
        G = Geostatistic(self.x, range_bins=bins)

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


