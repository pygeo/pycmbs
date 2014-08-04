# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from unittest import TestCase
import unittest

from pycmbs.data import Data
import os
import scipy as sc
import matplotlib.pylab as pl
import numpy as np
from scipy import stats
from dateutil.rrule import rrule
from dateutil.rrule import MONTHLY
import datetime

import tempfile

from nose.tools import assert_raises



class TestData(unittest.TestCase):

    def setUp(self):
        pass

    def test_rasterize_init(self):
        x = Data(None, None)
        x._init_sample_object(ny=1, nx=272)
        x.lon = np.random.random(272)*10. + 5.  # 5 ... 15
        x.lat = np.random.random(272)*20. + 0.  # 0 ... 20

        lon = np.random.random((10,20))
        lat = np.random.random((30,20))

        with self.assertRaises(ValueError):
            x._rasterize(lon, lat, radius=0.1)
        lon = np.random.random((10,20))
        lat = np.random.random((10,20))

        with self.assertRaises(ValueError):
            x._rasterize(lon, lat, radius=None)

    def test_rasterize_data(self):
        """
        testdataset

        +---+---+---+
        |1.2|2.3|   |
        +---+---+---+
        |   |   |0.7|
        +---+---+---+
        |   |5.2|   |
        +---+---+---+
        """
        x = Data(None, None)
        x._init_sample_object(ny=1, nx=272)

        x.lon = np.asarray([2.25, 2.45, 1.8, 3.6])
        x.lat = np.asarray([11.9, 10.1, 10.2, 11.3])
        x.data = np.asarray([5.2, 2.3, 1.2, 0.7])

        # target grid
        lon = np.asarray([1.5, 2.5, 3.5])
        lat = np.asarray([10., 11., 12.])
        LON, LAT = np.meshgrid(lon, lat)

        # rasterize data

        # no valid data
        res = x._rasterize(LON, LAT, radius=0.000001, return_object=True)
        self.assertEqual(res.data.mask.sum(), np.prod(LON.shape))

        with self.assertRaises(ValueError):
            res = x._rasterize(LON, LAT, radius=0.000001, return_object=False)

        # check valid results
        res = x._rasterize(LON, LAT, radius=0.5, return_object=True)
        self.assertEqual(res.data[0,0], 1.2)
        self.assertEqual(res.data[0,1], 2.3)
        self.assertEqual(res.data[1,2], 0.7)
        self.assertEqual(res.ny*res.nx - res.data.mask.sum(), 4)


if __name__ == '__main__':
    unittest.main()
