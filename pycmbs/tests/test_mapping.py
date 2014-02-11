#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy
from pycmbs import mapping
import scipy as sc
from pycmbs.data import Data
import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from nose.tools import assert_raises
import os

class TestMapPlotGeneric(unittest.TestCase):
    def setUp(self):
        self.map_plot = mapping.MapPlotGeneric()

        # data for testing
        n=1000  # slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100.  # generate dummy data
        self.D = Data(None, None)
        ny=10
        nx=20
        d=np.ones((n, ny, nx))
        self.D.data = d
        for i in xrange(ny):
            for j in xrange(nx):
                self.D.data[:, i, j] = x[:]
        self.D.data = np.ma.array(self.D.data, mask=self.D.data != self.D.data)
        self.D.verbose = True
        self.D.unit = 'myunit'
        self.D.label = 'testlabel'
        self.D.filename = 'testinputfilename.nc'
        self.D.varname = 'testvarname'
        self.D.long_name = 'This is the longname'
        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01')
        self.D.time_str = "days since 0001-01-01 00:00:00"
        self.D.calendar = 'gregorian'
        self.D.cell_area = np.ones((ny, nx))
        self.D.lon = np.random.random((ny,nx))*10.
        self.D.lat = np.random.random((ny,nx))*20.

    def test_SingleMap_Init(self):
        # just test if things pass
        SM1 = mapping.SingleMap(self.D)
        SM2 = mapping.SingleMap(self.D)
        proj_prop = {'projection': 'robin'}
        SM3 = mapping.SingleMap(self.D, backend='basemap')
        SM4 = mapping.SingleMap(self.D, backend='cartopy')
        SM1.plot(show_zonal=True)
        SM2.plot(show_zonal=True, colorbar_orientation='horizontal')
        SM3.plot(show_zonal=True, colorbar_orientation='horizontal', proj_prop=proj_prop)
        SM4.plot(show_zonal=True, colorbar_orientation='horizontal', proj_prop=proj_prop)

    def test_SingleMap_WithPredefinedAxis(self):
        f = plt.figure()
        ax = f.add_subplot(2,1,1)
        SM1 = mapping.SingleMap(self.D, ax=ax)

    def test_SingleMap_WithPredefinedAxisButWhichIsNone(self):
        ax = None
        SM1 = mapping.SingleMap(self.D, ax=ax)

    def test_SingleMap_InitWithInvalidBackend(self):
        # just test if things pass
        with self.assertRaises(ValueError):
            SM = mapping.SingleMap(self.D, backend='some_invalid_backend')

    def test_SingleMap_InitWithMissingProjectionProperties(self):
        # just test if things pass
        with self.assertRaises(ValueError):
            SM = mapping.SingleMap(self.D, backend='cartopy')
            SM.plot()

    def test_SingleMap_Save(self):
        SM = mapping.SingleMap(self.D, savefile='my_test_save_file.nc')
        SM.save(save_mean=True, save_all=True)
        self.assertTrue(os.path.exists('my_test_save_file.nc_timmean.nc'))
        self.assertTrue(os.path.exists('my_test_save_file.nc_timmean.nc_all.nc'))

if __name__ == "__main__":
    unittest.main()

# vim: expandtab shiftwidth=4 softtabstop=4
