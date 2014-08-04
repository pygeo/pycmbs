# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import unittest
import numpy
from pycmbs import mapping
import scipy as sc
from pycmbs.data import Data
import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from nose.tools import assert_raises
from pycmbs.mapping import map_plot
import os

import tempfile

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

        self._tmpdir = tempfile.mkdtemp()

    def test_SingleMap_Init(self):
        try:
            import cartopy.crs as ccrs
        except:
            return True  #no testing if cartopy not installed

        # just test if things pass
        SM1 = mapping.SingleMap(self.D)
        SM2 = mapping.SingleMap(self.D, stat_type='sum')
        proj_prop = {'projection': 'robin'}
        SM3 = mapping.SingleMap(self.D, backend='basemap', stat_type='median')
        SM4 = mapping.SingleMap(self.D, backend='cartopy')
        SM1.plot(show_zonal=True)
        SM2.plot(show_zonal=True, colorbar_orientation='horizontal')
        SM3.plot(show_zonal=True, colorbar_orientation='horizontal', proj_prop=proj_prop)
        SM4.plot(show_zonal=True, colorbar_orientation='horizontal', proj_prop=proj_prop)

    def test_SingleMap_WithoutColorbar(self):
        SM = mapping.SingleMap(self.D)
        SM.plot(show_colorbar=False)

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
        SM = mapping.SingleMap(self.D, savefile=self._tmpdir + os.sep + 'my_test_save_file.nc')
        SM.save(save_mean=True, save_all=True)
        self.assertTrue(os.path.exists(self._tmpdir + os.sep + 'my_test_save_file.nc_timmean.nc'))
        self.assertTrue(os.path.exists(self._tmpdir + os.sep + 'my_test_save_file.nc_timmean.nc_all.nc'))
        if os.path.exists(self._tmpdir + os.sep + 'my_test_save_file.nc_timmean.nc'):
            os.remove(self._tmpdir + os.sep + 'my_test_save_file.nc_timmean.nc')
        if os.path.exists(self._tmpdir + os.sep + 'my_test_save_file.nc_timmean.nc_all.nc'):
            os.remove(self._tmpdir + os.sep + 'my_test_save_file.nc_timmean.nc_all.nc')

    @unittest.skip('skip as only for local testing')
    def test_SingleMap_add_cyclic(self):
        file='/home/m300028/shared/data/SEP/variables/land/Ta_2m/cru_ts_3_00.1901.2006.tmp_miss_t63.nc'
        ofile = 'world.png'
        if os.path.exists(ofile):
            os.remove(ofile)
        d=Data(file,'tmp',read=True)
        map_plot(d, use_basemap=True, savegraphicfile=ofile)
        if os.path.exists(ofile):
            os.remove(ofile)


if __name__ == "__main__":
    unittest.main()

# vim: expandtab shiftwidth=4 softtabstop=4
