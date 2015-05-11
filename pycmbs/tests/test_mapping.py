# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
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
        self.D = Data(None, None)
        self.D._init_sample_object(nt=1000, ny=10, nx=20)
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

    def test_invalid_colorbar_orientation(self):
        SM = mapping.SingleMap(self.D)
        with self.assertRaises(ValueError):
            SM.plot(colorbar_orientation='something')

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
