# -*- coding: utf-8 -*-

#
# TESTS with some test data
#
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""


from unittest import TestCase
import unittest

from pycmbs.data import Data
from pycmbs.examples import download

import os
import matplotlib.pylab as pl
import numpy as np

from nose.tools import assert_raises


class TestData(unittest.TestCase):
    def setUp(self):
        self.D = download.get_sample_file(name='air')
        self.file = download.get_sample_file(name='air', return_object=False)  # filename only
        self.areafile = self.file[:-3] + '_cell_area.nc'

    @unittest.skip('wait fixing of CDO problem')
    def test_cell_area(self):
        del self.D.cell_area
        self.D._set_cell_area()

        # check if cell area file is existing
        self.assertTrue(os.path.exists(self.areafile))
        self.assertTrue(hasattr(self.D, 'cell_area'))

        # now remove the cell area file to force calculation
        os.remove(self.areafile)
        del self.D.cell_area
        self.D._set_cell_area()
        self.assertTrue(os.path.exists(self.areafile))

    @unittest.skip('wait fixing of CDO problem')
    def test_cell_area_InvalidLatLon(self):
        self.D.lon = None
        self.D.lat = None
        self.D.cell_area = None
        os.remove(self.areafile)
        self.D._set_cell_area()
        self.assertTrue(np.all(self.D.cell_area == 1.))

    @unittest.skip('wait fixing of CDO problem')
    def test_zonal_mean(self):
        zm = self.D.get_zonal_mean()
        ZM = self.D.get_zonal_mean(return_object=True)
        self.assertTrue(np.all(np.abs(1.-zm / ZM.data.T) < 1.E-6 ))

    @unittest.skip('wait fixing of CDO problem')
    def test_get_deseasonalized_anomaly(self):
        ### TODO implement some reference solution
        c = self.D.get_deseasonalized_anomaly(base='current')
        c = self.D.get_deseasonalized_anomaly(base='all')

if __name__ == '__main__':
    unittest.main()



