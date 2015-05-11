# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""


from unittest import TestCase
import unittest

from pycmbs.data import Data
from pycmbs.netcdf import NetCDFHandler

import os
import numpy as np
import tempfile

from nose.tools import assert_raises

class TestData(unittest.TestCase):

    def setUp(self):
        self.nx = 20
        self.ny = 10
        self.tempfile = tempfile.mktemp(suffix='.nc')
        self.gfile1 = tempfile.mktemp(suffix='.nc')
        self.gfile2 = tempfile.mktemp(suffix='.nc')
        self.gfile3 = tempfile.mktemp(suffix='.nc')
        self.x = Data(None, None)
        self.x._init_sample_object(nt=10, ny=self.ny, nx=self.nx)
        self.x.save(self.tempfile, varname='myvar')

        # generate some arbitrary geometry file
        F = NetCDFHandler()
        F.open_file(self.gfile1, 'w')
        F.create_dimension('ny', size=self.ny)
        F.create_dimension('nx', size=self.nx)
        F.create_variable('lat', 'd', ('ny', 'nx'))
        F.create_variable('lon', 'd', ('ny', 'nx'))
        F.assign_value('lat', np.ones((self.ny,self.nx)) * 5.)
        F.assign_value('lon', np.ones((self.ny,self.nx)) * 3.)
        F.close()

        F = NetCDFHandler()
        F.open_file(self.gfile2, 'w')
        F.create_dimension('ny', size=self.ny)
        F.create_dimension('nx', size=self.nx)
        F.create_variable('latitude', 'd', ('ny', 'nx'))
        F.create_variable('longitude', 'd', ('ny', 'nx'))
        F.assign_value('latitude', np.ones((self.ny,self.nx)) * 7.)
        F.assign_value('longitude', np.ones((self.ny,self.nx)) * 8.)
        F.close()

        F = NetCDFHandler()
        F.open_file(self.gfile3, 'w')
        F.create_dimension('ny', size=self.ny*2)
        F.create_dimension('nx', size=self.nx*3)
        F.create_variable('latitude', 'd', ('ny', 'nx'))
        F.create_variable('longitude', 'd', ('ny', 'nx'))
        F.assign_value('latitude', np.ones((self.ny*2,self.nx*3)) * 7.)
        F.assign_value('longitude', np.ones((self.ny*2,self.nx*3)) * 8.)
        F.close()

    def test_read_coordinates(self):
        # read data normal
        x1 = Data(self.tempfile, 'myvar', read=True)
        self.assertEqual(x1.nx,self.nx)
        self.assertEqual(x1.ny,self.ny)

        # read data with separate geometry file 'lat', 'lon' names
        x2 = Data(self.tempfile, 'myvar', read=True, geometry_file=self.gfile1)
        self.assertTrue(np.all(x2.lat == 5.))
        self.assertTrue(np.all(x2.lon == 3.))

        # read data with separate geometry file 'latitude', 'longitude' names
        x3 = Data(self.tempfile, 'myvar', read=True, geometry_file=self.gfile2)
        self.assertTrue(np.all(x3.lat == 7.))
        self.assertTrue(np.all(x3.lon == 8.))

        # read data with separate geometry file 'lat', 'lon' names, invalid geometry
        with self.assertRaises(ValueError):
            x4 = Data(self.tempfile, 'myvar', read=True, geometry_file=self.gfile3)

