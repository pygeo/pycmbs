# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import unittest
from pycmbs import netcdf
import numpy as np

class TestPycmbsNetcdf(unittest.TestCase):

    def setUp(self):
        pass

    def test_netCDFHandlerinit_WithInvalidValue(self):
        with self.assertRaises(ValueError):
            cdf = netcdf.NetCDFHandler(netcdf_backend='invalid_backend')

    def test_netCDF_openFile_InvalidMode(self):
        with self.assertRaises(ValueError):
            cdf = netcdf.NetCDFHandler()
            cdf.open_file('nix.nc', 'xxx')

    def test_netCDF_openFile_InvalidFile(self):
        with self.assertRaises(ValueError):
            cdf = netcdf.NetCDFHandler()
            cdf.open_file('nopefile.nc', 'r')

    def test_netCDFHandlerinit_Default(self):
        cdf = netcdf.NetCDFHandler()


if __name__ == "__main__":
    unittest.main()

# vim: expandtab shiftwidth=4 softtabstop=4
