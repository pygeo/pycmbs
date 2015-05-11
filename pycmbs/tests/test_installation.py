# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

from unittest import TestCase
import unittest

import os
import tempfile

from nose.tools import assert_raises


class TestData(unittest.TestCase):

    def setUp(self):
        self.ogr_test = False  #specify if test shall be performed at all


    def test_gdal(self):
        if self.ogr_test:
            try:
                from osgeo import gdal
                sucess = True
            except:
                sucess = False

            self.assertTrue(sucess)

    def test_ogr(self):
        if self.ogr_test:
            try:
                from osgeo import ogr
                sucess = True
            except:
                sucess = False

            self.assertTrue(sucess)

if __name__ == '__main__':
    unittest.main()
