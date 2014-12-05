# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from unittest import TestCase
import unittest

import os
import tempfile

from nose.tools import assert_raises


class TestData(unittest.TestCase):

    def setUp(self):
        pass

    def test_gdal(self):
        try:
            from osgeo import gdal
            sucess = True
        except:
            sucess = False

        self.assertTrue(sucess)

    def test_ogr(self):
        try:
            from osgeo import ogr
            sucess = True
        except:
            sucess = False

        self.assertTrue(sucess)

if __name__ == '__main__':
    unittest.main()
