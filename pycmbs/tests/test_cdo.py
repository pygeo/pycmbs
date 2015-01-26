# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from nose.tools import assert_raises
from unittest import TestCase
import unittest
from cdo import Cdo
from pycmbs.examples import download
import os

import tempfile

class TestData(unittest.TestCase):
    def setUp(self):
        self.D = download.get_sample_file(name='air')
        self.file = download.get_sample_file(name='air', return_object=False)  # filename only
        self.areafile = self.file[:-3] + '_cell_area.nc'
        self._tmpdir = tempfile.mkdtemp()

    @unittest.skip('skip test as it is causing some installation related error on Travis')  # TODO
    def test_cdo_general(self):
        # test if cdos work in general
        cdo = Cdo()

