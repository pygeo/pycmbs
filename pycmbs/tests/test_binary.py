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
import numpy as np
import tempfile
import struct

from nose.tools import assert_raises

class TestData(unittest.TestCase):

    def setUp(self):
        self.x = np.zeros((7,9))
        self.x[2:4,5:7] = 1.

    def tearDown(self):
        pass

    def test_read_full_binary_file_double(self):
        # write binary test data
        fname = tempfile.mktemp()
        f = open(fname, 'w')
        f.write(self.x)
        f.close()

        D = Data(None, None)
        D.filename = fname
        ny, nx = self.x.shape
        D._read_binary_file(ny=ny, nx=nx, nt=1, dtype='double')


        # read data again
        #~ f = open(fname, 'r')
        #~ x = f.read()
        #~ f.close()

        # decode
        #~ a = np.asarray(struct.unpack('d'*np.prod(self.x.shape), x))
        #~ print a
        #~ a = a.reshape(self.x.shape)

        self.assertTrue(np.all(D.data-self.x == 0.))


















