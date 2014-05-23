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
        self.ymin = 2
        self.ymax = 6
        self.xmin = 3
        self.xmax = 8
        r = (np.random.random(self.x.shape)*100.).astype('int').astype('float')
        self.x[self.ymin:self.ymax,self.xmin:self.xmax] = r[self.ymin:self.ymax,self.xmin:self.xmax]
        ny, nx = self.x.shape
        self.lat = np.arange(ny)*10 + 10.
        self.lon = np.arange(nx)*20. - 150.


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
        self.assertTrue(np.all(D.data-self.x == 0.))

    def test_read_full_binary_file_int16(self):
        # write binary test data
        fname = tempfile.mktemp()
        self.x = np.asarray(self.x).astype('int16')
        f = open(fname, 'w')
        f.write(self.x)
        f.close()

        D = Data(None, None)
        D.filename = fname
        ny, nx = self.x.shape
        D._read_binary_file(ny=ny, nx=nx, nt=1, dtype='int16')
        self.assertTrue(np.all(D.data-self.x == 0.))

    def test_read_full_binary_file_int32(self):
        # write binary test data
        fname = tempfile.mktemp()
        self.x = np.asarray(self.x).astype('int32')
        f = open(fname, 'w')
        f.write(self.x)
        f.close()

        D = Data(None, None)
        D.filename = fname
        ny, nx = self.x.shape
        D._read_binary_file(ny=ny, nx=nx, nt=1, dtype='int32')
        self.assertTrue(np.all(D.data-self.x == 0.))



    def test_read_binary_subset_double(self):
        fname = tempfile.mktemp()
        f = open(fname, 'w')
        f.write(self.x)
        f.close()

        D = Data(None, None)
        f = open(fname, 'r')
        ny, nx = self.x.shape
        nt = 1

        # test 1: read entire file
        file_content = D._read_binary_subset2D(f, 8, ny=ny, nx=nx, xbeg=0, xend=nx, ybeg=0, yend=ny)
        d = np.reshape(np.asarray(struct.unpack('d'*ny*nx*nt, file_content)), (ny, nx))
        self.assertTrue(np.all(d-self.x == 0.))

        # test 2: read subset with 1-values only
        ny1 = self.ymax - self.ymin
        nx1 = self.xmax - self.xmin
        nt1 = 1

        file_content = D._read_binary_subset2D(f, 8, ny=ny, nx=nx, xbeg=self.xmin, xend=self.xmax, ybeg=self.ymin, yend=self.ymax)
        d1 = np.reshape(np.asarray(struct.unpack('d'*ny1*nx1*nt1, file_content)), (ny1, nx1))
        self.assertTrue(np.all(d1 - self.x[self.ymin:self.ymax, self.xmin:self.xmax] == 0.))

    def test_read_binary_subset_int(self):
        # INT16 = H
        fname = tempfile.mktemp()
        f = open(fname, 'w')
        ref = (self.x*10).astype('int16')
        f.write(ref)
        f.close()

        D = Data(None, None)
        f = open(fname, 'r')
        ny, nx = self.x.shape
        nt = 1

        # test 1: read entire file
        file_content = D._read_binary_subset2D(f, 2, ny=ny, nx=nx, xbeg=0, xend=nx, ybeg=0, yend=ny)
        d = np.reshape(np.asarray(struct.unpack('H'*ny*nx*nt, file_content)), (ny, nx))
        self.assertTrue(np.all(d-ref == 0.))

        # test 2: read subset with 1-values only
        ny1 = self.ymax - self.ymin
        nx1 = self.xmax - self.xmin
        nt1 = 1

        file_content = D._read_binary_subset2D(f, 2, ny=ny, nx=nx, xbeg=self.xmin, xend=self.xmax, ybeg=self.ymin, yend=self.ymax)
        d1 = np.reshape(np.asarray(struct.unpack('H'*ny1*nx1*nt1, file_content)), (ny1, nx1))
        self.assertTrue(np.all(d1 - ref[self.ymin:self.ymax, self.xmin:self.xmax] == 0.))

    def test_read_binary_subset_Data_double(self):
        # binary data from subset in Data object

        # write binary test data
        fname = tempfile.mktemp()
        f = open(fname, 'w')

        tmp = self.x
        f.write(tmp)
        f.close()

        D = Data(None, None)
        D.filename = fname
        ny, nx = self.x.shape

        latmin = self.lat[self.ymin]
        latmax = self.lat[self.ymax]
        lonmin = self.lon[self.xmin]
        lonmax = self.lon[self.xmax]

        D._read_binary_file(nt=1, dtype='double', latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, lat=self.lat, lon=self.lon)
        self.assertTrue(np.all(D.data-tmp[self.ymin:self.ymax+1,self.xmin:self.xmax+1] == 0.))

    def test_read_binary_subset_Data_int(self):
        # binary data from subset in Data object

        # write binary test data
        fname = tempfile.mktemp()
        f = open(fname, 'w')

        tmp = (np.random.random(self.x.shape)*100.).astype('int16')
        f.write(tmp)
        f.close()

        D = Data(None, None)
        D.filename = fname
        ny, nx = self.x.shape

        latmin = self.lat[self.ymin]
        latmax = self.lat[self.ymax]
        lonmin = self.lon[self.xmin]
        lonmax = self.lon[self.xmax]

        D._read_binary_file(nt=1, dtype='int16', latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, lat=self.lat, lon=self.lon)
        self.assertTrue(np.all(D.data-tmp[self.ymin:self.ymax+1,self.xmin:self.xmax+1] == 0.))














