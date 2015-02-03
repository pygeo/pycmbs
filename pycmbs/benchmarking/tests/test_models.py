# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import unittest
from pycmbs.benchmarking import models

from pycmbs.data import Data
import os
import scipy as sc
import matplotlib.pylab as pl
import numpy as np
from scipy import stats
from dateutil.rrule import rrule
from dateutil.rrule import MONTHLY
import datetime

from nose.tools import assert_raises

import tempfile

class TestPycmbsBenchmarkingModels(unittest.TestCase):

    def setUp(self):
        self.D = Data(None, None)
        self.D._init_sample_object(nt=1000, ny=1, nx=1)

        # generate dummy Model object
        data_dir = '.' + os.sep + 'test' + os.sep
        varmethods = {'albedo':'get_albedo()', 'sis': 'get_sis()'}
        self.model = models.Model(data_dir, varmethods, name='testmodel', intervals='monthly')

        sis = self.D.copy()
        sis.mulc(5., copy=False)
        sis.label='sisdummy'

        alb = self.D.copy()
        alb.label='albedodummy'

        # add some dummy data variable
        self.model.variables = {'albedo':alb, 'sis':sis}

    def test_save_prefix_missing(self):
        m = self.model
        odir = tempfile.mkdtemp() + os.sep
        with self.assertRaises(ValueError):
            m.save(odir)

    def test_save_create_odir(self):
        m = self.model
        odir = tempfile.mkdtemp() + os.sep
        if os.path.exists(odir):
            os.system('rm -rf ' + odir)
        m.save(odir, prefix='test')
        self.assertTrue(os.path.exists(odir))
        os.system('rm -rf ' + odir)

    def test_save(self):
        m = self.model
        odir = tempfile.mkdtemp() + os.sep

        sisfile = odir + 'testoutput_SIS.nc'
        albfile = odir + 'testoutput_ALBEDO.nc'
        if os.path.exists(sisfile):
            os.remove(sisfile)
        if os.path.exists(albfile):
            os.remove(albfile)

        m.save(odir, prefix='testoutput')
        self.assertTrue(os.path.exists(sisfile))
        self.assertTrue(os.path.exists(albfile))

        if os.path.exists(sisfile):
            os.remove(sisfile)
        if os.path.exists(albfile):
            os.remove(albfile)
        os.system('rm -rf ' + odir)

    def test_cmip5_init_singlemember(self):
        data_dir = tempfile.mkdtemp()

        # invalid model identifier
        with self.assertRaises(ValueError):
            M = models.CMIP5RAW_SINGLE(data_dir, 'MPI-M:MPI-ESM-LR1', 'amip', {}, intervals='monthly')
        with self.assertRaises(ValueError):
            M = models.CMIP5RAW_SINGLE(data_dir, 'MPI-M:MPI-ESM-LR#1#2', 'amip', {}, intervals='monthly')
        M1 = models.CMIP5RAW_SINGLE(data_dir, 'MPI-M:MPI-ESM-LR#1', 'amip', {}, intervals='monthly')
        M2 = models.CMIP5RAW_SINGLE(data_dir, 'MPI-M:MPI-ESM-LR#728', 'amip', {}, intervals='monthly')
        self.assertEqual(M1.ens_member, 1)
        self.assertEqual(M2.ens_member, 728)

    def test_cmip5_singlemember_filename(self):
        data_dir = tempfile.mkdtemp()

        # generate testfile
        testfile = data_dir + os.sep + 'MPI-M' + os.sep + 'MPI-ESM-LR' + os.sep + 'amip' + os.sep + 'mon' + os.sep + 'atmos' + os.sep + 'Amon' + os.sep + 'r1i1p1' + os.sep + 'ta' + os.sep + 'ta_Amon_MPI-ESM-LR_amip_r1i1p1_197901-200812.nc'
        os.makedirs(os.path.dirname(testfile))
        os.system('touch ' + testfile)
        self.assertTrue(os.path.exists(testfile))

        M = models.CMIP5RAW_SINGLE(data_dir, 'MPI-M:MPI-ESM-LR#1', 'amip', {}, intervals='monthly')
        kwargs = {'CMIP5RAWSINGLE' : {'mip' : 'Amon', 'realm' : 'atmos', 'temporal_resolution' : 'mon'}}
        f = M.get_raw_filename('ta', **kwargs)
        self.assertTrue(os.path.exists(f))
        self.assertEqual(f, testfile)

if __name__ == "__main__":
    unittest.main()

