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
        n=1000  # slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100.  # generate dummy data
        self.D = Data(None, None)
        d=np.ones((n, 1, 1))
        self.D.data = d
        self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data, mask=self.D.data != self.D.data)
        self.D.verbose = True
        self.D.unit = 'myunit'
        self.D.label = 'testlabel'
        self.D.filename = 'testinputfilename.nc'
        self.D.varname = 'testvarname'
        self.D.long_name = 'This is the longname'
        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01')
        self.D.time_str = "days since 0001-01-01 00:00:00"
        self.D.calendar = 'gregorian'
        self.D.oldtime=False

        # generate dummy Model object
        data_dir = './test/'
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


if __name__ == "__main__":
    unittest.main()

