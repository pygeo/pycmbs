# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""
from unittest import TestCase

from pycmbs.statistic import *
from pycmbs.data import Data

import numpy as np
import matplotlib.pyplot as plt



class TestStatistic(TestCase):
    def test_get_significance(self):
        # test get_significance routine
        self.assertAlmostEqual(get_significance(0.5, 100.), 1.18049202704e-07, delta=1.e-6)



class TestLomb(TestCase):
    # note that the tests are not 100percent stable!

    def setUp(self):
        self.t = np.arange(0, 365*3)   # time in days


    def tearDown(self):
        pass


    def test_lomb_basic(self):

        def _sample_data(t, w, A, B):
            e = np.random.random(len(t))*0.
            y = A * np.cos(w*self.t + B)
            return y, e

        def _test_ratio(x,y, thres=0.05):
            r = np.abs(1. - x / y)
            print r, x/y
            self.assertTrue(r <= thres) # accuracy of ration by 5%

        # test with single frequency
        p_ref = 10.
        w = 2.*np.pi / p_ref
        y, e = _sample_data(self.t, w, 5., 0.1)

        P = np.arange(2., 20., 2.)  # target period [days]
        Ar, Br = lomb_scargle_periodogram(self.t, P, y+e, corr=False)

        _test_ratio(Ar[4], 5.)
        _test_ratio(Br[4], 0.1)

        Ar, Br, Rr, Pr = lomb_scargle_periodogram(self.t, P, y)
        _test_ratio(Ar[4], 5.)
        _test_ratio(Br[4], 0.1)
        self.assertEqual(Rr[4], 1.)
        self.assertEqual(Pr[4], 0.)



        # test for functions with overlapping frequencies
        p_ref1 = 365.
        p_ref2 = 365.
        w1 = 2.*np.pi / p_ref1
        w2 = 2.*np.pi / p_ref2
        y1, e1 = _sample_data(self.t, w1, 4., 0.1)
        y2, e2 = _sample_data(self.t, w2, 3.6, 0.1)

        P = np.arange(1., 366., 1.)  # target period [days]
        Ar, Br = lomb_scargle_periodogram(self.t, P, y1+e1+y2+e2, corr=False)

        _test_ratio(Ar[-1], 7.6)
        _test_ratio(Br[-1], 0.1)

        # overlapping frequencies 2
        p_ref1 = 100.
        p_ref2 = 200.
        w1 = 2.*np.pi / p_ref1
        w2 = 2.*np.pi / p_ref2
        y1, e1 = _sample_data(self.t, w1, 2., np.pi*0.3)  # don't choose pi for phase, as this will result in an optimization with negative amplitude and zero phase (= sin)
        y2, e2 = _sample_data(self.t, w2, 3., np.pi*0.5)
        P = np.arange(1., 366., 1.)  # target period [days]
        hlp = y1+e1+y2+e2
        Ar, Br = lomb_scargle_periodogram(self.t, P, hlp, corr=False)

        # sample data object
        D = Data(None, None)
        D._init_sample_object(nt=len(y), ny=1, nx=1)
        D.data[:,0,0] = np.ma.array(hlp, mask=hlp!=hlp)
        D.time = self.t

        D_dummy = Data(None, None)
        D_dummy._init_sample_object(nt=len(y), ny=1, nx=1)
        with self.assertRaises(ValueError):
            D_dummy.time_str = 'hours since 2001-01-01'  # only days currently supported!
            xx, yy = D_dummy.lomb_scargle_periodogram(P, return_object=False)

        AD, BD = D.lomb_scargle_periodogram(P, return_object=False, corr=False)
        AD1, BD1 = D.lomb_scargle_periodogram(P, return_object=True, corr=False)
        self.assertEqual(AD.shape, BD.shape)
        self.assertEqual(D.ny, AD.shape[1])
        self.assertEqual(D.nx, AD.shape[2])

        _test_ratio(Ar[99], 2.)
        _test_ratio(AD[99,0,0], 2.)
        _test_ratio(AD1.data[99, 0,0], 2.)

        _test_ratio(Ar[199], 3.)
        _test_ratio(AD[199,0,0], 3.)
        _test_ratio(AD1.data[199,0,0], 3.)

        _test_ratio(Br[99], np.pi*0.3)
        _test_ratio(BD[99,0,0], np.pi*0.3)
        _test_ratio(BD1.data[99,0,0], np.pi*0.3)

        _test_ratio(Br[199], np.pi*0.5)
        _test_ratio(BD[199,0,0], np.pi*0.5)
        _test_ratio(BD1.data[199,0,0], np.pi*0.5)

        # test for data with gaps
        # tests are not very robust yet as results depend on noise applied!
        p_ref1 = 100.
        p_ref2 = 200.
        w1 = 2.*np.pi / p_ref1
        w2 = 2.*np.pi / p_ref2
        y1, e1 = _sample_data(self.t, w1, 2., np.pi*0.3)  # don't choose pi for phase, as this will result in an optimization with negative amplitude and zero phase (= sin)
        y2, e2 = _sample_data(self.t, w2, 3., np.pi*0.5)
        P = np.arange(1., 366., 1.)  # target period [days]

        ran = np.random.random(len(self.t))
        msk = ran > 0.1
        tmsk = self.t[msk]
        yref = y1+e1+y2+e2
        ymsk = yref[msk]

        Ar, Br = lomb_scargle_periodogram(tmsk, P, ymsk, corr=False)

        _test_ratio(Ar[99], 2., thres=0.1)
        _test_ratio(Ar[199], 3., thres=0.1)
        _test_ratio(Br[99], np.pi*0.3, thres=0.1)
        _test_ratio(Br[199], np.pi*0.5, thres=0.1)


    def test_lomb_normalize(self):
        # LOMB only works with zero mean data !!!!
        # normalization should be therefore implemented, but
        # isn't so far! This was not recognizable with the other
        # tests as these are zero mean anyway!

        self.assertTrue(False)

