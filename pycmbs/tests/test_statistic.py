# -*- coding: UTF-8 -*-
from unittest import TestCase

__author__ = 'm300028'

from pycmbs.statistic import *
import numpy as np
import matplotlib.pyplot as plt



class TestStatistic(TestCase):
    def test_get_significance(self):
        # test get_significance routine
        self.assertAlmostEqual(get_significance(0.5, 100.), 1.18049202704e-07, delta=1.e-6)



class TestLomb(TestCase):

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
        Ar, Br = lomb_scargle_periodogram(self.t, P, y+e)

        _test_ratio(Ar[4], 5.)
        _test_ratio(Br[4], 0.1)

        # test for functions with overlapping frequencies
        p_ref1 = 365.
        p_ref2 = 365.
        w1 = 2.*np.pi / p_ref1
        w2 = 2.*np.pi / p_ref2
        y1, e1 = _sample_data(self.t, w1, 4., 0.1)
        y2, e2 = _sample_data(self.t, w2, 3.6, 0.1)

        P = np.arange(1., 366., 1.)  # target period [days]
        Ar, Br = lomb_scargle_periodogram(self.t, P, y1+e1+y2+e2)

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
        Ar, Br = lomb_scargle_periodogram(self.t, P, y1+e1+y2+e2)

        _test_ratio(Ar[99], 2.)
        _test_ratio(Ar[199], 3.)
        _test_ratio(Br[99], np.pi*0.3)
        _test_ratio(Br[199], np.pi*0.5)

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

        Ar, Br = lomb_scargle_periodogram(tmsk, P, ymsk)

        _test_ratio(Ar[99], 2., thres=0.1)
        _test_ratio(Ar[199], 3., thres=0.1)
        _test_ratio(Br[99], np.pi*0.3, thres=0.1)
        _test_ratio(Br[199], np.pi*0.5, thres=0.1)

