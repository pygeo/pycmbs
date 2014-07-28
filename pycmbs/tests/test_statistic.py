# -*- coding: UTF-8 -*-
from unittest import TestCase

__author__ = 'm300028'

from pycmbs.statistic import *
import numpy as np



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
            e = np.random.random(len(t))
            y = A * np.cos(w*self.t + phase)
            return y, e

        # test with single frequency
        p_ref = 10.
        w = 2.*np.pi / p_ref
        phase = 0.
        amplitude = 5.
        y, e = _sample_data(self.t, w, amplitude, phase)

        P = np.arange(0., 20., 2.)  # target period [days]
        Ar, bR = lomb_scargle_periodogram(self.t, P, y+e)
        self.assertEqual(Ar[5], amplitude)
        self.assertEqual(bR[5], phase)

        # test for data with phase offset
        phase = 1.5
        amplitude = 3.
        y, e = _sample_data(self.t, w, amplitude, phase)

        P = np.arange(0., 20., 2.)  # target period [days]
        Ar, bR = lomb_scargle_periodogram(self.t, P, y+e)
        self.assertEqual(Ar[5], amplitude)
        self.assertEqual(bR[5], phase)


        # test for functions with overlapping frequencies
        p_ref1 = 10.
        p_ref2 = 14.
        w1 = 2.*np.pi / p_ref1
        w2 = 2.*np.pi / p_ref2
        y1, e1 = _sample_data(self.t, w1, 4., 0.)
        y2, e2 = _sample_data(self.t, w2, 2., np.pi*0.5)
        P = np.arange(0., 20., 2.)  # target period [days]
        Ar, bR = lomb_scargle_periodogram(self.t, P, y1+e1+y2+e2)

        self.assertTrue(False)

        # test for data with gaps


