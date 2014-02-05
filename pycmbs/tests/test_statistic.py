from unittest import TestCase

__author__ = 'm300028'

from pyCMBS.statistic import *

class TestStatistic(TestCase):
    def test_get_significance(self):
        #test get_significance routine
        self.assertAlmostEqual(get_significance(0.5,100.),1.18049202704e-07,delta=1.e-6)

