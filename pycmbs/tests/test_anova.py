# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import numpy
import unittest
from pycmbs import anova

class TestAnova1(unittest.TestCase):
    def setUp(self):
        x_2d_array = numpy.ones((10, 10))
        self.anova = anova.Anova1(x_2d_array)

    def test_AnovaImportWorks(self):
        self.assertEqual(1,1)

if __name__ == "__main__":
    unittest.main()

# vim: expandtab shiftwidth=4 softtabstop=4
