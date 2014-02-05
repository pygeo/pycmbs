#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy
from pycmbs import mapping 

class TestMapPlotGeneric(unittest.TestCase):
    def setUp(self):
        self.map_plot = mapping.MapPlotGeneric() 
    
    def test_StubTest(self):
        self.assertEqual(1, 1)

if __name__ == "__main__":
    unittest.main()

# vim: expandtab shiftwidth=4 softtabstop=4
