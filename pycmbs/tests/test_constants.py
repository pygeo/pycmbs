#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
from pycmbs import constants

class TestConstants(unittest.TestCase):
    
    def setUp(self):
        pass

    def testPycmbsConstants(self):
        test_earth_radius = constants.EarthRadius
        ref_earth_radius = 6371000. 
        self.assertEqual(test_earth_radius, ref_earth_radius)

if __name__ == "__main__":
  unittest.main()
# vim: expandtab shiftwidth=4 softtabstop=4
    


