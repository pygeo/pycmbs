# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the files
LICENSE.md and COPYRIGHT.md
"""

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
    


