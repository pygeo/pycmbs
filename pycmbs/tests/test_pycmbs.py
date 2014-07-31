# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import unittest
import pycmbs

class testPycmbs(unittest.TestCase):
    def setUp(self):
        pass

    def testPycmbsName(self):
        import pycmbs as module
        ref_name = 'pycmbs'
        test_name = module.__name__
        self.assertEqual(ref_name, test_name)


if __name__ == "__main__":
    unittest.main()
