# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import unittest
import os
from pycmbs.benchmarking import preprocessor

class TestPreprocessor(unittest.TestCase):

    def setUp(self):
        self.preprocessor = preprocessor.CMIP5Preprocessor('.' + os.sep, 'test_file', 'nc', 'model', 'experiment', institute='MPI', mip='Amon', realm='atmos')

    def test_StubTest(self):
        self.assertEqual(1, 1)

if __name__ == "__main__":
    unittest.main()


