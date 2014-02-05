# -*- coding: utf-8 -*-

import unittest
from pycmbs.benchmarking import preprocessor

class TestPreprocessor(unittest.TestCase):

    def setUp(self):
        self.preprocessor = preprocessor.CMIP5Preprocessor('./', 'test_file', 'nc', 'model', 'experiment')

    def test_StubTest(self):
        self.assertEqual(1, 1)

if __name__ == "__main__":
    unittest.main()

# vim: expandtab shiftwidth=4 softtabstop=4
