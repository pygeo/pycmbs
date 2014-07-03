# -*- coding: utf-8 -*-

import unittest
from pycmbs.plots.violin import _classic_example as Violin_example
import tempfile

class TestPycmbsPlots(unittest.TestCase):

    def setUp(self):
        pass

    def test_violin_plot(self):
        Violin_example()

if __name__ == "__main__":
    unittest.main()

