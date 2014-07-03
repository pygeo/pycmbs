# -*- coding: utf-8 -*-

import unittest
#~ from pycmbs import plots
#~ from pycmbs.data import Data
#~ from pycmbs.plots import ReichlerPlot, ScatterPlot, LinePlot, HistogrammPlot, ZonalPlot
#~ from pycmbs.plots import map_difference, map_season, GlecklerPlot
#~ from pycmbs.plots import xx_map_plot, HstackTimeseries, HovmoellerPlot
#~ from pycmbs.plots import rotate_ticks, CorrelationAnalysis
from pycmbs.plots.violin import _classic_example as Violin_example

#~ import scipy
#~ import os
#~ import numpy as np
#~ import matplotlib.pylab as pl
#~ import matplotlib.pyplot as plt
import tempfile

class TestPycmbsPlots(unittest.TestCase):

    def setUp(self):
        pass


    def test_violin_plot(self):
        Violin_example()




if __name__ == "__main__":
    unittest.main()

