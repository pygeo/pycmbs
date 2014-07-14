# -*- coding: utf-8 -*-

import unittest
from pycmbs import plots
from pycmbs.data import Data
from pycmbs.plots import ReichlerPlot, ScatterPlot, LinePlot, HistogrammPlot, ZonalPlot
from pycmbs.plots import map_difference, map_season, GlecklerPlot
from pycmbs.plots import xx_map_plot, HstackTimeseries, HovmoellerPlot
from pycmbs.plots import rotate_ticks, CorrelationAnalysis
from pycmbs.plots.violin import _classic_example as Violin_example
from pycmbs.plots import GlobalMeanPlot

import scipy
import os
import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import tempfile

class TestPycmbsPlots(unittest.TestCase):

    def setUp(self):
        self.D = Data(None, None)
        self.D._init_sample_object(nt=1000, ny=1, nx=1)
        self._tmpdir = tempfile.mkdtemp()

    def test_ReichlerPlotGeneral(self):
        RP = ReichlerPlot()
        for i in xrange(10):
            RP.add([i*12.], 'test'+str(i))
        RP.simple_plot()
        RP.bar(title='some title', vmin=-10., vmax=10.)
        #~ RP.circle_plot()

    def test_ScatterPlot_General(self):
        x = self.D
        S = ScatterPlot(x)
        S.plot(x)
        S.legend()

    def test_rotate_ticks(self):
        f = plt.figure()
        ax=f.add_subplot(111)
        ax.plot(np.random.random(1000))
        rotate_ticks(ax, 20.)

    def test_correlation_analysis(self):
        x = self.D
        y = self.D
        C = CorrelationAnalysis(x, y)
        C.do_analysis()

    def test_ScatterPlot_GeneralWithNormalization(self):
        x = self.D
        S = ScatterPlot(x, normalize_data=True)
        S.plot(x)
        S.legend()

    def test_ScatterPlot_FldemeanFalse(self):
        x = self.D
        S = ScatterPlot(x)
        S.plot(x, fldmean=False)
        S.legend()

    def test_ScatterPlot_InvalidShape(self):
        x = self.D
        S = ScatterPlot(x)
        y = self.D.copy()
        y.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            S.plot(y, fldmean=False)

    def test_LinePlot_General(self):
        x = self.D
        L = LinePlot()
        L1 = LinePlot(regress=True)
        L.plot(x)
        L1.plot(x)

    def test_LinePlot_WithAxis(self):
        x = self.D
        f = plt.figure()
        ax = f.add_subplot(111)
        L = LinePlot(ax=ax)
        L.plot(x)

    def test_HistogrammPlot_General(self):
        H = HistogrammPlot(normalize=True)
        H.plot(self.D, bins=10, shown=True)

    def test_ZonalPlot(self):
        Z = ZonalPlot()
        Z.plot(self.D)

    def test_map_difference_General(self):
        map_difference(self.D, self.D)


    def test_GlecklerPlot_InvalidNumberOfObservations(self):
        G = GlecklerPlot()
        G.add_model('echam5')
        G.add_model('mpi-esm')
        G.add_variable('ta')
        G.add_data('ta', 'echam5', 0.5,pos=1)
        G.add_data('ta', 'echam5',0.25,pos=2)
        G.add_data('ta', 'echam5',-0.25,pos=3)
        G.add_data('ta', 'mpi-esm',-0.25,pos=4)
        G.add_data('ta', 'mpi-esm',-0.25,pos=5)
        with self.assertRaises(ValueError):
            G.plot()


    def test_GlecklerPlot_4obs(self):
        G = GlecklerPlot()
        G.add_model('echam5')
        G.add_model('mpi-esm')
        G.add_variable('ta')
        G.add_variable('P')
        G.add_data('P', 'echam5', 0.5,pos=1)
        G.add_data('P', 'echam5',0.25,pos=2)
        G.add_data('P', 'echam5',-0.25,pos=3)
        G.add_data('P', 'mpi-esm',-0.25,pos=4)
        G.plot()

    def test_GlecklerPlot(self):
        G = GlecklerPlot()
        G.add_model('echam5')
        G.add_model('mpi-esm')
        G.add_variable('ta')
        G.add_variable('P')
        G.add_data('ta', 'echam5', 0.5,pos=1)
        G.add_data('P', 'echam5',0.25,pos=1)
        G.add_data('P', 'echam5',-0.25,pos=2)
        G.add_data('P', 'mpi-esm',-0.25,pos=1)
        G.plot()

        G.plot_model_error('ta')
        G.plot_model_ranking('ta')
        G.write_ranking_table('ta', self._tmpdir + os.sep + 'nix.tex', fmt='latex')
        self.assertTrue(os.path.exists(self._tmpdir + os.sep + 'nix.tex'))
        if os.path.exists(self._tmpdir + os.sep + 'nix.tex'):
            os.remove(self._tmpdir + os.sep + 'nix.tex')
        G.write_ranking_table('ta', self._tmpdir + os.sep + 'nix1', fmt='latex')
        self.assertTrue(os.path.exists(self._tmpdir + os.sep + 'nix1.tex'))
        if os.path.exists(self._tmpdir + os.sep + 'nix1.tex'):
            os.remove(self._tmpdir + os.sep + 'nix1.tex')

    def test_old_map_plot(self):
        xx_map_plot(self.D)
        #~ xx_map_plot(self.D, use_basemap=True)


    def test_HstackTimeSeries(self):
        HT = HstackTimeseries()
        for i in xrange(15):
            x = np.random.random(100)*2.-1.
            HT.add_data(x, 'model' + str(i).zfill(3) )
        HT.plot(cmap='RdBu_r', interpolation='nearest', vmin=-1., vmax=1., nclasses=15, title='Testtitle')


    def test_Hovmoeller(self):
        H = HovmoellerPlot(self.D)
        with self.assertRaises(ValueError):
            H.plot()
        H.plot(climits=[0., 1.])

    def test_violin_plot(self):
        Violin_example()


    def test_globalmeanplot(self):
        G = GlobalMeanPlot()
        with self.assertRaises(ValueError):
            G.plot(self.D, stat_type='no_stat_type')
        G.plot(self.D, show_std=True)
        G.plot_mean_result()

# map_season






if __name__ == "__main__":
    unittest.main()
# vim: expandtab shiftwidth=4 softtabstop=4
