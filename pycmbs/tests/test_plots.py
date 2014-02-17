# -*- coding: utf-8 -*-

import unittest
from pycmbs import plots
from pycmbs.data import Data
from pycmbs.plots import ReichlerPlot, ScatterPlot, LinePlot, HistogrammPlot, ZonalPlot
from pycmbs.plots import map_difference, map_season, GlecklerPlot
from pycmbs.plots import xx_map_plot

import scipy
import os
import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as plt

class TestPycmbsPlots(unittest.TestCase):

    def setUp(self):
        n=1000  # slows down significantly! constraint is percentile  test
        x = scipy.randn(n)*100.  # generate dummy data
        self.D = Data(None, None)
        d=np.ones((n, 1, 1))
        self.D.data = d
        self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data, mask=self.D.data != self.D.data)
        self.D.verbose = True
        self.D.unit = 'myunit'
        self.D.label = 'testlabel'
        self.D.filename = 'testinputfilename.nc'
        self.D.varname = 'testvarname'
        self.D.long_name = 'This is the longname'
        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01')
        self.D.time_str = "days since 0001-01-01 00:00:00"
        self.D.calendar = 'gregorian'
        self.D.cell_area = np.ones((1,1))
        self.D.lon=np.ones((1,1)) * 20.
        self.D.lat=np.ones((1,1)) * 10.

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
        G.write_ranking_table('ta', 'nix.tex', fmt='latex')
        self.assertTrue(os.path.exists('nix.tex'))
        if os.path.exists('nix.tex'):
            os.remove('nix.tex')
        G.write_ranking_table('ta', 'nix1', fmt='latex')
        self.assertTrue(os.path.exists('nix1.tex'))
        if os.path.exists('nix1.tex'):
            os.remove('nix1.tex')

    def test_old_map_plot(self):
        xx_map_plot(self.D)
        #~ xx_map_plot(self.D, use_basemap=True)

# map_season






if __name__ == "__main__":
    unittest.main()
# vim: expandtab shiftwidth=4 softtabstop=4
