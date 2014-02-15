# -*- coding: utf-8 -*-

import unittest
from pycmbs import plots
from pycmbs.data import Data
from pycmbs.plots import ReichlerPlot, ScatterPlot
import scipy
import numpy as np
import matplotlib.pylab as pl

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

    def test_ReichlerPlotGeneral(self):
        RP = ReichlerPlot()
        for i in xrange(10):
            RP.add([i*12.], 'test'+str(i))
        RP.simple_plot()
        RP.bar(title='some title', vmin=-10., vmax=10.)

    def test_ScatterPlotGeneral(self):
        x = self.D
        S = ScatterPlot(x)
        S.plot(x)



    def test_DummyTest(self):
        pass

if __name__ == "__main__":
    unittest.main()
# vim: expandtab shiftwidth=4 softtabstop=4
