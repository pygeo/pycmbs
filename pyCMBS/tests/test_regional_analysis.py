from unittest import TestCase
import unittest

from pyCMBS.data import *
from pyCMBS.diagnostic import RegionalAnalysis
import scipy as sc
#import pylab as pl
import numpy as np
#from scipy import stats
#from pyCMBS.netcdf import *
#from dateutil.rrule import *

class TestData(TestCase):

    def setUp(self):
        # init Data object for testing
        n=1000  # slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100.  # generate dummy data
        self.D = Data(None,None)
        d=np.ones((n,1,1))
        self.D.data = d
        self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data,mask=self.D.data != self.D.data)
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

    def test_regional_analysis(self):

        # generate two datasets
        ny = 2
        nx = 6
        nt = 20
        # regional mask looks like the following
        #
        # | 1 | 2 | 2 | 3 | 4 | 3 |
        # | 1 | 2 | 2 | 4 | 3 | 4 |
        m = np.zeros((2,6))
        m[0, 0] = 1.
        m[0, 1] = 2.
        m[0, 2] = 2.
        m[0, 3] = 3.
        m[0, 4] = 4.
        m[0, 5] = 3.
        m[1, 0] = 1.
        m[1, 1] = 2.
        m[1, 2] = 2.
        m[1, 3] = 4.
        m[1, 4] = 3.
        m[1, 5] = 4.

        cell_area = np.ones_like(m)

        # generate mask
        x = self.D.copy()

        tmp = np.random.random((nt,ny,nx))
        x.data = np.ma.array(tmp, mask = tmp != tmp)
        x.cell_area = cell_area.copy()
        x.time = x.time[0:nt]
        del tmp

        y = self.D.copy()
        tmp = np.random.random((nt,ny,nx))
        y.data = np.ma.array(tmp, mask = tmp != tmp)
        y.cell_area = cell_area.copy()
        y.time = y.time[0:nt]
        del tmp

        # todo unittest extension for different area weighting !!!!

        reg = Data(None,None)
        reg.data = m
        REGSTAT = RegionalAnalysis(x, y, reg)
        REGSTAT.calculate()
        print REGSTAT.statistics.keys()




        #self.assertEqual(s1,'2001-01-05 00:00:00+00:00')
        #self.assertEqual(s2,'2001-05-05 00:00:00+00:00')
