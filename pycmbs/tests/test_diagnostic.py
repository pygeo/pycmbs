from unittest import TestCase

__author__ = 'm300028'


from pyCMBS.data import *
from pyCMBS.diagnostic import *
from pyCMBS.plots import *
import scipy as sc
import pylab as pl
import numpy as np
from pyCMBS.region import *
####from pyCMBS import *

class TestData(TestCase):

    def setUp(self):
        # init Data object for testing
        n=1000  # slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100.  # generate dummy data
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
        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01') - 1
        self.D.time_str = "days since 0001-01-01 00:00:00"
        self.D.calendar = 'gregorian'

    def test_pattern_correlation(self):
        """
        test pattern correlation function
        """
        x = self.D.copy()

        # correlation with random values
        y = self.D.copy()
        tmp = np.random.random(y.shape)
        y.data = np.ma.array(tmp, mask=tmp != tmp)
        P2 = PatternCorrelation(x, y)
        P2._correlate()
        self.assertEqual(x.nt,len(P2.r_value))
        self.assertEqual(x.nt,len(P2.t))

        for i in xrange(x.nt):
            slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(x.data[i,:,:].flatten(),y.data[i,:,:].flatten())
            self.assertEqual(P2.r_value[i], r_value)
            self.assertEqual(P2.p_value[i], p_value)
            self.assertEqual(P2.slope[i], slope)
            self.assertEqual(P2.intercept[i], intercept)
            self.assertEqual(P2.std_err[i], std_err)



    def test_gleckler_index(self):
        """
        test Reichler index/Gleckler plot
        """

        # generate sample data
        # sample data
        tmp = np.zeros((5, 3, 1))
        tmp[:,0,0] = np.ones(5)*1.
        tmp[:,1,0] = np.ones(5)*2.
        tmp[:,2,0] = np.ones(5)*5.

        # The data is like ...
        #| 1 | 2 | 5 |
        #| 1 | 2 | 5 |
        #| 1 | 2 | 5 |
        #| 1 | 2 | 5 |
        #| 1 | 2 | 5 |

        x = self.D.copy()
        x._temporal_subsetting(0, 4)

        x.data = np.ma.array(tmp, mask=tmp!=tmp)
        x.std = np.ones(x.data.shape)
        x.time[0] = pl.datestr2num('2000-02-15')
        x.time[1] = pl.datestr2num('2000-03-15')
        x.time[2] = pl.datestr2num('2000-04-15')
        x.time[3] = pl.datestr2num('2000-05-15')
        x.time[4] = pl.datestr2num('2000-06-15')

        y = self.D.copy()
        y._temporal_subsetting(0, 4)
        tmp = np.ones(x.data.shape)  # sample data 2
        y.data = np.ma.array(tmp, mask=tmp!=tmp)
        y.time[0] = pl.datestr2num('2000-02-15')
        y.time[1] = pl.datestr2num('2000-03-15')
        y.time[2] = pl.datestr2num('2000-04-15')
        y.time[3] = pl.datestr2num('2000-05-15')
        y.time[4] = pl.datestr2num('2000-06-15')

        # Case 1: same area weights
        # cell area
        tmp = np.ones((3, 1))
        x.cell_area = tmp*1.

        #| 1-1 | 2-1 | 5-1 |
        #| 1-1 | 2-1 | 5-1 |
        #| 1-1 | 2-1 | 5-1 |
        #| 1-1 | 2-1 | 5-1 |
        #| 1-1 | 2-1 | 5-1 |
        #===================
        #| 0   | 5   | 5*4**2=5*16. = 80 |
        #==> E2 = sqrt(85./(15.))
        D = GlecklerPlot()
        r = D.calc_index(x, y, 'a', 'b', time_weighting=False)

        wt = np.ones(5) / 5.
        ref = np.sqrt(((85./15.) * wt).sum())
        t = np.abs(1. - r / ref)
        self.assertLess(t, 0.000001)  # relative error

        D = GlecklerPlot()
        r = D.calc_index(x, y, 'a', 'b')

        wt = np.asarray([29., 31., 30., 31., 30.])
        wt = wt / wt.sum()
        ref = np.sqrt(((85./15.) * wt).sum())
        t = np.abs(1. - r / ref)
        self.assertLess(t, 0.000001)  # relative error



        # Case 2: Different area weights
        # cell area
        tmp = np.ones((3, 1))
        tmp[1, 0] = 2.
        x.cell_area = tmp*1.

        #| 1-1=0 | 2-1=1 | 5-1=16 |
        #| 1-1=0 | 2-1=1 | 5-1=16 |
        #| 1-1=0 | 2-1=1 | 5-1=16 |
        #| 1-1=0 | 2-1=1 | 5-1=16 |
        #| 1-1=0 | 2-1=1 | 5-1=16 |
        #--------------------------
        # w = 0.25 w = 0.5  w=0.25|
        #--------------------------

        # 0.25*0 + 0.5 * 1 + 0.25 * 16 = 0 + 0.5 + 4 = 4.5
        # the mean of that is 4.5 for each timestep
        # mean because the overall weights are calculated as such that
        # they give a total weight if 1

        #  diagnostic
        D = GlecklerPlot()
        r = D.calc_index(x, y, 'a', 'b', time_weighting=False)

        wt = np.ones(5) / 5.
        ref = np.sqrt((4.5 * wt).sum())
        t = np.abs(1. - r / ref)
        self.assertLess(t, 0.000001)  # relative error

        wt = np.asarray([29., 31., 30., 31., 30.])
        wt = wt / wt.sum()
        ref = np.sqrt((4.5 * wt).sum())
        t = np.abs(1. - r / ref)
        self.assertLess(t, 0.000001)  # relative error

        # Case 3: use different std
        x.std = np.ones(x.data.shape)
        x.std[:, 2, 0] = 0.5

        #| 1-1=0 | 2-1=1  | 5-1=16 / 0.5 |
        #| 1-1=0 | 2-1=1  | 5-1=16 / 0.5 |
        #| 1-1=0 | 2-1=1  | 5-1=16 / 0.5 |
        #| 1-1=0 | 2-1=1  | 5-1=16 / 0.5 |
        #| 1-1=0 | 2-1=1  | 5-1=16 / 0.5 |
        #--------------------------------
        # w = 0.25    w = 0.5    w=0.25|
        #    0    +   0.5   +  0.25*32 = 0.5 + 8 = 8.5

        D = GlecklerPlot()
        r = D.calc_index(x, y, 'a', 'b', time_weighting=False)

        wt = np.ones(5) / 5.
        ref = np.sqrt((8.5 * wt).sum())
        t = np.abs(1. - r / ref)
        self.assertLess(t, 0.000001)  # relative error

        wt = np.asarray([29., 31., 30., 31., 30.])
        wt = wt / wt.sum()
        ref = np.sqrt((8.5 * wt).sum())
        t = np.abs(1. - r / ref)
        self.assertLess(t, 0.000001)  # relative error

