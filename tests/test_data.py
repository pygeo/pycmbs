from unittest import TestCase

__author__ = 'm300028'

from data import *
import scipy as sc
import pylab as pl
import numpy as np
from scipy import stats

class TestData(TestCase):

    def setUp(self):
        #init Data object for testing
        n=1000000
        x = sc.randn(n)*100. #generate dummy data
        self.D = Data(None,None)
        d=np.ones((n,1,1))
        self.D.data = d
        self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data,mask=self.D.data != self.D.data)
        self.D.verbose = True
        self.D.unit = 'myunit'

        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01') - 1

    def test_get_time_indices(self):
        d1 = pl.num2date(pl.datestr2num('2001-01-05'))
        d2 = pl.num2date(pl.datestr2num('2001-05-05'))
        i1,i2 = self.D._get_time_indices(d1,d2)
        s1 = str(pl.num2date(self.D.time[i1]))
        s2 = str(pl.num2date(self.D.time[i2]))

        print s1, i1
        print s2, i2

        self.assertEqual(s1,'2001-01-05 00:00:00+00:00')
        self.assertEqual(s2,'2001-05-05 00:00:00+00:00')

    def test_addc(self):
        #test with data copy
        r1 = self.D.addc(5.,copy=True)
        self.assertEqual(r1.data[4,0,0]-5.,self.D.data[4,0,0])
        #test without data copy
        ref = self.D.data[5,0,0]
        self.D.addc(666.,copy=False)
        self.assertEqual(ref+666.,self.D.data[5,0,0])

    def test_get_percentile(self):
        #print self.D.data
        r = self.D.get_percentile(0.5,return_object = False)[0,0]
        self.assertAlmostEqual(r,0.,delta = 0.5)

    def test_correlate(self):
        #test for correlation calculations
        r,p = self.D.correlate(self.D,pthres=1.01) #1) correlation with itself (returns data objects)
        self.assertEqual(r.data[0,0],1.)
        self.assertEqual(p.data[0,0],0.)

    def test_set_time(self):
        #self.D.time_str = 'day as %Y%m%d.%f'
        #self.D.time = [20120503,19750101]
        #self.D.set_time()
        #print self.time
        self.D.time_str = "days since 1992-01-01 00:00:00"
        self.D.time = np.array([6633, 6725, 6817, 6908])
        self.D.set_time()
        if self.D.verbose:
            print self.D.time[0]
            print self.D.time[3]
        self.assertEqual(self.D.time[0],733831.)
        self.assertEqual(self.D.time[3],734106.)

        print ''
        print '---'
        self.D.time_str = "days since 0001-01-01 00:00:00" #behaviour look a bit strange here, as num2date gives number of days PLUS one (see num2date docstring)
        self.D.time = np.array([0.])
        self.D.set_time()
        if self.D.verbose:
            print 'Results: ', self.D.time[0], pl.num2date(self.D.time[0])
        self.assertEqual(self.D.time[0],1.)

        print ''
        print '---'
        self.D.time_str = "days since 1997-10-21 00:00:00" #o.k.
        self.D.time = np.array([0., 2.])
        self.D.set_time()
        if self.D.verbose:
            print self.D.time[0],pl.num2date(self.D.time[0])
            print self.D.time[1],pl.num2date(self.D.time[1])
        self.assertEqual(self.D.time[0],729318.0)
        self.assertEqual(self.D.time[1],729320.0)


    def test_temporal_trend(self):
        y = np.arange(len(self.D.time))*2.+8.
        self.D.data[:,0,0] = y

        #reference solution
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.D.time,y)

        #calculate temporal correlation WITHOUT normalization of time
        R,S,I,P = self.D.temporal_trend() #no object is returned (default)
        self.assertEqual(R[0,0],r_value); self.assertEqual(S[0,0],slope)
        self.assertEqual(I[0,0],intercept)






