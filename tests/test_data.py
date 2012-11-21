from unittest import TestCase

__author__ = 'm300028'

#identify pyCMBS path and add it to pythonpath, as otherwise the modules are not found properly!

#import os
#p=os.path.split(os.getcwd())[0]
#print p
#os.environ['PYTHONPATH']=os.environ['PYTHONPATH']+':'+p+':'
#
#print os.environ['PYTHONPATH']

from data import *
#from diagnostic import *
import scipy as sc
import pylab as pl
import numpy as np
from scipy import stats

class TestData(TestCase):

    def setUp(self):
        #init Data object for testing
        n=1000 #slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100. #generate dummy data
        self.D = Data(None,None)
        d=np.ones((n,1,1))
        self.D.data = d
        self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data,mask=self.D.data != self.D.data)
        self.D.verbose = True
        self.D.unit = 'myunit'
        self.D.label = 'testlabel'

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

#    def test_get_percentile(self):
#        #print self.D.data
#        r = self.D.get_percentile(0.5,return_object = False)[0,0]
#        self.assertAlmostEqual(r,0.,delta = 0.5)

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

    def test_get_yearmean(self):
        #check get_yeartime
        D = self.D.copy()
        t1 = pl.datestr2num('2001-01-01') + np.arange(4) #year 2001
        t2 = pl.datestr2num('2005-05-15') + np.arange(4) #year 2005
        t3 = pl.datestr2num('2010-07-15') + np.arange(4) #year 2010
        D.time = np.asarray([t1,t2,t3]).flatten()
        data = pl.rand(len(D.time),1,1)
        data[8:,0,0] = np.nan
        D.data = np.ma.array(data,mask=np.isnan(data))       #generate random data
        r1 = np.mean(D.data[0:4]); r2 = np.mean(D.data[4:8]); r3=np.mean(D.data[8:])
        #print 'Reference results: ', r1, r2, r3
        years, res = D.get_yearmean()
        #print 'Result: ', res
        self.assertEqual(years[0],2001); self.assertEqual(years[1],2005)
        self.assertEqual(res[0,0,0],r1); self.assertEqual(res[1,0,0],r2)
        self.assertEqual(res[2,0,0].mask,r3.mask)

#    def test_diagnostic__get_valid_timeseries(self):
#        #test _get_valid_timeseries() of diagnostic tool
#        D = self.D.copy()
#
#        S = Diagnostic(D,D)
#        d,m = S._get_valid_timeseries(S.x)
#        print d
#        print m
#        stop

    def test_weighting_matrix(self):
        D = self.D.copy()
        D.cell_area = np.ones(D.data[0,:,:].shape)*0.333
        r = D._get_weighting_matrix()
        self.assertFalse(np.any(r != 1.))

    def test_adjust_time(self):
        D = self.D.copy()
        D.adjust_time(day=17)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).day,17)
        D.adjust_time(month=10)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).month,10)


    def test_diff(self):
        #test diff() function

        D = self.D.copy()

        D.data = np.ma.array(np.zeros((1000,1,2)),mask=np.zeros((1000,1,2)).astype('bool')  )
        D.data[:,0,0] = sc.randn(len(D.time))
        D.data[:,0,1] = sc.randn(len(D.time))

        A=D.copy()
        A.data[:,0,0] = sc.randn(len(D.time))
        A.data[:,0,1] = sc.randn(len(D.time))

        D.label='test2'

        x=D.data[:,0,0]; y=A.data[:,0,0]
        x1=D.data[:,0,1]; y1=A.data[:,0,1]
        t,p = stats.ttest_ind(x,y,axis=0)
        t1,p1 = stats.ttest_ind(x1,y1,axis=0)

        s  = A.diff(D,pthres=0.05)
        s1 = D.diff(D,pthres=0.05) #test with the same data

        #checks

#        print 'P: ', s.p_value, p, 1.-p
#        print s.p_value.shape
#        print s.p_value[0,0]
        self.assertAlmostEqual(s.p_value[0,0],1.-p,places=8)
        self.assertAlmostEqual(s.p_value[0,1],1.-p1,places=8)
        if p <= 0.05:
            self.assertEqual(s.p_mask[0,0],True)
        else:
            self.assertEqual(s.p_mask[0,0],False)

        #test for same data
        self.assertEqual(s1.p_value[0,0],0.)
        self.assertEqual(s1.p_value[0,1],0.)


        #another test of the t-test, taken from http://web.mst.edu/~psyworld/texample.htm
        x = np.asarray([5.,7.,5.,3.,5.,3.,3.,9.])
        y = np.asarray([8.,1.,4.,6.,6.,4.,1.,2.])

        A=self.D.copy(); B=self.D.copy()
        X = np.zeros((len(x),1,1)); Y = np.zeros((len(y),1,1))
        X[:,0,0] = x; Y[:,0,0] = y
        A.data = np.ma.array(X,mask=X!=X); B.data = np.ma.array(Y,mask=Y!=Y)

        u = A.diff(B,pthres=0.05)
        self.assertAlmostEqual(u.t_value[0,0],0.847,places=3)
        self.assertEqual(u.data[0,0],1.)









