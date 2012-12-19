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
import Nio

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
        self.D.filename = 'testinputfilename.nc'
        self.D.varname = 'testvarname'
        self.D.long_name = 'This is the longname'

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
        '''
        test the adjust_time function
        @return:
        '''
        '''
        @return:
        '''
        D = self.D.copy()
        D.adjust_time(day=17)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).day,17)
        D.adjust_time(month=10)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).month,10)
        D.adjust_time(year=2025)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).year,2025)

    def test_timsort(self):
        D=self.D.copy()
        D.adjust_time(day=15)

        #- generate some sample data
        D.time = pl.datestr2num('2001-05-03') + np.arange(5)
        D.data = D.data[0:5,:,:]
        D.data[:,0,0] = np.arange(5)

        D.std = D.data.copy()+2.2


        #- reshuffle the data
        t1=D.time[1]*1.; t2=D.time[3]*1.
        D.time[3] = t1; D.time[1] = t2

        #save reference solutions before sorting
        y = D.data[:,0,0]*1.
        t = D.time*1.

        s = np.argsort(t)
        y1 = y[s]

        print 'Time BEFORE sorting: ', D.time
        print 'Data BEFORE sorting: ', D.data[:,0,0]


        #sort data
        D.timsort()


        print 'Time AFTER sorting: ', D.time
        print 'Data AFTER sorting: ', D.data[:,0,0]
        print '                    ', y1

        #/// checks

        #a) check if time is sorted
        self.assertTrue(np.all(np.diff(D.time) > 0))

        #b) check if data was sorted also appropriately
        self.assertTrue(   np.all(y1-D.data[:,0,0]) == 0.   )
        self.assertTrue(   np.all(y1+2.2-D.std [:,0,0]) == 0.   )


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


    def test_save_netCDF(self):
        """
        test netCDF save routine
        """
        testfile='mytestfile.nc'
        self.D.save(testfile,varname='testvar',format='nc',delete=True)

        #read data again
        F = Data(testfile,'testvar',read=True,verbose=False)

        self.assertEqual(len(F.time),len(self.D.time))
        self.assertFalse(np.any(self.D.data-F.data) != 0. )

        del F

        #read data from default, this should then have the same variable name as self.D
        self.D.save(testfile,format='nc',delete=True)
        F = Data(testfile,'testvarname',read=True,verbose=False)

        self.assertEqual(len(F.time),len(self.D.time))
        self.assertFalse(np.any(self.D.data-F.data) != 0. )

        os.remove(testfile)


    def test_interp_time(self):
        D = self.D

        #time is from 2001-01-01 for 1000 days as default

        #case 1: interpolate to half daily values for a small timeperiod
        tref = pl.datestr2num('2001-05-05') + np.arange(200)*0.5+0.25
        #print 'Minimum/Maximum ' \
        #      ' date target: ', pl.num2date(tref.min()), pl.num2date(tref.max())
        #print 'Minimum/Maximum  date source: ', pl.num2date(D.time.min()), pl.num2date(D.time.max())

        #... interpolate data object for time period specified by tref
        I = D.interp_time(tref)

        #... original data
        y = D.data[:,0,0]

        #... generate reference solution using numpy
        yy = np.interp(tref, D.time, y)

        #... optional: plotting (good for validation of test routine)
        if False:
            pl.figure()
            pl.plot(D.time,y,color='blue')
            pl.plot(I.time,I.data[:,0,0],color='red',label='interpolated')
            pl.plot(tref,yy,color='green',label='reference interp',linestyle='--')
            pl.legend()
            pl.show()

        d = yy - I.data[:,0,0]
        self.assertFalse(np.any(np.abs(d[0:-1]) > 1.E-10 ) ) #boundary effects at end of period, therefore last value not used

    def test_div(self):
        #unittest for division

        D = self.D.copy()
        R = D.div(D)

        self.assertTrue(np.all(R.data == 1.))

    def test_divc(self):

        D = self.D.copy()
        R = D.divc(2.)

        d = D.data[:,0,0] *0.5
        self.assertTrue(np.all(d-R.data[:,0,0]) == 0.)

    def test_partial_correlation(self):
        """Test of partial correlation """

        x = self.D;
        nt,ny,nx = x.data.shape
        y = x.copy(); y.data = y.data + pl.rand(nt,ny,nx)*1000.
        z = x.copy(); z.data = z.data * pl.rand(nt,ny,nx)*100.

        res = x.partial_correlation(y,z)

        #generate reference solution
        slope, intercept, rxy, p_value, std_err = stats.linregress(x.data[:,0,0],y.data[:,0,0])
        slope, intercept, rxz, p_value, std_err = stats.linregress(x.data[:,0,0],z.data[:,0,0])
        slope, intercept, rzy, p_value, std_err = stats.linregress(z.data[:,0,0],y.data[:,0,0])

        ref = (rxy - rxz*rzy) / (np.sqrt(1.-rxz*rxz)*np.sqrt(1.-rzy*rzy))

        self.assertAlmostEqual(ref,res.data[0,0],places=5)


    def test__equal_lon(self):
        D=self.D

        #1) not equal longitudes
        D.lon = pl.rand(100,200)
        self.assertFalse(D._equal_lon())

        #2) equal longitudes
        x=np.arange(100)
        D.lon = np.zeros((2,100))
        D.lon[0,:] = x
        D.lon[1,:] = x
        self.assertTrue(D._equal_lon())

    def test__get_unique_lon(self):
        D = self.D
        #equal longitudes
        x=np.arange(100)
        D.lon = np.zeros((2,100))
        D.lon[0,:] = x; D.lon[1,:] = x

        r = D._get_unique_lon()
        self.assertTrue(np.all((x-r) == 0.))



    def generate_tuple(self,n=None,mask=True):
        #generate perturbed tuple of data
        x = self.D.copy(); y = self.D.copy()
        nt,ny,nx = x.data.shape
        z = pl.randn(nt,ny,nx)
        y.data = y.data*z
        if mask:
            y.data = np.ma.array(y.data,mask=z>0.5) #mask some data so we have data with different masks
        else:
            y.data = np.ma.array(y.data,mask=y.data != y.data) #mask some data so we have data with different masks

        if n != None:
            if n < len(x.data)-1:
                x._temporal_subsetting(0,n)
                y._temporal_subsetting(0,n)

        return x,y


    #-----------------------------------------------------------

    def test_correlate(self):
        #test correlation

        for n in [None,100,10,5]: #different size

            x,y = self.generate_tuple(n=n,mask=True)
            x1=x.data[:,0,0]; y1=y.data[:,0,0]
            msk = (x1.mask == False) & (y1.mask == False)
            x2 = x1[msk]; y2 = y1[msk] #this is only the valid data

            print 'Number of masked pixels: ', sum(y.data.mask), n

            ##################################################################
            # PEARSON CORRELATION
            ##################################################################
            slope, intercept, r_value1, p_value1, std_err = stats.mstats.linregress(x1,y1) #masked
            slope, intercept, r_value2, p_value2, std_err = stats.linregress(x2,y2) #not masked
            r,p = x.correlate(y)

            #1) test if scipy functions return similar results
            self.assertAlmostEqual(r_value1,r_value2,places=15)
            #self.assertAlmostEqual(p_value1,p_value2,places=15) #not used, as BUG in stats.mstats.linregress!

            #2) test data.correlate() results
            self.assertAlmostEqual(r.data[0,0],r_value2,places=10) #results from stats.linregress are used, as mstats is BUGGY!!
            self.assertAlmostEqual(p.data[0,0],p_value2,places=10)


            ##################################################################
            # SPEARMAN RANK CORRELATION
            ##################################################################

            # 1) test if scipy functions return similar results for masked/not masked arrays
            r_value1, p_value1 = stats.mstats.spearmanr(x1,y1) #masked
            r_value2, p_value2 = stats.spearmanr(x2,y2) #not masked

            self.assertAlmostEqual(r_value1,r_value2,places=10)
            self.assertAlmostEqual(p_value1,p_value2,places=10)

            #2) test data.correlate() function
            r,p = x.correlate(y,spearman=True)
            self.assertAlmostEqual(r.data[0,0],r_value1,places=10)
            self.assertAlmostEqual(p.data[0,0],p_value1,places=10)
            self.assertAlmostEqual(r.data[0,0],r_value2,places=10)
            self.assertAlmostEqual(p.data[0,0],p_value2,places=10)

























