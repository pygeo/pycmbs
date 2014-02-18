from unittest import TestCase
import unittest

from pycmbs.data import Data
import os
import scipy as sc
import pylab as pl
import numpy as np
from scipy import stats
from dateutil.rrule import rrule
from dateutil.rrule import MONTHLY
import pylab as pl

from nose.tools import assert_raises



class TestData(unittest.TestCase):

    def setUp(self):
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
        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01')
        self.D.time_str = "days since 0001-01-01 00:00:00"
        self.D.calendar = 'gregorian'
        self.D.oldtime=False

    def test_DataInitLabelNotNone(self):
        d = Data(None,None, label='testlabel')
        self.assertEqual(d.label, 'testlabel')

    def test_DataInitUnitNotNone(self):
        d = Data(None,None, unit='myfunkyunit')
        self.assertEqual(d.unit, 'myfunkyunit')

    def test_DataInitTimeCycleNotNone(self):
        d = Data(None,None, time_cycle=24)
        self.assertEqual(d.time_cycle, 24)

    def test_get_time_indices(self):
        d1 = pl.num2date(pl.datestr2num('2001-01-05'))
        d2 = pl.num2date(pl.datestr2num('2001-05-05'))
        self.D._oldtime = True
        i1,i2 = self.D._get_time_indices(d1,d2)
        s1 = str(pl.num2date(self.D.time[i1]))
        s2 = str(pl.num2date(self.D.time[i2]))
        self.assertEqual(s1,'2001-01-05 00:00:00+00:00')
        self.assertEqual(s2,'2001-05-05 00:00:00+00:00')

    def test_get_time_indices_InvalidDates(self):
        i1, i2 = self.D._get_time_indices(None, None)
        self.assertEqual(i1, 0)
        self.assertEqual(i2, len(self.D.time)-1)

    def test_sub_sample_InvalidGeometry(self):
        x = self.D.copy()
        x.data = np.random.random((3,4,5,6))
        with self.assertRaises(ValueError):
            x._sub_sample(3)

    def test_squeeze(self):
        self.D.squeezed = False
        self.D._squeeze()
        self.assertTrue(self.D.squeezed)



    def test_sub_sample(self):
        x = self.D.copy()

        # 3D data
        nt_org = self.D.nt
        tmp = np.random.random((nt_org, 50, 80))
        x.data = np.ma.array(tmp, mask=tmp != tmp)
        x._sub_sample(10)
        nt,ny,nx = x.shape
        # check geometry of results first
        self.assertEqual(nt, nt_org)
        self.assertEqual(ny, 5)
        self.assertEqual(nx, 8)
        # now check values
        self.assertEqual(x.data[0,0,0], tmp[0,0,0])
        self.assertEqual(x.data[172,0,0], tmp[172,0,0])

        #todo continue here with checks of values!
        #~ print x.data[10,0,0]
        #~ print tmp[10,8:12,8:12]
        #~ self.assertEqual(x.data[10,0,0], tmp[10,11,11])

        # 2D
        tmp = np.random.random((73, 92))
        x.data = np.ma.array(tmp, mask=tmp != tmp)
        x._sub_sample(5)
        # check geometry of results first
        ny,nx = x.shape
        self.assertEqual(ny, 14+1)
        self.assertEqual(nx, 18+1)


    def test_oldtimeoffset_Invalid(self):
        d = self.D.copy()
        del d.time_str
        with self.assertRaises(ValueError):
            d._oldtimeoffset()

    def test_oldtimeoffset_InvalidTimeStr(self):
        d = self.D.copy()
        d.time_str = 'no_time_str'
        with self.assertRaises(ValueError):
            d._oldtimeoffset()

    def test_oldtimeoffset_Invalid(self):
        d = self.D.copy()
        d.time_str = 'hours'
        self.assertEqual(d._oldtimeoffset(), 24.)
        d.time_str = 'seconds'
        self.assertEqual(d._oldtimeoffset(), 86400.)
        d.time_str = 'days'
        self.assertEqual(d._oldtimeoffset(), 1.)



    def test_get_temporal_mask(self):
        x = self.D.copy()

        # test monthly mask
        mm = x.get_temporal_mask([1,5,11],mtype='monthly')
        d = x.date[mm]
        for t in d:
            self.assertTrue(t.month in [1,5,11])

        # yearly mask
        ym = x.get_temporal_mask([2002],mtype='yearly')
        d = x.date[ym]
        for t in d:
            self.assertTrue(t.year in [2002])
        ym = x.get_temporal_mask([2003],mtype='yearly')
        d = x.date[ym]
        for t in d:
            self.assertTrue(t.year in [2003])

    def test__get_date_from_month(self):
        x = self.D.copy()
        x.time_str = 'months since 1983-05-01 00:00:00'
        d1 = x._get_date_from_month(2)
        self.assertEqual(d1.year,1983)
        self.assertEqual(d1.month,7)
        self.assertEqual(d1.day,1)

        x.time_str = 'months since 1987-07-13 00:00:00'
        d1 = x._get_date_from_month(5)
        self.assertEqual(d1.year,1987)
        self.assertEqual(d1.month,12)
        self.assertEqual(d1.day,13)

        x.time_str = 'months since 1987-08-22 00:00:00'
        d1 = x._get_date_from_month(9)
        self.assertEqual(d1.year,1988)
        self.assertEqual(d1.month,5)
        self.assertEqual(d1.day,22)


    def test_get_climatology(self):
        x = self.D.copy()

        # timecycle = 1
        x.time_cycle = 1
        r = x.data.mean(axis=0)
        c = x.get_climatology()
        d = np.abs(1.-r/c)
        self.assertTrue(np.all(d < 1.E-6))

        # ... same, but with object returned
        c = x.get_climatology(return_object=True)
        d = np.abs(1.-r/c.data)
        self.assertTrue(np.all(d < 1.E-6))

        # varying timecycles
        for time_cycle in [1,5,12,23]:
            x.time_cycle=time_cycle
            c = x.get_climatology()
            nt,ny,nx = x.shape
            r = np.zeros((time_cycle,ny,nx))
            n = np.zeros((time_cycle,ny,nx))
            cnt = 0
            for i in xrange(nt):
                if cnt % time_cycle == 0:
                    cnt = 0
                r[cnt,:,:] = r[cnt,:,:] + x.data[i,:,:]
                n[cnt,:,:] = n[cnt,:,:] + (~x.data.mask[i,:,:]).astype('int')
                cnt +=1
            res = r / n  # reference mean
            d = np.abs(1.-res/c)
            self.assertTrue(np.all(d < 1.E-6))


    def test_set_valid_range(self):
        x = self.D.copy()
        tmp = np.random.random((100,200,300)) * 10. - 5.
        x.data = np.ma.array(tmp, mask=tmp != tmp)
        x._set_valid_range(-2., 2.)
        self.assertTrue(np.all(x.data >=-2.))
        self.assertTrue(np.all(x.data <=2.))

        x._set_valid_range(-0.5, 1.)
        self.assertTrue(np.all(x.data >=-0.5))
        self.assertTrue(np.all(x.data <=1.))

    def test_is_monthly(self):
        a = self.D.copy()
        b = self.D.copy()
        t=[]
        x = pl.datestr2num('2001-01-15')
        for i in xrange(20):
            t.append(x)
            x += 30
        t = np.asarray(t)
        b.time = t
        self.assertEqual(a._is_monthly(), False)
        self.assertEqual(b._is_monthly(), True)

    def test_is_monthly_MissingTime(self):
        x = self.D.copy()
        del x.time
        self.assertFalse(x._is_monthly())


    def test_detrend_InvalidGeometry(self):
        x = self.D.copy()
        x.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            x.detrend()

    def test_cut_bounding_box(self):
        x = self.D.copy()
        # sample data with invalid boundaries
        t = np.random.random((6,6,5))
        #... left border 1 pix
        t[:,:,0] = np.nan
        #... top border 2pix
        t[:,0,:] = np.nan
        t[:,1,:] = np.nan
        #... right border only some pixels invalid
        t[:,0:4,-1] = np.nan

        x.data = np.ma.array(t, mask = np.isnan(t))
        y = x.cut_bounding_box(return_object=True)

        # left border
        self.assertEqual(y.data[0,0,0], x.data[0,2,1])
        self.assertEqual(y.data[2,0,0], x.data[2,2,1])

        # right border
        #todo
        #~ print x.data[0,2,:]
        #~ print y.data[0,0,:]
        #~ self.assertEqual(y.data[0,0,-2], x.data[0,2,-2])


    def test_add(self):
        x = self.D.copy()
        y = self.D.copy()
        y.data += 3.
        c = x.add(y)
        self.assertEqual(c.data[0,0,0], x.data[0,0,0]*2.+3.)
        self.assertEqual(c.data[100,0,0], x.data[100,0,0]*2.+3.)

    def test_sub(self):
        x = self.D.copy()
        y = self.D.copy()
        y.data += 3.
        c = x.sub(y)
        self.assertTrue(np.abs(1.-c.data[0,0,0]/-3.) < 1.E-6)
        self.assertTrue(np.abs(1.-c.data[100,0,0]/-3.) < 1.E-6)

    def test_addc(self):
        r1 = self.D.addc(5.,copy=True)
        self.assertAlmostEqual(r1.data[4,0,0]-5.,self.D.data[4,0,0], 8)

    def testAddcWithoutDataCopy(self):
        ref = self.D.data[5,0,0]
        self.D.addc(666.,copy=False)
        self.assertEqual(ref+666.,self.D.data[5,0,0])

    def test_get_percentile(self):
        for p in [0.05, 0.5, 0.95]:
            r = self.D.get_percentile(p, return_object = False)[0,0]
            res = stats.mstats.scoreatpercentile(self.D.data[:,0,0], p * 100.)
            self.assertAlmostEqual(r, res)

            r = self.D.get_percentile(p, return_object = True)
            self.assertAlmostEqual(r.data[0,0], res)

    def test_timn(self):
        A = self.D.copy()
        B = self.D.copy()
        me = A.timmean()
        su = B.timsum()
        an = A.timn()  # ndarray
        bn = B.timn(return_object=True)  # Data object
        r = su/me
        self.assertEqual(r[0,0], an[0,0])
        self.assertEqual(r[0,0], bn.data[0,0])

    def test_flipud(self):
        x = self.D.copy()
        y = x.copy()
        y._flipud()
        self.assertEqual(x.data[0,0,0], y.data[0,-1,0])
        self.assertEqual(x.data[0,-1,0], y.data[0,0,0])


    def test_flipud_InvalidGeometry(self):
        x = self.D.copy()
        x.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            x._flipud()
        x.data = np.random.random((10,))
        with self.assertRaises(ValueError):
            x._flipud()


    @unittest.skip('wait for bugfree scipy')
    def test_correlate1(self):
        #test for correlation calculations
        r,p = self.D.correlate(self.D, pthres=1.01) #1) correlation with itself (returns data objects)
        self.assertEqual(r.data[0,0], 1.)
        self.assertEqual(p.data[0,0], 0.)

    def test_correlate_WithInvalidGeometries(self):
        x = self.D.copy()
        y = self.D.copy()
        y.data = np.random.random((10, 200, 273))
        with self.assertRaises(ValueError):
            x.correlate(y)

    def get_date_from_months_WithInvalidTimestr(self):
        d = self.D
        d.time_str = 'some_invalid_str'
        with self.assertRaises(ValueError):
            d._get_date_form_month(10)

    def test_set_time(self):
        # NB: num2date gives number of days PLUS one (see num2date docstring)
        self.D.time_str = "days since 0001-01-01 00:00:00"
        self.D.time = np.array([1.])
        self.D.set_time()
        self.assertEqual(self.D.time[0],1.)

    def testTemporalTrendNoTimeNormalization(self):
        y = np.arange(len(self.D.time))*2.+8.
        self.D.data[:, 0, 0] = y

        # reference solution
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.D.time,y)

        # calculate temporal correlation WITHOUT normalization of time
        R, S, I, P = self.D.temporal_trend()  # no object is returned (default)
        self.assertEqual(R[0,0], r_value)
        self.assertEqual(S[0,0], slope)

        R, S, I, P = self.D.temporal_trend(return_object=True)  # TODO further tests for slope and significance
        self.assertEqual(R.data[0, 0], r_value)
        self.assertEqual(S.data[0,0], slope)

    def test_timmean_InvalidDimension(self):
        with self.assertRaises(ValueError):
            d = self.D.copy()
            d.data = np.random.random((10, 20, 30, 40))
            d.timmean()

    def test_timstd_InvalidDimension(self):
        with self.assertRaises(ValueError):
            d = self.D.copy()
            d.data = np.random.random((10, 20, 30, 40))
            d.timstd()

    def test_timstd_2D(self):
        d = self.D.copy()
        d.data = np.random.random((10, 20))
        r = d.timstd()
        self.assertEqual(r, None)

    def test_timvar_2D(self):
        d = self.D.copy()
        d.data = np.random.random((10, 20))
        r = d.timvar()
        self.assertEqual(r, None)

    def test_timstd_2D(self):
        d = self.D.copy()
        d.data = np.random.random((10, 20))
        r = d.timstd()
        self.assertEqual(r, None)

    def test_timmean_timvar_consistency(self):
        d = self.D.copy()
        s = d.timstd(return_object=False)
        v = d.timvar(return_object=False)
        r = np.abs(1.- v / (s*s))
        self.assertTrue(np.all(r < 1.E-6))

    def test_timmin_InvalidDimension(self):
        with self.assertRaises(ValueError):
            d = self.D.copy()
            d.data = np.random.random((10, 20, 30, 40))
            d.timmin()

    def test_timmin_2D(self):
        d = self.D.copy()
        d.data = d.data[0,:,:]
        r= d.timmin(return_object=False)
        self.assertEqual(self.D.data[0,0,0], r[0,0])

    def test_timmax_InvalidDimension(self):
        with self.assertRaises(ValueError):
            d = self.D.copy()
            d.data = np.random.random((10, 20, 30, 40))
            d.timmax()

    def test_timmax_2D(self):
        d = self.D.copy()
        d.data = d.data[0,:,:]
        r= d.timmax(return_object=False)
        self.assertEqual(self.D.data[0,0,0], r[0,0])



    def test_get_yearmean(self):
        #check get_yeartime
        D = self.D.copy()
        t1 = pl.datestr2num('2001-01-01') + np.arange(4)
        t2 = pl.datestr2num('2005-05-15') + np.arange(4)
        t3 = pl.datestr2num('2010-07-15') + np.arange(4)
        D.time = np.asarray([t1,t2,t3]).flatten()
        D._oldtime = True
        data = pl.rand(len(D.time), 1, 1)
        data[8:, 0, 0] = np.nan
        D.data = np.ma.array(data,mask=np.isnan(data))
        r1 = np.mean(D.data[0:4])
        r2 = np.mean(D.data[4:8])
        r3=np.mean(D.data[8:])

        years, res = D.get_yearmean()

        self.assertEqual(years[0],2001)
        self.assertEqual(years[1],2005)
        self.assertEqual(res[0,0,0],r1)
        self.assertEqual(res[1,0,0],r2)
        self.assertEqual(res[2,0,0].mask,r3.mask)

        R = D.get_yearmean(return_data=True)
        self.assertEqual(R.date[0].year, 2001)
        self.assertEqual(R.date[1].year, 2005)
        self.assertEqual(R.data[0,0,0], r1)
        self.assertEqual(R.data[1,0,0], r2)
        self.assertEqual(R.data[2,0,0].mask, r3.mask)


    def test_get_yearsum(self):
        #check get_yeartime
        D = self.D.copy()
        t1 = pl.datestr2num('2001-01-01') + np.arange(4) #year 2001
        t2 = pl.datestr2num('2005-05-15') + np.arange(4) #year 2005
        t3 = pl.datestr2num('2010-07-15') + np.arange(4) #year 2010
        D.time = np.asarray([t1,t2,t3]).flatten()
        D._oldtime = True #use old python pylab time definition to be compliant with the test results here
        data = pl.rand(len(D.time), 1, 1)
        data[8:, 0, 0] = np.nan
        D.data = np.ma.array(data,mask=np.isnan(data))       #generate random data
        r1 = np.sum(D.data[0:4])
        r2 = np.sum(D.data[4:8])
        r3 = np.sum(D.data[8:])
        years, res = D.get_yearsum()
        resobj = D.get_yearsum(return_data=True)

        self.assertEqual(years[0],2001)
        self.assertEqual(resobj.date[0].year,2001)

        self.assertEqual(years[1],2005)
        self.assertEqual(resobj.date[1].year,2005)
        self.assertEqual(res[0,0,0],r1)
        self.assertEqual(resobj.data[0,0,0],r1)
        self.assertEqual(res[1,0,0],r2)
        self.assertEqual(resobj.data[1,0,0],r2)
        #~ self.assertEqual(res[2,0,0].mask,r3)





#    def test_diagnostic__get_valid_timeseries(self):
#        #test _get_valid_timeseries() of diagnostic tool
#        D = self.D.copy()
#
#        S = Diagnostic(D,D)
#        d,m = S._get_valid_timeseries(S.x)
#        print d
#        print m
#        stop


    def test_weighting_matrix_InvalidType(self):
        d = self.D.copy()
        d.weighting_type = 'invalid_value'
        with self.assertRaises(ValueError):
            d._get_weighting_matrix()

    def test_weighting_matrix(self):
        D = self.D.copy()  # single pixel

        x = np.ones((10,2,1))
        D.data=np.ma.array(x,mask=x == 0.)

        # case 1: valid data for all timestep
        D.cell_area = np.ones(D.data[0,:,:].shape)
        D.cell_area[0,0] = 75.; D.cell_area[1,0] = 25. #3/4 ; 1/4
        r = D._get_weighting_matrix()
        self.assertFalse(np.any(r[:,0,0] != 0.75))
        self.assertFalse(np.any(r[:,1,0] != 0.25))

        # case 2: invalid data for some timesteps
        D.data.mask[0,0,0] = True #mask one data as invalid
        r = D._get_weighting_matrix()
        self.assertFalse(np.any(r[1:,0,0] != 0.75))
        self.assertFalse(np.any(r[1:,1,0] != 0.25))
        self.assertFalse(r[0,1,0] != 1.)
        self.assertFalse(r.mask[0,0,0] == False)

        #case 3: invalid data, but normalization for whole area!
        D.weighting_type = 'all'
        r = D._get_weighting_matrix()
        self.assertFalse(r[0,1,0] != 0.25)

    def test_adjust_time(self):
        D = self.D.copy()
        D._oldtime = True #use old time convention to be compliant with test routines here
        D.adjust_time(day=17)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).day, 17)
        D.adjust_time(month=10)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).month, 10)
        D.adjust_time(year=2025)
        for i in xrange(len(D.time)):
            self.assertEqual(pl.num2date(D.time[i]).year, 2025)

    def test_timstat(self):
        """
        test temporal statistic functions
        @return:
        """
        D = self.D.copy()

        me = D.data.mean(axis=0)
        ME = D.timmean(return_object=True)
        self.assertEquals(me[0],ME.data[0])

        su = D.data.sum(axis=0)
        SU = D.timsum(return_object=True)
        self.assertEquals(su[0],SU.data[0])

        st = D.data.std(axis=0)
        ST = D.timstd(return_object=True)
        self.assertEquals(st[0],ST.data[0])

        cv = st/me
        CV = D.timcv(return_object=True)
        self.assertEquals(cv[0],CV.data[0])

        cv = st/me
        CV = D.timcv(return_object=False)
        self.assertEquals(cv[0],CV[0])

        va = D.data.var(axis=0)
        VA = D.timvar(return_object=True)
        self.assertEquals(va[0],VA.data[0])

        mi = D.data.min(axis=0)
        MI = D.timmin(return_object=True)
        self.assertEquals(mi[0],MI.data[0])

        ma = D.data.max(axis=0)
        MA = D.timmax(return_object=True)
        self.assertEquals(ma[0],MA.data[0])


    def test_get_years(self):
        d = self.D.date
        y = self.D._get_years()
        for i in xrange(self.D.nt):
            self.assertEqual(d[i].year, y[i])

    def test_get_months(self):
        d = self.D.date
        y = self.D._get_months()
        for i in xrange(self.D.nt):
            self.assertEqual(d[i].month, y[i])


    def test_days_per_month(self):
        ref = {1:[31],2:[28,29],3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
        x = self.D.copy()
        days = x._days_per_month()
        for i in xrange(x.nt):
            d = x.date[i]
            if d.month == 2:
                if d.year % 4 == 0:
                    self.assertEqual(days[i],29)
                else:
                    self.assertEqual(days[i],28)
            else:
                self.assertEqual(ref[d.month], days[i])

    def test_get_dateboundaries(self):
        # check mindate/maxdate functions
        x = self.D.copy()
        ma_date = x.date.max()
        mi_date = x.date.min()

        self.assertEqual(x._get_maxdate(), ma_date)
        self.assertEqual(x._get_mindate(), mi_date)

        self.assertEqual(x._get_maxdate(base='day').hour, 23)
        self.assertEqual(x._get_maxdate(base='day').minute, 59)
        self.assertEqual(x._get_maxdate(base='day').second, 59)

        self.assertEqual(x._get_mindate(base='day').hour, 0)
        self.assertEqual(x._get_mindate(base='day').minute, 0)
        self.assertEqual(x._get_mindate(base='day').second, 0)

        self.assertEqual(x._get_maxdate(base='month').hour, 23)
        self.assertEqual(x._get_maxdate(base='month').minute, 59)
        self.assertEqual(x._get_maxdate(base='month').second, 59)
        #~ self.assertEqual(x._get_maxdate(base='month').day, 1)

        self.assertEqual(x._get_mindate(base='month').hour, 0)
        self.assertEqual(x._get_mindate(base='month').minute, 0)
        self.assertEqual(x._get_mindate(base='month').second, 0)
        self.assertEqual(x._get_mindate(base='month').day, 1)

        self.assertEqual(x._get_maxdate(base='year').hour, 23)
        self.assertEqual(x._get_maxdate(base='year').minute, 59)
        self.assertEqual(x._get_maxdate(base='year').second, 59)
        self.assertEqual(x._get_maxdate(base='year').day, 31)
        self.assertEqual(x._get_maxdate(base='year').month, 12)

        self.assertEqual(x._get_mindate(base='year').hour, 0)
        self.assertEqual(x._get_mindate(base='year').minute, 0)
        self.assertEqual(x._get_mindate(base='year').second, 0)
        self.assertEqual(x._get_mindate(base='year').day, 1)
        self.assertEqual(x._get_mindate(base='year').month, 1)






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

        #print 'Time BEFORE sorting: ', D.time
        #print 'Data BEFORE sorting: ', D.data[:,0,0]


        #sort data
        D.timsort()


        #print 'Time AFTER sorting: ', D.time
        #print 'Data AFTER sorting: ', D.data[:,0,0]
        #print '                    ', y1

        #/// checks

        #a) check if time is sorted
        self.assertTrue(np.all(np.diff(D.time) > 0))

        #b) check if data was sorted also appropriately
        self.assertTrue(   np.all(y1-D.data[:,0,0]) == 0.   )
        self.assertTrue(   np.all(y1+2.2-D.std [:,0,0]) == 0.   )


    @unittest.skip('wait for bugfree scipy')
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

        x=D.data[:,0,0]
        y=A.data[:,0,0]
        x1=D.data[:,0,1]
        y1=A.data[:,0,1]
        t,p = stats.ttest_ind(x,y,axis=0)
        t1,p1 = stats.ttest_ind(x1,y1,axis=0)

        s = A.diff(D, pthres=0.05)
        s1 = D.diff(D, pthres=0.05)  # test with the same data

        #checks
        self.assertAlmostEqual(s.p_value[0,0], 1.-p, places=8)
        self.assertAlmostEqual(s.p_value[0,1], 1.-p1, places=8)
        if p <= 0.05:
            self.assertEqual(s.p_mask[0,0], True)
        else:
            self.assertEqual(s.p_mask[0,0], False)

        #test for same data
        self.assertEqual(s1.p_value[0,0], 0.)
        self.assertEqual(s1.p_value[0,1], 0.)


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
        testfile = './mytestfile.nc'
        self.D.save(testfile, varname='testvar', format='nc', delete=True)

        # read data again
        F = Data(testfile, 'testvar', read=True, verbose=False)

        self.assertEqual(len(F.time),len(self.D.time))
        self.assertFalse(np.any(self.D.data-F.data) != 0. )
        self.assertFalse(np.any(self.D.time-F.time) != 0. )

        del F

        # read data from default, this should then have the same variable name as self.D
        self.D.save(testfile, format='nc', delete=True)
        F = Data(testfile, 'testvarname', read=True, verbose=False)

        self.assertEqual(len(F.time), len(self.D.time))
        self.assertFalse(np.any(self.D.data-F.data) != 0. )

        os.remove(testfile)


    def test_interp_time_InvalidMethod(self):
        tref = self.D.num2date(pl.datestr2num('2001-05-05') + np.arange(200)*0.5+0.25)
        with self.assertRaises(ValueError):
            self.D.interp_time(tref, method='invalid_method')

    def test_interp_time_InvalidGeometry(self):
        d = self.D.copy()
        d.data = np.random.random((10,20,30,40))
        tref = d.num2date(pl.datestr2num('2001-05-05') + np.arange(200)*0.5+0.25)
        with self.assertRaises(ValueError):
            d.interp_time(tref)

    def test_interp_time_TimeNotAscending(self):
        d = self.D.copy()
        d.time = np.random.random(d.nt)
        tref = d.num2date(pl.datestr2num('2001-05-05') + np.arange(200)*0.5+0.25)
        with self.assertRaises(ValueError):
            d.interp_time(tref)

    def test_interp_time_TrefNotAscending(self):
        d = self.D.copy()
        tref = d.num2date(np.random.random(200))
        with self.assertRaises(ValueError):
            d.interp_time(tref)

    def test_interp_time(self):
        D = self.D.copy()
        import datetime

        start_date = datetime.datetime(2001,6,5)
        stop_date = datetime.datetime(2001,7,31)
        D.apply_temporal_subsetting(start_date, stop_date)

        #time is from 2001-01-01 for 1000 days as default

        #case 1: interpolate to half daily values for a small timeperiod
        tref = D.num2date(pl.datestr2num('2001-07-05') + np.arange(20)*0.5+0.25)
        # 5.July ... 14.July
        #~ print tref

        #... interpolate data object for time period specified by tref
        I = D.interp_time(tref)

        #... original data
        y = D.data[:, 0, 0]
        #... generate reference solution using numpy
        yy = np.interp(D.date2num(tref), D.time, y)

        #... optional: plotting (good for validation of test routine)
        if True:
            pl.figure()
            pl.plot(D.date, y, color='blue', label='original data')
            pl.plot(I.date, I.data[:,0,0], color='red', label='interpolated')
            pl.plot(tref, yy, color='green',label='reference interp',linestyle='--')
            pl.legend()
            pl.show()

        d = yy - I.data[:, 0, 0]
        self.assertFalse(np.any(np.abs(d[0:-1]) > 1.E-10 ) )  # boundary effects at end of period, therefore last value not used


    def test_date2num_NoTimeStr(self):
        del self.D.time_str
        t = np.arange(10).astype('float')
        with self.assertRaises(ValueError):
            self.D.date2num(t)

    def test_date2num_InvalidTimeStr(self):
        self.D.time_str=None
        t = np.arange(10).astype('float')
        with self.assertRaises(ValueError):
            self.D.date2num(t)

    def test_num2date_NoTimeStr(self):
        del self.D.time_str
        t = np.arange(10).astype('float')
        with self.assertRaises(ValueError):
            self.D.num2date(t)

    def test_num2date_InvalidTimeStr(self):
        self.D.time_str=None
        t = np.arange(10).astype('float')
        with self.assertRaises(ValueError):
            self.D.num2date(t)

    def test_save_ascii(self):
        self.D._save_ascii('testexport.txt', delete=True)
        self.assertTrue(os.path.exists('testexport.txt'))
        os.remove('testexport.txt')

    def test_save_ascii_FileExistingAlreadyDelete(self):
        if not os.path.exists('testexport.txt'):
            os.system('touch testexport.txt')
        self.D._save_ascii('testexport.txt', delete=True)
        self.assertTrue(os.path.exists('testexport.txt'))
        os.remove('testexport.txt')

    def test_save_ascii_FileExistingAlreadyNoDelete(self):
        if not os.path.exists('testexport.txt'):
            os.system('touch testexport.txt')
        with self.assertRaises(ValueError):
            self.D._save_ascii('testexport.txt', delete=False)
        os.remove('testexport.txt')

    def test_div_Default(self):
        D = self.D.copy()
        R = D.div(D)
        self.assertTrue(np.all(R.data == 1.))

    def test_div_InvalidGeometry(self):
        D = self.D.copy()
        B = D.copy()
        B.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            R = D.div(B)

    def test_mul_InvalidGeometry(self):
        D = self.D.copy()
        B = D.copy()
        B.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            R = D.mul(B)

    def test_add_InvalidGeometry(self):
        D = self.D.copy()
        B = D.copy()
        B.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            R = D.add(B)

    def test_sub_InvalidGeometry(self):
        D = self.D.copy()
        B = D.copy()
        B.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            R = D.sub(B)

    def test_mul_CopyFalse(self):
        D = self.D.copy()
        ref = D.data*D.data
        D.mul(D, copy=False)
        self.assertTrue(np.all(ref-D.data < 1.E-6))

    def test_mul_CopyTrue(self):
        D = self.D.copy()
        ref = D.data*D.data
        R = D.mul(D)
        self.assertTrue(np.all(ref-R.data < 1.E-6))

    def test_divc_Default(self):
        D = self.D.copy()
        R = D.divc(2.)
        d = D.data[:,0,0] *0.5
        self.assertTrue(np.all(d-R.data[:,0,0]) == 0.)


    def test_subc(self):
        D = self.D.copy()
        R = D.subc(10.)
        d = D.data - 10.
        self.assertTrue(np.all(d-R.data) == 0.)

    def test_mulc(self):
        D = self.D.copy()
        R = D.mulc(2.)
        d = D.data[:,0,0] * 2.
        self.assertTrue(np.all(d-R.data[:,0,0]) == 0.)

    def test_ConvertMonthlyTimeSeries_RaisesValueErrorForInvalidCalendar(self):
        data_object= self.D
        data_object.calendar = 'nothing_calendar'
        with self.assertRaises(ValueError):
            data_object._convert_monthly_timeseries()

    def test_apply_temporal_mask_WithInvalidGeometryForMask(self):
        data_object= self.D
        with self.assertRaises(ValueError):
            data_object._apply_temporal_mask(np.random.random((10,20)))

    def test_apply_temporal_mask_WithInvalidMaskTimes(self):
        data_object= self.D
        with self.assertRaises(ValueError):
            data_object._apply_temporal_mask(np.random.random(data_object.nt + 1))


    def test_getunit_ForEmptyUnit(self):
        d = self.D
        d.unit = None
        self.assertEqual(d._get_unit(), '')

    def test_getunit_ForValidUnit(self):
        d = self.D
        d.unit = 'mm/h'
        self.assertEqual(d._get_unit(), '[mm/h]')


    def test_get_percentile_ForInvalidGeometry(self):
        d = self.D
        d.data = np.random.random((10, 20))
        with self.assertRaises(ValueError):
            d.get_percentile(0.5, return_object=True)

    def test_getmindate_ForInvalidBase(self):
        data_object= self.D
        with self.assertRaises(ValueError):
            data_object._get_mindate(base='something')

    def test_getmaxdate_ForInvalidBase(self):
        data_object= self.D
        with self.assertRaises(ValueError):
            data_object._get_maxdate(base='something')

    def test_partial_correlation(self):
        x = self.D
        nt,ny,nx = x.data.shape
        y = x.copy(); y.data = y.data + pl.rand(nt,ny,nx)*1000.
        z = x.copy(); z.data = z.data * pl.rand(nt,ny,nx)*100.

        res  = x.partial_correlation(y, z)
        resarr  = x.partial_correlation(y, z, return_object=False)
        res1 = x.partial_correlation(y, z, ZY=z)  # test with second condition

        #generate reference solution
        slope, intercept, rxy, p_value, std_err = stats.linregress(x.data[:,0,0],y.data[:,0,0])
        slope, intercept, rxz, p_value, std_err = stats.linregress(x.data[:,0,0],z.data[:,0,0])
        slope, intercept, rzy, p_value, std_err = stats.linregress(z.data[:,0,0],y.data[:,0,0])

        ref = (rxy - rxz*rzy) / (np.sqrt(1.-rxz*rxz)*np.sqrt(1.-rzy*rzy))

        self.assertAlmostEqual(ref,res.data[0, 0], places=5)
        self.assertAlmostEqual(ref,resarr[0, 0], places=5)
        self.assertAlmostEqual(ref,res1.data[0, 0], places=5)
        self.assertAlmostEqual(ref,resarr[0, 0], places=5)

    def test_equal_lon(self):
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
        D = self.D.copy()
        # equal longitudes
        x=np.arange(100)
        D.lon = np.zeros((2,100))
        D.lon[0,:] = x
        D.lon[1,:] = x

        r = D._get_unique_lon()
        self.assertTrue(np.all((x-r) == 0.))

    def test_get_unique_lon_Invalid(self):
        D = self.D.copy()
        D.lon = None
        with self.assertRaises(ValueError):
            r = D._get_unique_lon()

    def test_get_unique_lon_InvalidDimension(self):
        D = self.D.copy()
        D.lon = np.random.random((10,20,30))
        with self.assertRaises(ValueError):
            r = D._get_unique_lon()

    def test_get_unique_lon_InvalidLons(self):
        D = self.D.copy()
        D.lon = np.random.random((10,20))
        with self.assertRaises(ValueError):
            r = D._get_unique_lon()



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
                x._temporal_subsetting(0, n)
                y._temporal_subsetting(0, n)

        return x,y


    #-----------------------------------------------------------

    def test_corr_single(self):
        x = self.D.copy()
        y = x.data[:,0,0].copy()*2.
        y += np.random.random(len(y))

        #--- pearson
        slope, intercept, r, prob, sterrest = stats.linregress(y, x.data[:,0,0])
        Rout, Sout, Iout, Pout, Cout = x.corr_single(y)

        self.assertAlmostEqual(r,Rout.data[0,0], 8)
        self.assertAlmostEqual(slope,Sout.data[0,0],8)
        self.assertAlmostEqual(intercept,Iout.data[0,0], 8)
        self.assertAlmostEqual(prob,Pout.data[0,0], 8)

        #--- spearman
        y = x.data[:,0,0].copy()*5.
        y += np.random.random(len(y))*3.

        rho, prob = stats.mstats.spearmanr(y, x.data[:,0,0])
        Rout, Sout, Iout, Pout, Cout = x.corr_single(y, method='spearman')
        self.assertAlmostEqual(r,Rout.data[0,0],5)
        self.assertAlmostEqual(prob,Pout.data[0,0],8)

    def test_corr_single_InvalidGeometry(self):
        x = self.D.copy()
        y = x.data[:,0,0].copy()*2.
        y += np.random.random(len(y))
        x.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            x.corr_single(y)

    def test_corr_single_InvalidMethod(self):
        x = self.D.copy()
        y = x.data[:,0,0].copy()*2.
        y += np.random.random(len(y))
        with self.assertRaises(ValueError):
            x.corr_single(y, method='some_funky_method')

    @unittest.skip('wait for bugfree scipy')
    def test_correlate(self):
        for n in [None,100,10,5]:  # different size
            x,y = self.generate_tuple(n=n,mask=True)
            x1=x.data[:, 0, 0]
            y1=y.data[:, 0, 0]
            msk = (x1.mask == False) & (y1.mask == False)
            x2 = x1[msk]
            y2 = y1[msk]  # this is only the valid data

            #print 'Number of masked pixels: ', sum(y.data.mask), n

            ##################################################################
            # PEARSON CORRELATION
            ##################################################################
            slope, intercept, r_value1, p_value1, std_err = stats.mstats.linregress(x1,y1) #masked
            slope, intercept, r_value2, p_value2, std_err = stats.linregress(x2,y2) #not masked
            r,p = x.correlate(y)

            #1) test if scipy functions return similar results
            self.assertAlmostEqual(r_value1,r_value2,places=10)

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

        #/// linear detrending of data ///
        x = self.D.copy()
        tmp = np.arange(len(x.time))
        tmp = np.ma.array(tmp, mask = tmp != tmp)
        x.data[:,0,0] = np.ma.array(tmp, mask=tmp!=tmp)
        y = x.copy()
        y.data = y.data * 1.2 + 3.
        #~ y.data = np.ma.array(y.data, mask = y.data != y.data)

        r,p = x.correlate(y)
        self.assertAlmostEqual(r.data[0,0], 1., 10)

        #--- detrending ---
        r,p = x.correlate(y, detrend=True)
        self.assertEquals(r.data[0,0], 0.)

    @unittest.skip('wait for bugfree scipy')
    def test_detrend(self):
        x = self.D.copy()
        t = np.arange(len(x.time))
        r = np.random.random(len(x.time))
        y = t * 10.73 + 5.39 + r
        x.data[:,0,0] = np.ma.array(y, mask=y!=y)

        # return object
        xd = x.detrend()
        slope, intercept, r_value, p_value, std_err = stats.linregress(t,y)
        ref = y - (slope*t+intercept)
        d = np.abs(1.-xd.data[:,0,0]/ref)
        self.assertTrue(np.all(d < 1.E-10))

        # no object
        x = self.D.copy()
        x.data[:,0,0] = np.ma.array(y, mask=y!=y)
        x.detrend(return_object=False)
        slope, intercept, r_value, p_value, std_err = stats.linregress(t,y)
        ref = y - (slope*t+intercept)
        d = np.abs(1.-x.data[:,0,0]/ref)
        self.assertTrue(np.all(d < 1.E-10))

    def test_normalize(self):
        x = self.D.copy()
        d = x.data[:,0,0].copy()

        r = (d - d.mean()) / d.std()
        x.normalize(return_object=False)
        dif = np.abs(x.data[:,0,0]-r)
        self.assertTrue(np.all(dif == 0.))

        x = self.D.copy()
        y=x.normalize(return_object=True)
        dif = np.abs(y.data[:, 0, 0]-r)
        self.assertTrue(np.all(dif == 0.))

    def test_normalize_InvalidGeometry(self):
        x = self.D.copy()
        x.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            x.normalize(return_object=False)

    def test_condstat(self):
        """
        conditional statistics unittest
        """

        #sample data
        D = self.D.copy()
        D.data = pl.randn(100,3,1) #some sample data
        msk = np.asarray([[1,1,3],]).T #sample mask

        #calculate conditional statistics
        res = D.condstat(msk)

        #test for mask value == 1 (2 pixels)
        rm = 0.5*(D.data[:,0,0] + D.data[:,1,0])
        rs = (D.data[:,0,0] + D.data[:,1,0])

        self.assertTrue(np.all((res[1]['mean']-rm) == 0. ))
        self.assertTrue(np.all((res[1]['sum']-rs) == 0. ))

        #test for mask value == 3 (1 pixel)
        rm = rs = D.data[:,2,0]
        self.assertTrue(np.all( (res[3]['mean']-rm) == 0. ))
        self.assertTrue(np.all( (res[3]['sum']-rs) == 0. ))


    def test_condstat_InvalidGeometry(self):
        D = self.D.copy()
        D.data = np.random.random((10,20,30,40))
        msk = np.asarray([[1,1,3],]).T
        with self.assertRaises(ValueError):
            res = D.condstat(msk)

    def test_condstat_InvalidGeometryMask(self):
        D = self.D.copy()
        D.data = pl.randn(100,3,1)
        msk = np.asarray([[1,2,4,5],]).T
        with self.assertRaises(ValueError):
            res = D.condstat(msk)

    def test_temporal_subsettingInvalidGeometry(self):
        x = self.D.copy()
        x.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            x._temporal_subsetting(2, 5)

    def test_apply_temporal_subsetting(self):
        # checks only if the right time is subsetted
        import datetime
        x = self.D.copy()

        start_date = datetime.datetime(2003,3,1)
        stop_date = datetime.datetime(2003,5,28)
        x.apply_temporal_subsetting(start_date, stop_date)

        d = x.date
        self.assertEqual(d[0].year,2003)
        self.assertEqual(d[0].month,3)
        self.assertEqual(d[0].day,1)
        self.assertEqual(d[-1].year,2003)
        self.assertEqual(d[-1].month,5)
        self.assertEqual(d[-1].day,28)

    def test_apply_temporal_mask(self):
        D=self.D.copy()
        D.data[:,:,:]=1.
        m = np.zeros(len(D.data)).astype('bool')
        m[1] = True; m[5]=True
        D._apply_temporal_mask(m)

    def test_bounding_box(self):
        D = self.D.copy()
        D.data = np.ma.array(pl.rand(10,5,8),mask=np.zeros((10,5,8)).astype('bool'))

        #generate some sample data with known bounding box
        D.data.mask[:,:,0] = True
        D.data.mask[:,:,7] = True
        D.data.mask[:,0,:] = True
        D.data.mask[:,4,:] = True

        #validate function
        i1,i2,j1,j2 = D.get_bounding_box()
        self.assertEqual(i1,1)
        self.assertEqual(i2,3)
        self.assertEqual(j1,1)
        self.assertEqual(j2,6)

    @unittest.skip('wait for bug free scipy')
    def test_fldmean(self):
        """
        unittest for fldmean() function
        @return:
        """

        # define testdata
        D = self.D
        x = np.ones((1,3,1))
        for i in [0]:
            x [i,0,0] = 5.
            x [i,1,0] = 10.
            x [i,2,0] = 20.
        D.data = np.ma.array(x,mask=x!=x)
        y = np.ones((3,1))
        y[0,0] = 75.
        y[1,0] = 25.
        y[2,0] = 25.
        D.cell_area = y

        D1=D.copy()  # 2D version
        xx = np.ones((3,1))
        xx[0,0]=5.
        xx[1,0]=10.
        xx[2,0]=20.
        D1.data = np.ma.array(xx,mask=xx!=xx)

        # do test
        r1 = D.fldmean()[0]  # with weights
        r1a = D1.fldmean()[0]

        self.assertEqual(r1, 9.)
        self.assertEqual(r1a, 9.)

        r2 = D.fldmean(apply_weights=False)  # without weights
        r2a = D1.fldmean(apply_weights=False)
        self.assertEqual(r2[0],x.mean())
        self.assertEqual(r2a[0],xx.mean())

        # 2D case
        D=self.D.copy()
        x=np.ones((1,4))
        x[0,1] = 1.
        x[0,2] = 5.
        D.data = np.ma.array(x,mask=x==0.)
        ny,nx = x.shape
        ca = np.ones((ny,nx))
        D.cell_area = np.ma.array(ca,mask=ca < 0.)
        r = D.fldmean()[0]
        self.assertEquals(r, 2.)

        # now test against results from CDO
        D.save('tmp_data.nc', delete=True, varname='test')
        cmd = 'cdo -f nc fldmean tmp_data.nc tmp_fldmean.nc'
        os.system(cmd)
        T = Data('tmp_fldmean.nc', 'test', read=True)
        self.assertEquals(r, T.data[0,0])
        self.assertEquals(2., T.data[0,0])
        os.remove('tmp_fldmean.nc')
        os.remove('tmp_data.nc')

        # testcase where some of data is not valid and different weighting approaches are applied
        D = self.D.copy()

        x=np.ones((1,1,4))
        D.data = np.ma.array(x,mask=x==0.)
        nt,ny,nx = x.shape
        ca = np.ones((ny,nx))
        D.cell_area = np.ma.array(ca,mask=ca < 0.)

        D.weighting_type='valid'
        r = D.fldmean()[0]
        self.assertEquals(r,1.)
        x[:,0,0] = np.nan
        D.data = np.ma.array(x,mask=np.isnan(x))
        r = D.fldmean()[0]
        self.assertEquals(r, 1.)

        #... now check what happens if normalization factor is for ALL pixels and not only the valid ones! --> should give 0.75
        D.weighting_type='all'
        r = D.fldmean()[0]
        self.assertEquals(r, 0.75)


    def test_fldstd(self):
        #define testdata
        D = self.D
        x = np.ones((1,3,1))
        for i in [0]:
            x [i,0,0] = 5.; x [i,1,0] = 10.; x [i,2,0] = 20.
        D.data = np.ma.array(x,mask=x!=x)
        y = np.ones((3,1))
        y[0,0] = 75.; y[1,0] = 25.; y[2,0] = 25.
        D.cell_area = y

        # define testcase described under http://en.wikipedia.org/wiki/Weighted_mean#Weighted_sample_variance
        # For example, if values \{2, 2, 4, 5, 5, 5\} are drawn from
        # the same distribution, then we can treat this set as an
        # unweighted sample, or we can treat it as the weighted
        # sample \{2, 4, 5\} with corresponding weights
        # \{2, 1, 3\}, and we should get the same results

        xdat = np.asarray([2., 2., 4., 5., 5., 5.])

        ### 2D data ###

        # 1) no weighting
        A = self.D.copy()
        x = np.ones((1, len(xdat)))
        x[0,:] = xdat*1.
        y = np.ones_like(x)*3.  # cell area dummy
        A.cell_area = y*1.
        A.data = np.ma.array(x, mask=x != x)

        # ddof = 0
        r = A.fldstd(apply_weights=False, ddof=0)
        self.assertEqual(r,xdat.std(ddof=0))
        # ddof = 1
        #~ r = A.fldstd(apply_weights=False, ddof=1)
        #~ self.assertEqual(r,xdat.std(ddof=1))

        # 2) weighting

        # a) same cell size
        r = A.fldstd(apply_weights=True, ddof=0)
        self.assertAlmostEqual(r, xdat.std(ddof=0), 10)

        #~ r = A.fldstd(apply_weights=True, ddof=1)
        #~ self.assertAlmostEqual(r, xdat.std(ddof=1),10)

        # b) different cell size
        refdat = np.asarray([2.,4.,5.])
        x = np.ones((1,3))
        x[0,0] = 2.
        x[0,1] = 4.
        x[0,2] = 5.
        A.data = np.ma.array(x, mask=x != x)
        y = np.ones_like(x)
        y[0,0] = 2. # weight in acordance with the number of
        y[0,1] = 1. # occurences in xdat (se above)
        y[0,2] = 3.
        y = y * 10. # scale cell sizes still a bit
        A.cell_area = y*1.

        # ddof = 0
        r = A.fldstd(apply_weights=True, ddof=0)
        self.assertAlmostEqual(r, xdat.std(ddof=0), 10)

        # ddof = 1
        #~ r = A.fldstd(apply_weights=True, ddof=1)
        #~ self.assertAlmostEqual(r, xdat.std(ddof=1), 10)




        ### 3D data ###
        del A
        A = self.D.copy()
        x = np.ones((3,6,1))
        y = np.ones((6,1))
        x[0,:,0] = xdat*1.
        x[1,:,0] = xdat*1.
        x[2,:,0] = xdat*1.
        A.data = np.ma.array(x, mask= x!=x)
        A.cell_area = y*1.

        #1) no weighting

        # ddof = 0
        r = A.fldstd(apply_weights=False, ddof=0)
        self.assertEqual(r[0],xdat.std(ddof=0))
        self.assertEqual(r[1],xdat.std(ddof=0))
        self.assertEqual(r[2],xdat.std(ddof=0))
        # ddof = 1
        #~ r = A.fldstd(apply_weights=False, ddof=1)
        #~ self.assertEqual(r[0],xdat.std(ddof=1))
        #~ self.assertEqual(r[1],xdat.std(ddof=1))
        #~ self.assertEqual(r[2],xdat.std(ddof=1))

        #2) weighting

        # a) same size
        r = A.fldstd(apply_weights=True, ddof=0)
        #ddof = 0
        self.assertAlmostEqual(r[0], xdat.std(ddof=0),10)
        self.assertAlmostEqual(r[1], xdat.std(ddof=0),10)
        self.assertAlmostEqual(r[2], xdat.std(ddof=0),10)
        #ddof = 1
        #~ r = A.fldstd(apply_weights=True, ddof=1)
        #~ self.assertAlmostEqual(r[0], xdat.std(ddof=1),10) todo does not work, but not sure it std(ddof=1) is the proper reference!
        #~ self.assertAlmostEqual(r[1], xdat.std(ddof=1),10)
        #~ self.assertAlmostEqual(r[2], xdat.std(ddof=1),10)


        # b) different cell size
        B = self.D.copy()
        refdat = np.asarray([2.,4.,5.])
        x = np.ones((3,3,1))
        x[0,:,0] = refdat*1.
        x[1,:,0] = refdat*1.
        x[2,:,0] = refdat*1.

        y = np.ones((3,1))
        y[0,0] = 2.
        y[1,0] = 1.
        y[2,0] = 3.
        B.data = np.ma.array(x, mask= x!=x)
        B.cell_area = y*1.

        # ddof = 0
        r = B.fldstd(apply_weights=True, ddof=0)
        self.assertAlmostEqual(r[0], xdat.std(ddof=0),10)
        self.assertAlmostEqual(r[1], xdat.std(ddof=0),10)
        self.assertAlmostEqual(r[2], xdat.std(ddof=0),10)

        #ddof = 1
        #~ r = B.fldstd(apply_weights=True, ddof=1)
        #~ self.assertAlmostEqual(r[0], xdat.std(ddof=1),10) todo does not work, but not sure it std(ddof=1) is the proper reference!
        #~ self.assertAlmostEqual(r[1], xdat.std(ddof=1),10)
        #~ self.assertAlmostEqual(r[2], xdat.std(ddof=1),10)

        # now test against results from CDO
        #~ D1.save('tmp_data.nc', delete=True, varname='test')
        #~ cmd = 'cdo -f nc fldstd tmp_data.nc tmp_fldstd.nc'
        #~ os.system(cmd)
        #~ T = Data('tmp_fldstd.nc', 'test', read=True)
        #~ print T.data, r1, r1a, ref
        #~ #stop
        #~ self.assertEquals(r1a, T.data[0,0])


    def test_areasum(self):
        """
        unittest for areasum() function
        """

        # define testdata
        D = self.D
        x = np.ones((1,3,1))
        for i in [0]:
            x [i,0,0] = 5.
            x [i,1,0] = 10.
            x [i,2,0] = 20.
        D.data = np.ma.array(x, mask=x!=x)
        y = np.ones((3,1))
        y[0,0] = 75.
        y[1,0] = 25.
        y[2,0] = 25.  # total area = 125.
        D.cell_area = y

        D1=D.copy()  # 2D version
        xx = np.ones((3,1))
        xx[0,0]=5.
        xx[1,0]=10.
        xx[2,0]=20.
        D1.data = np.ma.array(xx, mask=xx!=xx)

        # do test
        r1  = D .areasum()[0] #result should be 5.*75. + 10.*25. + 20.*25.
        r1a = D1.areasum()[0]
        r1d  = D .areasum(return_data=True).data[0]
        r1ad = D1.areasum(return_data=True).data[0]

        self.assertEqual(r1, 5.*75. + 10.*25. + 20.*25.)
        self.assertEqual(r1a, 5.*75. + 10.*25. + 20.*25.)
        self.assertEqual(r1d, 5.*75. + 10.*25. + 20.*25.)
        self.assertEqual(r1ad, 5.*75. + 10.*25. + 20.*25.)

        r2 = D.areasum(apply_weights=False) #without weights
        r2a = D1.areasum(apply_weights=False)
        self.assertEqual(r2[0], x.sum())
        self.assertEqual(r2a[0], xx.sum())


        # 2D case
        D=self.D.copy()
        x=np.ones((1,4))
        x[0,1]=1.; x[0,2] = 5.
        D.data = np.ma.array(x,mask=x==0.)
        ny,nx = x.shape
        ca = np.ones((ny,nx))
        D.cell_area = np.ma.array(ca,mask=ca < 0.)
        r = D.areasum()[0]
        self.assertEquals(r,8.)

        # testcase where some of data is not valid and different weighting approaches are applied
        D = self.D.copy()

        x=np.ones((1,1,4))
        D.data = np.ma.array(x, mask=x==0.)
        nt,ny,nx = x.shape
        ca = np.ones((ny, nx))
        D.cell_area = np.ma.array(ca, mask=ca < 0.)

        D.weighting_type='valid'
        r = D.areasum()[0]
        self.assertEquals(r, 4.)
        x[:,0,0] = np.nan
        D.data = np.ma.array(x,mask=np.isnan(x))
        r = D.areasum()[0]
        self.assertEquals(r, 3.)

        # ... now check what happens if normalization factor is for ALL pixels and not only the valid ones! --> should give 0.75
        D.weighting_type='all'
        r = D.areasum()[0]
        self.assertEquals(r, 3.)

    def test_set_timecycle(self):
        D = self.D

        # set some monthly timeseries
        s_start_time = '2003-01-01'
        s_stop_time  = '2005-12-31'
        start_time = pl.num2date(pl.datestr2num(s_start_time))
        stop_time  = pl.num2date(pl.datestr2num(s_stop_time ))
        tref = rrule(MONTHLY, dtstart = start_time).between(start_time, stop_time, inc=True) #monthly timeseries
        D.time = pl.date2num(tref)

        #1) a perfect monthly timeseries
        #check that that timeseries is based on monthly data
        self.assertTrue(D._is_monthly())
        D._set_timecycle()
        self.assertEquals(D.time_cycle,12)

        #2) some timeseries that is not monthly
        D.time_cycle=None
        D.time[2]=pl.datestr2num('2010-05-01')
        D._set_timecycle()
        self.assertFalse(D._is_monthly())
        self.assertEquals(D.time_cycle,None)

        #3) some timeseries that is has increasing months, but wrong years!
        D.time_cycle=None
        D.time = pl.date2num(tref)
        t = pl.num2date(D.time[2])
        D.time[2]=pl.datestr2num('2010-' + str(t.month).zfill(2) + '-' + str(t.day).zfill(2))
        D._set_timecycle()
        self.assertFalse(D._is_monthly())
        self.assertEquals(D.time_cycle,None)

    def test_get_valid_mask(self):
        D = self.D.copy()

        #case 1: 2D data
        x = np.ones((1,2))
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask()
        self.assertTrue(m[0,0]==True)
        self.assertTrue(m[0,1]==True)

        #case 2: 3D with all valid data
        D = self.D.copy()
        x = np.ones((50,1,2))
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask()
        self.assertTrue(m[0,0]==True)
        self.assertTrue(m[0,1]==True)

        #case 3: some invalid data at one pixel (frac=1=default)
        D = self.D.copy()
        x = np.ones((50,1,2))
        x[0:25,0,0] = 0.
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask()
        self.assertTrue(m[0,0]==False)
        self.assertTrue(m[0,1]==True)

        #case 4 exactly 50% invalid
        D = self.D.copy()
        x = np.ones((50,1,2))
        x[0:25,0,0] = 0.
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask(frac=0.5)
        self.assertTrue(m[0,0]==True)
        self.assertTrue(m[0,1]==True)

        #case 5: <50% valid
        D = self.D.copy()
        x = np.ones((50,1,2))
        x[0:26,0,0] = 0.
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask(frac=0.5)
        self.assertTrue(m[0,0]==False)
        self.assertTrue(m[0,1]==True)

        #case 6: 1D data (all valid)
        x = np.ones(100)
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask()
        self.assertTrue(m[0,0]==True)

        #case 7: 1D data (51% invalid)
        x = np.ones(100)
        x[0:51] = 0.
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask(frac=0.5)
        self.assertTrue(m[0,0]==False)

        #case 7: 1D data (50% invalid)
        x = np.ones(100)
        x[0:50] = 0.
        D.data = np.ma.array(x,mask=x == 0)
        m = D.get_valid_mask(frac=0.5)
        self.assertTrue(m[0,0]==True)

    def test_time_conversion(self):
        x = self.D.copy()
        t = x.time
        dref = pl.num2date(t)
        t2=pl.date2num(dref)
        d1 = x.num2date(t)  # convert time to datetime object
        t1 = x.date2num(d1)  # convert back

        d = t-t1
        self.assertTrue(np.all(d == 0.))
        d = t-t2
        self.assertTrue(np.all(d == 0.))

    def test_align(self):
        """test temporal alignment of two datasets"""
        x = self.D.copy()
        y = self.D.copy()
        y.subc(2.5, copy=False)
        y._temporal_subsetting(500, 750)  # generate a subset dataset

        x1, y1 = x.align(y, base='day')  # do temporal alignment
        d = x1.sub(y1).divc(2.5).subc(1.)  # should give small diff.

        # check dates
        self.assertEqual(x1.date[0], y1.date[0])
        self.assertEqual(x1.date[-1], y1.date[-1])

        # check that really the same data is used
        self.assertTrue(np.all(np.abs(d.data) < 0.00000001))

        #... and the other way round
        y1, x1 = y.align(x, base='day')
        d = x1.sub(y1).divc(2.5).subc(1.)
        self.assertEqual(x1.date[0], y1.date[0])
        self.assertEqual(x1.date[-1], y1.date[-1])
        self.assertTrue(np.all(np.abs(d.data) < 0.00000001))

        # test monthly with invalid data
        with self.assertRaises(ValueError):
            y1, x1 = y.align(x, base='month')

    def test_align_unsorted_X(self):
        x = self.D.copy()
        y = self.D.copy()
        y.subc(2.5, copy=False)
        y._temporal_subsetting(500, 750)  # generate a subset dataset
        x.time = np.random.random(x.nt)  # generate unsorted time
        with self.assertRaises(ValueError):
            x1, y1 = x.align(y, base='day')  # do temporal alignment

    def test_align_unsorted_X(self):
        x = self.D.copy()
        y = self.D.copy()
        y.subc(2.5, copy=False)
        y._temporal_subsetting(500, 750)  # generate a subset dataset
        y.time = np.random.random(y.nt)  # generate unsorted time
        with self.assertRaises(ValueError):
            x1, y1 = x.align(y, base='day')  # do temporal alignment

    def test_align_InvalidBase(self):
        """test temporal alignment of two datasets"""
        x = self.D.copy()
        y = self.D.copy()
        y.subc(2.5, copy=False)
        y._temporal_subsetting(500, 750)  # generate a subset dataset
        with self.assertRaises(ValueError):
            x1, y1 = x.align(y, base='some_invalid_base')
        with self.assertRaises(ValueError):
            x1, y1 = x.align(y, base=None)




    def test_is_daily(self):
        x = self.D.copy()  # is already daily
        self.assertTrue(x._is_daily())

        x = self.D.copy()
        x.time[2] += 4.
        self.assertFalse(x._is_daily())

    def test_is_daily_WithOutTime(self):
        x = self.D.copy()  # is already daily
        del x.time
        self.assertFalse(x._is_daily())

    def test_is_sorted(self):
        x = self.D.copy()
        self.assertTrue(x._is_sorted())

        x.time[5] += 10.
        self.assertFalse(x._is_sorted())

    def test_days_per_month(self):
        x = self.D.copy()
        x.time[0] = pl.datestr2num('1999-03-15')
        x.time[1] = pl.datestr2num('2000-03-15')
        x.time[2] = pl.datestr2num('2002-09-15')
        x.time[3] = pl.datestr2num('2000-02-15')
        x.time[4] = pl.datestr2num('2003-02-15')

        d = x._days_per_month()
        self.assertEqual(d[0], 31)
        self.assertEqual(d[1], 31)
        self.assertEqual(d[2], 30)
        self.assertEqual(d[3], 29)
        self.assertEqual(d[4], 28)

    def test_temporal_smooth(self):
        """
        test smooth routine
        """
        x = self.D.copy()

        #--- TEST for 1D data ---
        tmp = np.random.random(1000)
        x.data = np.ma.array(tmp, mask=tmp!=tmp)


        # windowsize 2
        with self.assertRaises(ValueError):
            xxx = x.temporal_smooth(2)

        # windowsize 3
        y3a = x.temporal_smooth(3)
        y3b = x.temporal_smooth(3, return_object=False)
        self.assertEqual(y3a.data[10], y3b[10])
        self.assertAlmostEqual(tmp[10:13].sum()/3., y3a.data[11], 8)

        # windowsize 5
        y5a = x.temporal_smooth(5)
        y5b = x.temporal_smooth(5, return_object=False)
        self.assertEqual(y5a.data[20], y5b[20])
        self.assertAlmostEqual(tmp[30:35].sum()/5., y5a.data[32], 8)

        #--- TEST FOR 3D data ---
        tmp = np.random.random((100, 2, 3))
        x.data = np.ma.array(tmp, mask=tmp!=tmp)
        y3a = x.temporal_smooth(3)
        y3b = x.temporal_smooth(3, return_object=False)
        self.assertEqual(y3a.data[10,0,0], y3b[10,0,0])
        self.assertAlmostEqual(tmp[10:13,0,0].sum()/3., y3a.data[11,0,0], 8)
        self.assertAlmostEqual(tmp[10:13,1,1].sum()/3., y3a.data[11,1,1], 8)
        self.assertAlmostEqual(tmp[10:13,1,0].sum()/3., y3a.data[11,1,0], 8)


    def test_hp_filter_InvalidLambda(self):
        with self.assertRaises(ValueError):
            self.D.hp_filter(-10., return_object=True)

    def test_hp_filter_Invalid2D(self):
        x = self.D.copy()
        x.data = np.random.random((20, 30))
        with self.assertRaises(ValueError):
            x.hp_filter(100, return_object=True)


    def test_areasum_InvalidGeometry(self):
        x = self.D.copy()
        x.data = np.random.random((10,20,30,40))
        with self.assertRaises(ValueError):
            x.areasum()

if __name__ == '__main__':
    unittest.main()
