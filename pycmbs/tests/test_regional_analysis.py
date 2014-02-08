from unittest import TestCase
import unittest

from pycmbs.data import *
from pycmbs.diagnostic import RegionalAnalysis
import scipy as sc
import numpy as np


class TestData(TestCase):

    def setUp(self):
        # init Data object for testing
        n = 1000  # slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100.  # generate dummy data
        self.D = Data(None, None)
        d = np.ones((n, 1, 1))
        self.D.data = d
        self.D.data[:, 0, 0] = x
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
        self.D.cell_area = np.ones((1, 1))

    def test_regional_analysis(self):

        # generate two datasets
        ny = 2
        nx = 6
        nt = 500
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

        tmp = np.random.random((nt, ny, nx))
        x.data = np.ma.array(tmp, mask=tmp != tmp)
        x.cell_area = cell_area.copy()
        x.time = x.time[0:nt]
        del tmp

        y = self.D.copy()
        tmp = np.random.random((nt, ny, nx))
        y.data = np.ma.array(tmp, mask=tmp != tmp)
        y.cell_area = cell_area.copy()
        y.time = y.time[0:nt]
        del tmp

        # todo unittest extension for different area weighting !!!!

        reg = Data(None,None)
        reg.data = m
        REGSTAT = RegionalAnalysis(x, y, reg)
        REGSTAT.calculate()

        #//////////////////////////////////////////
        # Simple first and second order statistics
        #//////////////////////////////////////////

        # mask = 1
        refx = (x.data[:, 0, 0] + x.data[:, 1, 0]) * 0.5
        refy = (y.data[:, 0, 0] + y.data[:, 1, 0]) * 0.5
        self.assertTrue(np.all((refx - REGSTAT.statistics['xstat'][1]['mean']) == 0.))
        self.assertTrue(np.all((refy - REGSTAT.statistics['ystat'][1]['mean']) == 0.))
        del refx, refy

        # mask = 4
        # mean
        refx = (x.data[:, 1, 3] + x.data[:, 1, 5] + x.data[:, 0, 4]) / 3.
        refy = (y.data[:, 1, 3] + y.data[:, 1, 5] + y.data[:, 0, 4]) / 3.
        self.assertTrue(np.all(np.abs(1. - refx / REGSTAT.statistics['xstat'][4]['mean']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - refy / REGSTAT.statistics['ystat'][4]['mean']) < 0.000001))

        # std
        stdx = np.asarray([x.data[:, 1, 3], x.data[:, 1, 5], x.data[:, 0, 4]]).std(axis=0)
        stdy = np.asarray([y.data[:, 1, 3], y.data[:, 1, 5], y.data[:, 0, 4]]).std(axis=0)
        self.assertTrue(np.all(np.abs(1. - stdx / REGSTAT.statistics['xstat'][4]['std']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - stdy / REGSTAT.statistics['ystat'][4]['std']) < 0.000001))

        maxx = np.asarray([x.data[:, 1, 3], x.data[:, 1, 5], x.data[:, 0, 4]]).max(axis=0)
        maxy = np.asarray([y.data[:, 1, 3], y.data[:, 1, 5], y.data[:, 0, 4]]).max(axis=0)
        self.assertTrue(np.all(np.abs(1. - maxx / REGSTAT.statistics['xstat'][4]['max']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - maxy / REGSTAT.statistics['ystat'][4]['max']) < 0.000001))

        minx = np.asarray([x.data[:, 1, 3], x.data[:, 1, 5], x.data[:, 0, 4]]).min(axis=0)
        miny = np.asarray([y.data[:, 1, 3], y.data[:, 1, 5], y.data[:, 0, 4]]).min(axis=0)
        self.assertTrue(np.all(np.abs(1. - minx / REGSTAT.statistics['xstat'][4]['min']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - miny / REGSTAT.statistics['ystat'][4]['min']) < 0.000001))

        #//////////////////////////////////////////
        # Correlation statistic
        # Three approaches
        #//////////////////////////////////////////

        # A) calculate once correlation and then calculate regional statistics
        #    as average of the statistical scores
        corrstat = REGSTAT.statistics['corrstat']
        r, p = x.correlate(y)

        # mask = 2
        self.assertAlmostEqual(r.data[:, 1:3].mean(), corrstat['analysis_A'][2]['mean'], 8)
        self.assertAlmostEqual(r.data[:, 1:3].std(), corrstat['analysis_A'][2]['std'], 8)
        self.assertAlmostEqual(r.data[:, 1:3].min(), corrstat['analysis_A'][2]['min'], 8)
        self.assertAlmostEqual(r.data[:, 1:3].max(), corrstat['analysis_A'][2]['max'], 8)
        self.assertAlmostEqual(r.data[:, 1:3].sum(), corrstat['analysis_A'][2]['sum'], 8)

        # B) calculate regional statistics based on entire dataset for a region
        #    This means that all data from all timesteps and locations is used
        #    an a single vector is built, which is then used for comparison

        # mask = 3
        xvec = []
        xvec.append(x.data[:, 0, 3])
        xvec.append(x.data[:, 0, 5])
        xvec.append(x.data[:, 1, 4])
        yvec = []
        yvec.append(y.data[:, 0, 3])
        yvec.append(y.data[:, 0, 5])
        yvec.append(y.data[:, 1, 4])

        xvec = np.asarray(xvec).flatten()
        yvec = np.asarray(yvec).flatten()
        slope, intercept, r_value, p_value, std_err = sc.stats.linregress(xvec, yvec)

        self.assertLess(abs(1. - slope / corrstat['analysis_B'][3]['slope']), 0.000000000001)
        self.assertLess(abs(1. - intercept / corrstat['analysis_B'][3]['intercept']), 0.000000000001)
        self.assertLess(abs(1. - r_value / corrstat['analysis_B'][3]['correlation']), 0.000000000001)

        # todo: note that it is currently not usefull to compare pvalues, due to the insufficient
        # implementation of the p-value in mstats
        # see this issue: https://github.com/scipy/scipy/pull/3084
        # self.assertLess(abs(1. - p_value / corrstat['analysis_B'][3]['pvalue']), 0.000000000001)

        # C) fldmean() for each region and then correlate
        #    Calculate first the mean time series vecotr for each region and then
        #    do correlation analysis

        # mask = 4
        xvec = (x.data[:, 0, 4] + x.data[:, 1, 3] + x.data[:, 1, 5]) / 3.
        yvec = (y.data[:, 0, 4] + y.data[:, 1, 3] + y.data[:, 1, 5]) / 3.

        slope, intercept, r_value, p_value, std_err = sc.stats.linregress(xvec, yvec)
        self.assertLess(abs(1. - slope / corrstat['analysis_C'][4]['slope']), 0.000000000001)
        self.assertLess(abs(1. - intercept / corrstat['analysis_C'][4]['intercept']), 0.000000000001)
        self.assertLess(abs(1. - r_value / corrstat['analysis_C'][4]['correlation']), 0.000000000001)

        # todo: note that it is currently not usefull to compare pvalues, due to the insufficient
        # see above !

        ##############################
        # SAVE
        REGSTAT.save('testprefix', format='pkl', dir='./xxxtest1')  # save as PKL
        REGSTAT.save('testprefix', format='txt', dir='./xxxtest2')  # save as ASCII

        # ... now check if saved data is o.k
        #1) standard statistics
        fname = './xxxtest2/testprefix_regional_statistics_standard_' + str(3).zfill(16) + '.txt'
        d = np.loadtxt(fname, skiprows=1)

        self.assertTrue(np.all(np.abs(1. - d[:,1] / REGSTAT.statistics['xstat'][3]['mean']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,2] / REGSTAT.statistics['ystat'][3]['mean']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,3] / REGSTAT.statistics['xstat'][3]['std']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,4] / REGSTAT.statistics['ystat'][3]['std']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,5] / REGSTAT.statistics['xstat'][3]['min']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,6] / REGSTAT.statistics['ystat'][3]['min']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,7] / REGSTAT.statistics['xstat'][3]['max']) < 0.000001))
        self.assertTrue(np.all(np.abs(1. - d[:,8] / REGSTAT.statistics['ystat'][3]['max']) < 0.000001))
        del d

        #2) correlation statistics: A
        fname = './xxxtest2/testprefix_regional_statistics_correlation_A.txt'
        d = np.loadtxt(fname, skiprows=1)  # | id | rmean | rstd | rsum | rmin | rmax |
        ids = d[:, 0]
        m = ids == 2
        rmean = d[:, 1][m][0]
        rstd = d[:, 2][m][0]
        rsum = d[:, 3][m][0]
        rmin = d[:, 4][m][0]
        rmax = d[:, 5][m][0]

        self.assertLess(np.abs(1. - rmean / REGSTAT.statistics['corrstat']['analysis_A'][2]['mean'][0]), 0.0000000001)
        self.assertLess(np.abs(1. - rstd / REGSTAT.statistics['corrstat']['analysis_A'][2]['std'][0]), 0.0000000001)
        self.assertLess(np.abs(1. - rsum / REGSTAT.statistics['corrstat']['analysis_A'][2]['sum'][0]), 0.0000000001)
        self.assertLess(np.abs(1. - rmin / REGSTAT.statistics['corrstat']['analysis_A'][2]['min'][0]), 0.0000000001)
        self.assertLess(np.abs(1. - rmax / REGSTAT.statistics['corrstat']['analysis_A'][2]['max'][0]), 0.0000000001)

        # correlation statistics: B
        fname = './xxxtest2/testprefix_regional_statistics_correlation_B.txt'
        d = np.loadtxt(fname, skiprows=1)  # | id | slope | intercept | correlation | pvalue |
        ids = d[:, 0]
        m = ids == 4
        slope = d[:, 1][m][0]
        intercept = d[:, 2][m][0]
        correlation = d[:, 3][m][0]
        # pvalue = d[:, 4][m][0] #todo

        self.assertLess(np.abs(1. - slope / REGSTAT.statistics['corrstat']['analysis_B'][4]['slope']), 0.0000000001)
        self.assertLess(np.abs(1. - intercept / REGSTAT.statistics['corrstat']['analysis_B'][4]['intercept']), 0.0000000001)
        self.assertLess(np.abs(1. - correlation / REGSTAT.statistics['corrstat']['analysis_B'][4]['correlation']), 0.0000000001)
        del d

        # correlation statistics: C
        fname = './xxxtest2/testprefix_regional_statistics_correlation_C.txt'
        d = np.loadtxt(fname, skiprows=1)  # | id | slope | intercept | correlation | pvalue |
        ids = d[:, 0]
        m = ids == 3
        slope = d[:, 1][m][0]
        intercept = d[:, 2][m][0]
        correlation = d[:, 3][m][0]
        # pvalue = d[:, 4][m][0] #todo

        self.assertLess(np.abs(1. - slope / REGSTAT.statistics['corrstat']['analysis_C'][3]['slope']), 0.0000000001)
        self.assertLess(np.abs(1. - intercept / REGSTAT.statistics['corrstat']['analysis_C'][3]['intercept']), 0.0000000001)
        self.assertLess(np.abs(1. - correlation / REGSTAT.statistics['corrstat']['analysis_C'][3]['correlation']), 0.0000000001)

        # Clean up
        os.system('rm -r ./xxxtest1')
        os.system('rm -r ./xxxtest2')
