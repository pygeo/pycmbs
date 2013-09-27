from unittest import TestCase

__author__ = 'm300028'


from data import *
from diagnostic import *
from plots import *
import scipy as sc
import pylab as pl
import numpy as np
from region import *
####from pyCMBS import *

class TestData(TestCase):

    def setUp(self):
        #init Data object for testing
        n=1000 #slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100. #generate dummy data
        self.D = Data(None,None)
        d=np.ones((n,1,1))
        self.D.data = d; self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data,mask=self.D.data != self.D.data)
        self.D.verbose = True; self.D.unit = 'myunit'
        self.D.label = 'testlabel'; self.D.filename = 'testinputfilename.nc'
        self.D.varname = 'testvarname'; self.D.long_name = 'This is the longname'
        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01') - 1


    def test_gleckler_index(self):
        """
        test Reichler index/Gleckler plot
        @return:
        """

        #--- generate sample data ---
        tmp = np.zeros((5,3,1)); tmp[:,0,0] = np.ones(5)*1.
        tmp[:,1,0] = np.ones(5)*2.; tmp[:,2,0] = np.ones(5)*5.

        x = self.D.copy()
        x.data = np.ma.array(tmp,mask=tmp!=tmp)
        x.std  = np.ones(x.data.shape)

        tmp = np.ones((3,1))
        tmp[1,0] = 2.; x.cell_area = tmp*1.

        y = self.D.copy()
        tmp = np.ones(x.data.shape)
        y.data = np.ma.array(tmp,mask=tmp!=tmp)

        #--- diagnostic ---
        D = GlecklerPlot()
        r=D.calc_index(x,y,'a','b')
        self.assertEqual(r,np.sqrt(22.5))

        #... use different std
        x.std  = np.ones(x.data.shape)
        x.std[:,2,0] = 0.5
        r=D.calc_index(x,y,'a','b')
        self.assertEqual(r,np.sqrt(42.5))


