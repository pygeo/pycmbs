"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

from unittest import TestCase
from pycmbs.data import Data
#~ from pycmbs.diagnostic import *
import scipy as sc
import pylab as pl
import numpy as np

__author__ = 'm300028'

class TestEOF(TestCase):

    def setUp(self):
        #init Data object for testing
        n=4 #slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100. #generate dummy data
        self.D = Data(None,None)
        d=np.ones((n,1,2))
        self.D.data = d
        self.D.data[:,0,0]=x
        self.D.data = np.ma.array(self.D.data,mask=self.D.data != self.D.data)
        self.D.verbose = True
        self.D.unit = 'myunit'

        self.D.time = np.arange(n) + pl.datestr2num('2001-01-01') - 1


    def test_eof_analysis(self):
        #test of EOF
        #example taken from:
        #   http://www.atmos.washington.edu/~dennis/552_Notes_4.pdf , page 80

        #assign sampled data
        D = self.D.copy()
#        d1 = np.array([2,4,-6,8])
#        d2 = np.array([1,2,-3,4])
        #d1 = np.array([2,1])
        #d2 = np.array([4,2])
        #d3 = np.array([-6,-3])
        #d4 = np.array([8,4])

        #D.data[:,0,0] = d1; D.data[:,0,1] = d2
        #D.data[:,0,2] = d3; D.data[:,0,3] = d4

#        D.data[:,0,0] = d1
#        D.data[:,0,1] = d2
#
#        E = EOF(D,cov_norm=False) #do not normalize covariance matrix, as this is also not done in example
#        print E.C #covariance matrix


#        shape of EOF is wrong !!!!!!!
#
#        something is not really working here !!!!
#
#
#        irgendwie ist hier ein problem
#        warum einmal 4x4 und einmal 2x2 ???

#        print 'eigval'
#        print E.eigval
#
#        print 'eigvec'
#        print E.eigvec






    pass
