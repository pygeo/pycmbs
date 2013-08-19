from unittest import TestCase
import unittest

__author__ = 'Walter Sauf'

#identify pyCMBS path and add it to pythonpath, as otherwise the modules are not found properly!

from data import *
from data4D import *
#from diagnostic import *
import scipy as sc
#import pylab as pl
import numpy as np
#from scipy import stats
#import Nio
#from dateutil.rrule import *

class TestData4D(TestCase):

    def setUp(self):
        #init Data object for testing
        n=1000 #slows down significantly! constraint is percentile  test
        x = sc.randn(n)*100. #generate dummy data
        d=np.ones((n,1,1))
        
        D0 = Data(None,None)        
        D0.data = d
        D0.data[:,0,0]=x
        D0.data = np.ma.array(D0.data,mask=D0.data != D0.data)
        
        self.D4d = Data4D(None,None)
        self.D4d.data4D.append(D0.data)
        self.D4d.levellist[1] = 0
        
        x = sc.randn(n)*100. #generate dummy data
        d=np.ones((n,1,1))
        D1 = Data(None,None)        
        D1.data = d
        D1.data[:,0,0]=x
        D1.data = np.ma.array(D1.data,mask=D1.data != D1.data)
        self.D4d.data4D.append(D1.data)
        self.D4d.levellist[2] = 1
        
        self.D4d.verbose = True
        self.D4d.unit = 'myunit'
        self.D4d.label = 'testlabel'
        self.D4d.filename = 'testinputfilename.nc'
        self.D4d.varname = 'testvarname'
        self.D4d.long_name = 'This is the longname'

        self.D4d.time = np.arange(n) + pl.datestr2num('2001-01-01') - 1
        self.D4d.calendar = 'standard'
        
#  Define a normal Data type array
        self.Da = Data(None,None)
        d=np.ones((n,1,1))
        self.Da.data = d
        self.Da.data[:,0,0]=x
        self.Da.data = np.ma.array(self.Da.data,mask=self.Da.data != self.Da.data)
        self.Da.verbose = True
        self.Da.unit = 'my_data_unit'
        self.Da.label = 'testlabel data'
        self.Da.filename = 'testinputfilename.nc'
        self.Da.varname = 'testvarname_data'
        self.Da.long_name = 'This is the longname for data'

        self.Da.time = np.arange(n) + pl.datestr2num('2001-01-01') - 1
        self.Da.calendar = 'standard'
        
# At the moment not included

    def test_get_time_indices(self):
        print "\nDoing test get_time_indices"
        d1 = pl.num2date(pl.datestr2num('2001-01-05'))
        d2 = pl.num2date(pl.datestr2num('2001-05-05'))
        i1,i2 = self.D4d._get_time_indices(d1,d2)
        s1 = str(pl.num2date(self.D4d.time[i1]))
        s2 = str(pl.num2date(self.D4d.time[i2]))

        print s1, i1
        print s2, i2

        self.assertEqual(s1,'2001-01-05 00:00:00+00:00')
        self.assertEqual(s2,'2001-05-05 00:00:00+00:00')


#================================================================================

    def test_add(self):
        print "\nDoing Test add"
        
        print "  - Check add with copy=True on two levels and both type are data4D"
        r1 = self.D4d.add(self.D4d,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0],self.D4d.data4D[0][4,0,0]*2.,places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0],self.D4d.data4D[1][4,0,0]*2.,places=7)
        
        print "  - Check add with copy=False on two levels and both type are data4D"
        ref1 = self.D4d.data4D[0][5,0,0]
        ref2 = self.D4d.data4D[1][5,0,0]
        self.D4d.add(self.D4d,copy=False)
        self.assertAlmostEqual(ref1*2,self.D4d.data4D[0][5,0,0],places=7)
        self.assertAlmostEqual(ref2*2,self.D4d.data4D[1][5,0,0],places=7)

        print "  - Check add with copy=True on two levels with a normal data type to da data4D type"
        r1 = self.D4d.add(self.Da,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0],self.D4d.data4D[0][4,0,0]+self.Da.data[4,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0],self.D4d.data4D[1][4,0,0]+self.Da.data[4,0,0],places=7)

        print "  - Check add with copy=False on two levels with normal data type to da data4D type"
        ref1 = self.D4d.data4D[0][7,0,0]
        ref2 = self.D4d.data4D[1][7,0,0]
        self.D4d.add(self.Da,copy=False)
        self.assertAlmostEqual(ref1+self.Da.data[7,0,0],self.D4d.data4D[0][7,0,0],places=7)
        self.assertAlmostEqual(ref2+self.Da.data[7,0,0],self.D4d.data4D[1][7,0,0],places=7)

#================================================================================

    def test_addc(self):
        print "\nDoing Test addc"
        
        print "  - Check addc with copy=True on two levels and both type are data4D"
        r1 = self.D4d.addc(5.,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0]-5.,self.D4d.data4D[0][4,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0]-5.,self.D4d.data4D[1][4,0,0],places=7)
        
        print "  - Check addc with copy=False on two levels and both type are data4D"
        ref1 = self.D4d.data4D[0][5,0,0]
        ref2 = self.D4d.data4D[1][5,0,0]
        self.D4d.addc(666.,copy=False)
        self.assertAlmostEqual(ref1+666.,self.D4d.data4D[0][5,0,0])
        self.assertAlmostEqual(ref2+666.,self.D4d.data4D[1][5,0,0])
        
#================================================================================

    def test_sub(self):
        print "\nDoing Test sub"
        
        print "  - Check sub with copy=True on two levels and both type are data4D"
        r1 = self.D4d.sub(self.D4d,copy=True)
        r2 = r1.sub(self.D4d,copy=True)
        self.assertAlmostEqual(r2.data4D[0][4,0,0],self.D4d.data4D[0][4,0,0]*-1.,places=7)
        self.assertAlmostEqual(r2.data4D[1][4,0,0],self.D4d.data4D[1][4,0,0]*-1.,places=7)
        
        print "  - Check sub with copy=False on two levels and both type are data4D"
        ref1 = self.D4d.data4D[0][5,0,0]
        ref2 = self.D4d.data4D[1][5,0,0]
        r1 = self.D4d.copy()
        self.D4d.sub(r1,copy=False)
        self.D4d.sub(r1,copy=False)
        self.assertAlmostEqual(ref1*-1,self.D4d.data4D[0][5,0,0],places=7)
        self.assertAlmostEqual(ref2*-1,self.D4d.data4D[1][5,0,0],places=7)
        
        print "  - Check sub with copy=True on two levels with a normal data type to da data4D type"
        r1 = self.D4d.sub(self.Da,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0],self.D4d.data4D[0][4,0,0]-self.Da.data[4,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0],self.D4d.data4D[1][4,0,0]-self.Da.data[4,0,0],places=7)
 
        print "  - Check sub with copy=False on two levels with normal data type to da data4D type"
        ref1 = self.D4d.data4D[0][7,0,0]
        ref2 = self.D4d.data4D[1][7,0,0]
        self.D4d.sub(self.Da,copy=False)
        self.assertAlmostEqual(ref1-self.Da.data[7,0,0],self.D4d.data4D[0][7,0,0],places=7)
        self.assertAlmostEqual(ref2-self.Da.data[7,0,0],self.D4d.data4D[1][7,0,0],places=7)
        
#================================================================================

    def test_subc(self):
        print "\nDoing Test subc"
        
        print "  - Check subc with copy=True on two levels"
        c = 8.
        r1 = self.D4d.subc(c,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0]+c,self.D4d.data4D[0][4,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0]+c,self.D4d.data4D[1][4,0,0],places=7)
        
        print "  - Check subc with copy=False on two levels"
        c = 45.
        ref1 = self.D4d.data4D[0][5,0,0]
        ref2 = self.D4d.data4D[1][5,0,0]
        self.D4d.subc(c,copy=False)
        self.assertAlmostEqual(ref1-c,self.D4d.data4D[0][5,0,0],places=7)
        self.assertAlmostEqual(ref2-c,self.D4d.data4D[1][5,0,0],places=7)
        

#================================================================================

    def test_mul(self):
        print "\nDoing Test mul"
        
        print "  - Check mul with copy=True on two levels and both type are data4D"
        r1 = self.D4d.mul(self.D4d,copy=True)
        self.assertAlmostEqual(r1.data4D[0][9,0,0],self.D4d.data4D[0][9,0,0]*self.D4d.data4D[0][9,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][9,0,0],self.D4d.data4D[1][9,0,0]*self.D4d.data4D[1][9,0,0],places=7)
        
        print "  - Check mul with copy=False on two levels and both type are data4D"
        ref1 = self.D4d.data4D[0][11,0,0]
        ref2 = self.D4d.data4D[1][11,0,0]
        self.D4d.mul(self.D4d,copy=False)
        self.assertAlmostEqual(ref1*ref1,self.D4d.data4D[0][11,0,0],places=7)
        self.assertAlmostEqual(ref2*ref2,self.D4d.data4D[1][11,0,0],places=7)

        print "  - Check mul with copy=True on two levels with a normal data type to da data4D type"
        r1 = self.D4d.mul(self.Da,copy=True)
        self.assertAlmostEqual(r1.data4D[0][1,0,0],self.D4d.data4D[0][1,0,0]*self.Da.data[1,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][1,0,0],self.D4d.data4D[1][1,0,0]*self.Da.data[1,0,0],places=7)

        print "  - Check mul with copy=False on two levels with normal data type to da data4D type"
        ref1 = self.D4d.data4D[0][7,0,0]
        ref2 = self.D4d.data4D[1][7,0,0]
        self.D4d.mul(self.Da,copy=False)
        self.assertAlmostEqual(ref1*self.Da.data[7,0,0],self.D4d.data4D[0][7,0,0],places=7)
        self.assertAlmostEqual(ref2*self.Da.data[7,0,0],self.D4d.data4D[1][7,0,0],places=7)
        
#================================================================================

    def test_mulc(self):
        print "\nDoing Test mulc"
        
        print "  - Check addc with copy=True on two levels and both type are data4D"
        r1 = self.D4d.addc(5.,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0]-5.,self.D4d.data4D[0][4,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0]-5.,self.D4d.data4D[1][4,0,0],places=7)
        
        print "  - Check addc with copy=False on two levels and both type are data4D"
        ref1 = self.D4d.data4D[0][5,0,0]
        ref2 = self.D4d.data4D[1][5,0,0]
        self.D4d.addc(666.,copy=False)
        self.assertAlmostEqual(ref1+666.,self.D4d.data4D[0][5,0,0],places=7)
        self.assertAlmostEqual(ref2+666.,self.D4d.data4D[1][5,0,0],places=7)
        
#================================================================================

    def test_div(self):
        print "\nDoing Test div"
        
        print "  - Check div with copy=True on two levels and both type are data4D"
        r1 = self.D4d.div(self.D4d,copy=True)
        self.assertAlmostEqual(r1.data4D[0][12,0,0],self.D4d.data4D[0][12,0,0]/self.D4d.data4D[0][12,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][12,0,0],self.D4d.data4D[1][12,0,0]/self.D4d.data4D[1][12,0,0],places=7)
        
        print "  - Check div with copy=False on two levels and both type are data4D"
        ref1 = self.D4d.data4D[0][11,0,0]
        ref2 = self.D4d.data4D[1][11,0,0]
        self.D4d.div(self.D4d,copy=False)
        self.assertAlmostEqual(ref1/ref1,self.D4d.data4D[0][11,0,0],places=7)
        self.assertAlmostEqual(ref2/ref2,self.D4d.data4D[1][11,0,0],places=7)

        print "  - Check div with copy=True on two levels with a normal data type to da data4D type"
        r1 = self.D4d.div(self.Da,copy=True)
        self.assertAlmostEqual(r1.data4D[0][1,0,0],self.D4d.data4D[0][1,0,0]/self.Da.data[1,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][1,0,0],self.D4d.data4D[1][1,0,0]/self.Da.data[1,0,0],places=7)

        print "  - Check div with copy=False on two levels with normal data type to da data4D type"
        ref1 = self.D4d.data4D[0][7,0,0]
        ref2 = self.D4d.data4D[1][7,0,0]
        self.D4d.div(self.Da,copy=False)
        self.assertAlmostEqual(ref1/self.Da.data[7,0,0],self.D4d.data4D[0][7,0,0],places=7)
        self.assertAlmostEqual(ref2/self.Da.data[7,0,0],self.D4d.data4D[1][7,0,0],places=7)
        
#================================================================================

    def test_divc(self):
        print "\nDoing Test divc"
        
        print "  - Check divc with copy=True on two levels and both type are data4D"
        r1 = self.D4d.divc(5.,copy=True)
        self.assertAlmostEqual(r1.data4D[0][4,0,0]*5.,self.D4d.data4D[0][4,0,0],places=7)
        self.assertAlmostEqual(r1.data4D[1][4,0,0]*5.,self.D4d.data4D[1][4,0,0],places=7)
        
        print "  - Check divc with copy=False on two levels and both type are data4D"
#        ref1 = self.D4d.data4D[0][5,0,0]
#        ref2 = self.D4d.data4D[1][5,0,0]
#        self.D4d.divc(666.,copy=False)
#        self.assertAlmostEqual(ref1+666.,self.D4d.data4D[0][5,0,0])
#        self.assertAlmostEqual(ref2+666.,self.D4d.data4D[1][5,0,0])
        
#================================================================================
    def test_setDataFromLevel(self):
        print "\nDoing Test setDataFromLevel"
        
	print "  - Check if Data  and Data4D are not equal at the start"  
        if self.D4d.data4D[0][17,0,0] == self.Da.data[17,0,0]:
	  print "  - Data and Data4D Elements are Equal. Please, try again"
	  
        print "  - Check if Data4d on level 1 and Data4D on level 2 are not equal at the start"  
        if self.D4d.data4D[0][17,0,0] == self.D4d.data4D[1][17,0,0]:
	  print "  - Data4D Elements on level 1 and 2 are Equal. Please, try again"
	  
	print "  - Check set Data to first level"  
	self.D4d.setDataFromLevel(self.Da,1)
        self.assertAlmostEqual(self.Da.data[17,0,0],self.D4d.data4D[0][17,0,0])
        
	print "  - Check set Data to second level"  
	self.D4d.setDataFromLevel(self.Da,2)
        self.assertAlmostEqual(self.Da.data[17,0,0],self.D4d.data4D[1][17,0,0])
        
	print "  - Check if Data on first level and second level equal"  
        self.assertAlmostEqual(self.D4d.data4D[0][17,0,0],self.D4d.data4D[1][17,0,0])

#================================================================================
 
    def test_getDataFromLevel(self):
        print "\nDoing Test getDataFromLevel"
        
        if self.D4d.data4D[0][33,0,0] == self.Da.data[33,0,0]:
	  print "  - Data and Data4D Elements are Equal. Please, try again"
	  
	print "  - Check getting Data from first level"  
	self.Da = self.D4d.getDataFromLevel(1)
        self.assertAlmostEqual(self.Da.data[33,0,0],self.D4d.data4D[0][33,0,0])
        
	print "  - Check getting Data from second level"  
	self.Da = self.D4d.getDataFromLevel(2)
        self.assertAlmostEqual(self.Da.data[33,0,0],self.D4d.data4D[1][33,0,0])
    
#================================================================================
 
    def test_copy(self):
        print "\nDoing Test copy"
        
        Da_1 = self.D4d.getDataFromLevel(1)
        Da_2 = self.D4d.getDataFromLevel(2)
        
        r_copy = self.D4d.copy()
        r_copy.mulc(2.,copy=False)
        
	print "  - Check if old values has not changed?"  
        self.assertAlmostEqual(Da_1.data[15,0,0],self.D4d.data4D[0][15,0,0])
        self.assertAlmostEqual(Da_2.data[14,0,0],self.D4d.data4D[1][14,0,0])
        
	print "  - Check if the values of the copy are correct?"  
        self.assertAlmostEqual(Da_1.data[29,0,0]*2.,r_copy.data4D[0][29,0,0])
        self.assertAlmostEqual(Da_2.data[24,0,0]*2.,r_copy.data4D[1][24,0,0])
       

#================================================================================
 
    def test_sum_data4D(self):
        print "\nDoing Test sum_data4D"
        
        Da_1 = self.D4d.getDataFromLevel(1)
        Da_2 = self.D4d.getDataFromLevel(2)
        
        D4d_sum = self.D4d.sum_data4D()
	print "  - Check if the sum across all levels are correct?"  
        self.assertAlmostEqual(Da_1.data[17,0,0]+Da_2.data[17,0,0],D4d_sum.data[17,0,0])
	

#================================================================================


if __name__ == '__main__':
    unittest.main()























