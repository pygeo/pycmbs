#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Walter Sauf"
__version__ = "0.1.4"
__date__ = "2013/04/16"
__email__ = "walter.sauf@zmaw.de"

'''
# Copyright (C) 2013 
# Alexander Loew, alexander.loew@mpimet.mpg.de
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
'''


'''
Data4D CLASS
class to work with the PFTs from JSBACH
The Data4D class has the Data Class as a base class. 
In addition to the Class Date has the Data4D class the attributs 

  from_level  : the number of the lowest level, starting by 0
  to_level    : the number of the highest level
  data4D      : is an array of the attribude 'data' from the Data-Class
  
becase all data are stored in the data4D array, the normal data attribude is deleted.

'''
from pylab import *

from data import Data
import matplotlib.colors as col
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from cdo import *


class Data4D(Data):

        
  
#-----------------------------------------------------------------------


    def read(self,shift_lon,start_time=None,stop_time=None,time_var='time',checklat=True):
        cdo = Cdo()

        if self.levellist == {}:
	  level = 0
	  L = cdo.showlevel(input='-selname,'+self.varname+' '+ self.filename )
	  levell = L[0].split(' ')
	  for k in levell:
            self.levellist[int(k)] = level
            level += 1
            
        self.level=int(self.levellist.get(self.levellist.keys()[0]))
        
        Data.read(self,shift_lon,start_time=start_time,stop_time=stop_time,time_var=time_var,checklat=checklat)

        
        del self.data
        for k in sorted(self.levellist.keys()):
          self.level=int(self.levellist[int(k)])
          self.data4D.append(Data.read_netcdf(self,self.varname))
        
          
#-----------------------------------------------------------------------

    def _copy_Data4D_Info_to_Data(self):
        """
        copies all Data-attribut information from the Data4D type to the Data type.
        @return C{Data} object
        """
        d = Data(None,None)


        for attr, value in self.__dict__.iteritems():
	    if attr != 'data4D' and attr != 'levellist':
#	      print " attr: " + attr + " value: "
              try:
                  #-copy (needed for arrays)
                  cmd = "d." + attr + " = self." + attr + '.copy()'
                  exec cmd
              except:
                  #-copy
                  cmd = "d." + attr + " = self." + attr
                  exec cmd

        return d

#-----------------------------------------------------------------------

    def copy(self):
        """
        copy complete C{Data4D} object including all attributes
        @return C{Data4D} object
        """
        d = Data4D(None,None)
#ws        print "copy of data4D.py"


        for attr, value in self.__dict__.iteritems():
#	  print attr
	  if not attr == "data4D":
#	    print 'if'
            try:
                #-copy (needed for arrays)
                cmd = "d." + attr + " = self." + attr + '.copy()'
                exec cmd
#                print "A: "+cmd
            except:
                #-copy
                cmd = "d." + attr + " = self." + attr
                exec cmd
#                print "B: "+cmd
          else:
	   #print 'else '+ str(self.levellist.keys())
	   #cmd = "d." + attr + " = self." + attr
           exec cmd
           for k in self.levellist.keys():
#	     print 'K: '+str(k)
	     d.data4D.append(self.data4D[self.levellist[k]].copy())

        return d


#-----------------------------------------------------------------------
    def getDataFromLevel(self,l):
        """
        returns a Data4D level as a Data object
        @param l level of data4D. Starting by 0:
        @return C{Data} object
        """
        if self.levellist.has_key(int(l)):
          ret = self._copy_Data4D_Info_to_Data()
          d = self
          #print int(l)
          #print self.levellist[int(l)]
          #print self.levellist
          ret.data = d.data4D[self.levellist[int(l)]]
          return ret
        else:
         raise ValueError, 'Given Level '+ str(l)+' not in data4D!'
	  

#-----------------------------------------------------------------------
    def setDataFromLevel(self,da,l):
        """
        saves the date from da into an data4D lavel. 
        The lebvel is given by parameter 'l'. 
        
        @param l level of Data4D. Starting by 0:
        @da  C{Data} object to be saved into Data4D
        """
      
        if self.levellist.has_key(int(l)):
          d = self
          d.data4D[self.levellist[int(l)]] = da.data
        else:
         raise ValueError, 'Given Level '+ str(l)+' not in data4D!'
          
        
#-----------------------------------------------------------------------
    def mulc(self,x,copy=True):
      
      if copy:
        d = self.copy()
      else:
        d = self

      for k in self.levellist.keys():
        d.data4D[self.levellist[k]] *= x
        
      return d
#-----------------------------------------------------------------------
    def mul(self,x,copy=True):
      
        if copy:
            d = self.copy()
        else:
            d = self

        for k in self.levellist.keys():
          if hasattr(x, 'data4D'):
	    d.data4D[self.levellist[k]] *= x.data4D[self.levellist[k]]
	  else:
	    d.data4D[self.levellist[k]] *= x.data
	  
#        d.label = self.label + ' * ' + x.label
        return d
#-----------------------------------------------------------------------
    def divc(self,x,copy=True):
      
      if copy:
          d = self.copy()
      else:
          d = self

      for k in self.levellist.keys():
        d.data4D[self.levellist[k]] /= x
        
      return d

#-----------------------------------------------------------------------
    def div(self,x,copy=True):
      
        if copy:
            d = self.copy()
        else:
            d = self

        for k in self.levellist.keys():
          if hasattr(x, 'data4D'):
	    d.data4D[self.levellist[k]] /= x.data4D[self.levellist[k]]
	  else:
	    d.data4D[self.levellist[k]] /= x.data
          
        return d

#-----------------------------------------------------------------------
    def addc(self,x,copy=True):
      
      if copy:
          d = self.copy()
      else:
          d = self

      for k in self.levellist.keys():
        d.data4D[self.levellist[k]] += x
        
      return d
#-----------------------------------------------------------------------
    def add(self,x,copy=True):
      
        if copy:
            d = self.copy()
        else:
            d = self
            
        #print "01 copy self "+str(self.data4D[0][4,0,0])+" d "+str(d.data4D[0][4,0,0])
        #print "a1 copy self "+self.label+" d "+str(d.label)
        #print "a1 copy self "+self.label+" d "+str(d.label)
#        d.data4D[0][4,0,0] = d.data4D[0][4,0,0] + 100.
        d.label = "myLabel"
#        print "a2 copy self "+self.label+" d "+str(d.label)
        for k in self.levellist.keys():
#          print "02 copy self "+str(self.data4D[0][4,0,0])+" d "+str(d.data4D[0][4,0,0])
          if hasattr(x, 'data4D'):
	    d.data4D[self.levellist[k]] += x.data4D[self.levellist[k]]
	  else:
	    d.data4D[self.levellist[k]] += x.data
          
#        print "03 copy self "+str(self.data4D[0][4,0,0])+" d "+str(d.data4D[0][4,0,0])
        return d

#-----------------------------------------------------------------------
    def subc(self,x,copy=True):
      
      if copy:
          d = self.copy()
      else:
          d = self

      for k in self.levellist.keys():
        d.data4D[self.levellist[k]] -= x
        
      return d
#-----------------------------------------------------------------------
    def sub(self,x,copy=True):
      
        if copy:
            d = self.copy()
        else:
            d = self

        for k in self.levellist.keys():
          if hasattr(x, 'data4D'):
	    d.data4D[self.levellist[k]] -= x.data4D[self.levellist[k]]
	  else:
	    d.data4D[self.levellist[k]] -= x.data
          
        return d

#-----------------------------------------------------------------------
    def sum_data4D(self):
        """
        summates all Data4D levels to one level
        @return C{Data} object
        """
      
        sum = self._copy_Data4D_Info_to_Data()
        d = self
            
        sum.data = 0.0
        
        for k in self.levellist.keys():
          sum.data += self.data4D[self.levellist[k]]
          
          
        return sum
#-----------------------------------------------------------------------------------------------------------------------

    def __init__(self,filename,varname,levellist=None,**kwargs):
        """
        Data4D class

        This class implements the functionality to generate Data4D objekt.


        EXAMPLES
        ========

        """
        self.levellist = {}
        
        if not levellist == None:
	  level = 0
	  for k in levellist:
            self.levellist[int(k)] = level
            level += 1
	  
        self.data4D = []
        
        Data.__init__(self,filename,varname, **kwargs)
        self.mulc(self.scale_factor,copy=False)
    