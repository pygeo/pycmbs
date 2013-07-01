__author__ = "Alexander Loew"
__version__ = "0.1.1"
__date__ = "2013/06/30"
__email__ = "alexander.loew@zmaw.de"

"""
Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
See COPYING file for copying and redistribution conditions.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
"""

"""
Main class for ICON data handling
"""
from pyCMBS.data import Data
import os
import Nio
import numpy as np


class Icon(Data):
    def __init__(self,filename,gridfile,varname,read=False,**kwargs):
        """
        @param filename: filename of data file
        @type filename : str

        @param gridfile: filename of grid definition file
        @type gridfile: str

        @param varname: name of variable to handle
        @type varname: str
        """
        Data.__init__(self,filename,varname,**kwargs)
        self.gridfile = gridfile

        self.gridtype = 'unstructured'

#---

    def read(self,time_var='time'):
        """
        This is a special routine for reading data from ICON structure
        a bit redundant to Data.read()
        """
        print 'Reading ICON data ...'

        if not os.path.exists(self.filename):
            raise ValueError, 'File not existing: ' + self.filename
        if not os.path.exists(self.gridfile):
            raise ValueError, 'File not existing: ' + self.gridfile

        #--- time variable
        self.time_var = time_var

        #--- data field
        self.data = self.read_netcdf(self.varname) #[time,ncell]
        nt,ncell = self.data.shape
        self.data = self.data.reshape((nt,1,ncell)) #reshape so we have a common 3D structure like always in pyCMBS
        if self.data == None:
            raise ValueError, 'The data in the file ' + self.filename + ' is not existing. This must not happen!'
        if self.scale_factor == None:
            raise ValueError, 'The scale_factor for file ' + self.filename + 'is NONE, this must not happen!'

        self.data = self.data * self.scale_factor

        #--- read lat/lon
        F=Nio.open_file(self.gridfile,'r')
        self.lon = F.variables['clon'].get_value()*180./np.pi #grid cell center coordinates
        self.lat = F.variables['clat'].get_value()*180./np.pi
        self.ncell = len(self.lon)

        self.vlon = F.variables['clon_vertices'].get_value()*180./np.pi #vertex coordinates [ncell,3]
        self.vlat = F.variables['clat_vertices'].get_value()*180./np.pi
        F.close()

        print 'ICON data was read!'








