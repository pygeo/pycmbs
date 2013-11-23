__author__ = "Alexander Loew"
__version__ = "0.1.4"
__date__ = "2013/06/30"
__email__ = "alexander.loew@mpimet.mpg.de"

"""
Copyright (C) 2012 Alexander Loew, alexander.loew@mpimet.mpg.de
See COPYING file for copying and redistribution conditions.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
"""


from pyCMBS.data import Data
import os
from netcdf import *
import numpy as np


class Icon(Data):
    """
    Main class for ICON data handling
    """

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
        File=NetCDFHandler()
        File.open_file(self.gridfile,'r')
        self.lon = File.get_variable('clon')*180./np.pi #grid cell center coordinates
        self.lat = File.get_variable('clat')*180./np.pi
        self.ncell = len(self.lon)

        self.vlon = File.get_variable('clon_vertices')*180./np.pi
        self.vlat = File.get_variable('clat_vertices')*180./np.pi

        File.close()

        #--- read time variable
        if self.time_var != None:
            self.time = self.read_netcdf(self.time_var) #returns either None or a masked array
            if hasattr(self.time,'mask'):
                self.time = self.time.data
            else:
                self.time = None
            if self.time != None:
                if self.time.ndim != 1:
                    self.time = self.time.flatten() #remove singletone dimensions
        else:
            self.time = None

        #--- determine time --> convert to python timestep
        if self.time != None:
            self.set_time()





