# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""

from pyCMBS.data import Data
import os
from pycmbs.netcdf import *
import numpy as np


class Icon(Data):
    """
    Main class for ICON data handling
    """

    def __init__(self, filename, gridfile, varname, read=False, **kwargs):
        """
        Parameters
        ----------

        filename : str
            filename of data file

        gridfile : str
            filename of grid definition file

        varname : str
            name of variable to handle

        read : bool
            specify if data should be read immediately
        """
        Data.__init__(self, filename, varname, **kwargs)
        self.gridfile = gridfile
        self.gridtype = 'unstructured'

#---

    def read(self, time_var='time'):
        """
        This is a special routine for reading data from ICON structure
        a bit redundant to Data.read()

        Parameters
        ----------

        time_var : str
            name of time variable (default='time')
        """
        print('Reading ICON data ...')

        if not os.path.exists(self.filename):
            raise ValueError('File not existing: %s' % self.filename)
        if not os.path.exists(self.gridfile):
            raise ValueError('File not existing: %s' % self.gridfile)

        #--- time variable
        self.time_var = time_var

        #--- data field
        # [time,ncell]
        self.data = self.read_netcdf(self.varname)
        nt, ncell = self.data.shape
        # reshape so we have a common 3D structure like always in pyCMBS
        self.data = self.data.reshape((nt, 1, ncell))
        if self.data is None:
            raise ValueError('The data in the file %s is not existing. \
                            This must not happen!' % self.filename)
        if self.scale_factor is None:
            raise ValueError('The scale_factor for file %s is NONE, \
                              this must not happen!' % self.filename)

        self.data *= self.scale_factor

        #--- read lat/lon
        File = NetCDFHandler()
        File.open_file(self.gridfile, 'r')
        # grid cell center coordinates
        self.lon = File.get_variable('clon') * 180. / np.pi
        self.lat = File.get_variable('clat') * 180. / np.pi
        self.ncell = len(self.lon)

        self.vlon = File.get_variable('clon_vertices') * 180. / np.pi
        self.vlat = File.get_variable('clat_vertices') * 180. / np.pi

        File.close()

        #--- read time variable
        if self.time_var is not None:
            # returns either None or a masked array
            self.time = self.read_netcdf(self.time_var)
            if hasattr(self.time, 'mask'):
                self.time = self.time.data
            else:
                self.time is None
            if self.time is not None:
                if self.time.ndim != 1:
                    # remove singletone dimensions
                    self.time = self.time.flatten()
        else:
            self.time = None

        #--- determine time --> convert to python timestep
        if self.time is not None:
            self.set_time()
