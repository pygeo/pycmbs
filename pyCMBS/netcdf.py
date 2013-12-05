# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""

"""
This module allows a flexible choice of the netCDF backend
"""
# specify backend type
netcdf_backend = 'netCDF4'
#~ netcdf_backend = 'Nio'

import os


class NetCDFHandler(object):
    def __init__(self):
        self.type = netcdf_backend
        if self.type.lower() == 'nio':
            import Nio as Cdf
        elif self.type.lower() == 'netcdf4':
            import netCDF4 as Cdf
        else:
            raise ValueError('Invalid netCDF backend!')
        self.handler = Cdf

    def open_file(self, filename, mode):
        """
        Open a netCDF file using the predefined backend

        Parameters
        ----------
        filename : str
            name of file to read
        mode : str
            specify read or write data access ['w','r']

        Returns
        -------
        F : file handler
            returns a file handler
        """
        if mode not in ['w', 'r']:
            raise ValueError('ERROR: Invalid mode! [w,r], %s' % mode)
        if mode == 'r':
            if not os.path.exists(filename):
                raise ValueError('ERROR: File not existing: %s' % filename)
        elif mode == 'w':
            if len(os.path.dirname(filename)) > 0:
                if not os.path.exists(os.path.dirname(filename)):
                    os.makedirs(os.path.dirname(filename))
                if not os.path.exists(os.path.dirname(filename)):
                    raise ValueError('ERROR: Directory for output \
                    file %s was not existing and could also \
                    not be generated!' % filename)

        if self.type.lower() == 'nio':
            self.F = self.handler.open_file(filename, mode=mode)
            self.create_dimension = self.F.create_dimension
            self.create_variable = self.F.create_variable
        else:
            self.F = self.handler.Dataset(filename, mode=mode)
            self.create_dimension = self.F.createDimension
            self.create_variables = self.F.createVariable

    def get_variable(self, varname):
        """
        Get data for a particular variable

        Parameters
        ----------
        varname : str
            variable name of the netcdf variable to read

        Returns
        -------
        data : ndarray
            returns data as a 2D,3D numpy array
        """
        if self.type.lower() == 'nio':
            return self.F.variables[varname].get_value().astype('float').copy()
        else:
            return self.F.variables[varname][:].astype('float').copy()

    def get_variable_handler(self, varname):
        """
        Get handler to a variable

        Returns
        -------

        """
        if self.type.lower() == 'nio':
            return self.F.variables[varname]
        else:
            return self.F.variables[varname]

    def _get_scale_factor(self, varname):
        if self.type.lower() == 'nio':
            try:
                return float(self.F.variables[varname].scale_factor)
            except:
                print('No scale factor for variable! % s' % varname)
                return 1.
        else:
            # netCDF4 library already applies the scaling factor!
            return 1.

    def _get_add_offset(self, varname):
        if self.type.lower() == 'nio':
            try:
                return float(self.F.variables[varname].add_offset)
            except:
                print('No offset for variable! % s' % varname)
                return 0.
        else:
            # netCDF4 library already applies the add_offset!
            return 0.

    def assign_value(self, varname, value):
        """
        assign a value to a variable to be written to a netCDF file
        """
        if self.type.lower() == 'nio':
            self.F.variables[varname].assign_value(value)
        else:
            self.F.variables[varname][:] = value[:]

    def create_dimension(self, dimname, size=None):
        """
        creates a new dimension in netCDF file

        dimname : str
            name of dimension
        size : int
            size of dimension; if None, then it is unlimited
        """
        if size is not None:
            assert(isinstance(size,int))

        if self.type.lower() == 'nio':
            self.F.create_dimension(dimname, size)
        else:
            self.F.createDimension(dimname, size=size)

    def create_variable(self, varname, dtype, dim):
        """
        create a new variable in a netCDF file

        varname : str
            name of variable
        dtype : str
            datatype of variable
        dim : tuple
            tuple specifying the dimensions of the variables
        """

        if self.type.lower() == 'nio':
            self.F.create_variable(varname, dtype, dim)
        else:
            self.F.createVariable(varname, dtype, dimensions=dim)

    def close(self):
        self.F.close()
