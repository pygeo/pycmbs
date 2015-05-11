# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""
"""
This module allows a flexible choice of the netCDF backend
"""

import os

valid_backends = ['netCDF4']


class NetCDFHandler(object):

    def __init__(self, netcdf_backend='netCDF4'):
        """
        general NetCDF handler
        Currently two different backends are supported

        netcdf_backend : str
            ['netCDF4'] only netCDF4 can be used so far
            backends

        """
        self.type = netcdf_backend

        if self.type not in valid_backends:
            raise ValueError('Invalid data backend!')

        if self.type.lower() == 'netcdf4':
            import netCDF4 as Cdf
        else:
            raise ValueError('Invalid netCDF backend!')
        self.handler = Cdf

    def open_file(self, filename, mode, format='NETCDF4'):
        """
        Open a netCDF file using the predefined backend

        Parameters
        ----------
        filename : str
            name of file to read
        mode : str
            specify read or write data access ['w','r']
        format : str
            output file format specifieralue
            ['NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC']

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

        if self.type.lower() == 'netcdf4':
            if mode == 'r':
                self.F = self.handler.Dataset(filename, mode=mode)
            elif mode == 'w':
                self.F = self.handler.Dataset(filename, mode=mode,
                                              format=format)  # TODO check format
            self.create_dimension = self.F.createDimension
            self.create_variables = self.F.createVariable
        else:
            print self.type, format
            raise ValueError('Something went wrong!')

    def get_variable_keys(self):
        """
        return a list of variable names in the file
        """
        if self.type.lower() == 'netcdf4':
            return self.F.variables.keys()
        else:
            raise ValueError('Something went wrong!')

    def _get_long_name(self, varname):
        """
        get longname of variable

        Parameters
        ----------
        varname : str
            variable name
        """
        var = self.get_variable_handler(varname)
        self.long_name = '-'
        if self.type.lower() == 'netcdf4':
            if hasattr(var, 'long_name'):
                self.long_name = var.long_name
        else:
            raise ValueError('Invalid backend!')

    def _get_unit(self, varname):
        """
        get unit of variable

        Parameters
        ----------
        varname : str
            variable name
        """
        var = self.get_variable_handler(varname)
        if self.type.lower() == 'netcdf4':
            if hasattr(var, 'units'):
                return var.units
        else:
            raise ValueError('Invalid backend!')

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
        if self.type.lower() == 'netcdf4':
            return self.F.variables[varname][:].astype('float').copy()
        else:
            raise ValueError('Something went wrong!')

    def set_attribute(self, varname, key, value):
        """
        set attributes for a field

        Parameters
        ----------
        varname : str
            name of variable
        key : str
            attribute key
        value : str/float/int
            attribute value
        """
        V = self.get_variable_handler(varname)
        if self.type.lower() == 'netcdf4':
            V.setncattr(key, value)
        else:
            raise ValueError('Something went wrong!')

    def get_variable_handler(self, varname):
        """
        Get handler to a variable
        """

        if self.type.lower() == 'netcdf4':
            return self.F.variables[varname]
        else:
            raise ValueError('Something went wrong!')

    def _get_scale_factor(self, varname):
        var = self.get_variable_handler(varname)
        if self.type.lower() == 'netcdf4':
            # netCDF4 library already applies the scaling factor!
            return 1.
        else:
            raise ValueError('Something went wrong!')

    def _get_fill_value(self, varname):
        """ return FillValue for variable if given """
        var = self.get_variable_handler(varname)
        if self.type.lower() == 'netcdf4':
            if hasattr(var, '_FillValue'):
                return float(var._FillValue)
            else:
                return None
        else:
            raise ValueError('Invalid backend!')

    def _get_add_offset(self, varname):
        var = self.get_variable_handler(varname)
        if self.type.lower() == 'netcdf4':
            # netCDF4 library already applies the add_offset!
            return 0.
        else:
            raise ValueError('Something went wrong!')

    def _get_time(self):
        att = self.F.get('Metadata').attrs
        if 'Date&Time' in att.keys():
            return att['Date&Time']
        else:
            return None

    def assign_value(self, varname, value):
        """
        assign a value to a variable to be written to a netCDF file
        """
        if self.type.lower() == 'netcdf4':
            if value.ndim == 1:
                self.F.variables[varname][:] = value[:]
            elif value.ndim == 2:
                self.F.variables[varname][:, :] = value[:, :]
            elif value.ndim == 3:
                self.F.variables[varname][:, :, :] = value[:, :, :]
            else:
                raise ValueError('Unsupported dimension!')
        else:
            raise ValueError('Something went wrong!')

    def create_dimension(self, dimname, size=None):
        """
        creates a new dimension in netCDF file

        dimname : str
            name of dimension
        size : int
            size of dimension; if None, then it is unlimited
        """
        if size is not None:
            assert(isinstance(size, int))

        if self.type.lower() == 'netcdf4':
            self.F.createDimension(dimname, size=size)
        else:
            raise ValueError('Something went wrong!')

    def create_variable(self, varname, dtype, dim, complevel=6, zlib=True, fill_value=None):
        """
        create a new variable in a netCDF file

        varname : str
            name of variable
        dtype : str
            datatype of variable
        dim : tuple
            tuple specifying the dimensions of the variables
        zlib : bool
            compress outout using zlib
        complevel : int
            compression level
        fill_value : float
            fill value for data
        """

        if self.type.lower() == 'netcdf4':
            if fill_value is not None:
                self.F.createVariable(varname, dtype, dimensions=dim,
                                      fill_value=fill_value, zlib=zlib, complevel=complevel)
            else:
                self.F.createVariable(
                    varname, dtype, dimensions=dim, zlib=zlib, complevel=complevel)
        else:
            raise ValueError('Something went wrong!')

    def close(self):
        self.F.close()
