"""
This module allows a flexible choice of the netCDF backend
"""

netcdf_backend='Nio'

class NetCDFHandler(object):
    def __init__(self):
        self.type = netcdf_backend
        if self.type.lower() == 'nio':
            import Nio as Cdf
        elif self.type.lower() == 'netcdf4':
            import netCDF4 as Cdf
        else:
            raise ValueError, 'Invalid netCDF backend!'
        self.handler = Cdf

    def open_file(self,filename,mode):
        """
        Returns
        -------
        F : file handler
            returns a file handler
        """
        if mode not in ['w','r']:
            raise ValueError, 'Invalid mode!'
        if self.type.lower() == 'nio':
            self.F = self.handler.open_file(filename,mode=mode)
        else:
            self.F = self.handler.Dataset(filename,mode=mode)


    def get_variable(self,varname):
        """
        Get data for a particular variable

        Returns
        -------
        data : ndarray
            returns data as a 2D,3D numpy array
        """
        if self.type.lower() == 'nio':
            return self.F.variables[varname].get_value().astype('float').copy()
        else:
            raise ValueError, 'TODO still needs to be implemented!'
            data = var[:].copy()

    def get_variable_handler(self,varname):
        """
        Get handler to a variable

        Returns
        -------

        """
        if self.type.lower() == 'nio':
            return self.F.variables[varname]
        else:
            raise ValueError, 'TODO still needs to be implemented!'


    def close(self):
        self.F.close()



