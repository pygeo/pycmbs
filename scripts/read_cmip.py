# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
test reading CMIP single data
"""

from pycmbs.benchmarking.models import CMIP5RAW_SINGLE
import datetime
from pycmbs.benchmarking import config

import matplotlib.pyplot as plt



plt.close('all')

model_dict = {'rain': {'CMIP5':
                       {
                           'variable': 'pr',
                           'unit': 'mm/day',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'model_prefix': 'Amon',
                           'file_format': 'nc',
                           'scale_factor': 86400.,
                           'valid_mask': 'ocean'
                       },


                       'JSBACH_RAW2':
                       {
                           'variable': 'precip_acc',
                           'unit': 'mm/day',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'file_format': 'nc',
                           'scale_factor': 86400.,
                           'valid_mask': 'global'
                       }
                       },


              'evap': {'CMIP5':
                       {
                           'variable': 'evspsbl',
                           'unit': 'mm/day',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'file_format': 'nc',
                           'model_prefix': 'Amon',
                           'scale_factor': 86400.,
                           'valid_mask': 'ocean'
                       }
                       },

              'twpa': {'CMIP5':
                       {
                           'variable': 'clwvi',
                           'unit': 'kg/m^2',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'file_format': 'nc',
                           'model_prefix': 'Amon',
                           'scale_factor': 1.,
                           'valid_mask': 'ocean'
                       }
                       },

              'wind': {'CMIP5':
                      {
                          'variable': 'sfcWind',
                          'unit': 'm/s',
                          'lat_name': 'lat',
                          'lon_name': 'lon',
                          'model_suffix': 'ensmean',
                          'file_format': 'nc',
                          'model_prefix': 'Amon',
                          'scale_factor': 1.,
                          'valid_mask': 'ocean'
                      }
              },

              'wvpa': {'CMIP5':
                       {
                           'variable': 'prw',
                           'unit': 'kg m^2',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'file_format': 'nc',
                           'model_prefix': 'Amon',
                           'scale_factor': 1,
                           'valid_mask': 'ocean'
                       }
                       },

              'late': {'CMIP5':
                       {
                           'variable': 'hfls',
                           'unit': 'W/m^2',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'file_format': 'nc',
                           'model_prefix': 'Amon',
                           'scale_factor': 1,
                           'valid_mask': 'ocean'
                       }
                       },

              'hair': {'CMIP5':
                       {
                           'variable': 'huss',
                           'unit': '$kg/kg^2$',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'file_format': 'nc',
                           'model_prefix': 'Amon',
                           'scale_factor': 1,
                           'valid_mask': 'ocean'
                       }
                       },

              'seaice_concentration': {'CMIP5':
                                       {
                                           'variable': 'sic',
                                           'unit': '-',
                                           'lat_name': 'lat',
                                           'lon_name': 'lon',
                                           'model_suffix': 'ens_mean_185001-200512',
                                           'file_format': 'nc',
                                           'model_prefix': 'OImon',
                                           'scale_factor': 1,
                                           'valid_mask': 'ocean',
                                           'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk'
                                       },

                                       'CMIP3':
                                       {
                                           'variable': 'SICOMO',
                                           'unit': '-',
                                           'lat_name': 'lat',
                                           'lon_name': 'lon',
                                           'model_suffix': '1860-2100.ext',
                                           'file_format': 'nc',
                                           'model_prefix': '',
                                           'scale_factor': 100.,
                                           'valid_mask': 'ocean',
                                           'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk',
                                           'level': 0
                                       },
                                       },

              'seaice_extent': {'CMIP5':
                               {
                                'variable': 'sic',
                                'unit': '-',
                                'lat_name': 'lat',
                                'lon_name': 'lon',
                                'model_suffix': 'ens_mean_185001-200512',
                                'file_format': 'nc',
                                'model_prefix': 'OImon',
                                'scale_factor': 1,
                                'valid_mask': 'ocean',
                                'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk'
                                },

                                'CMIP3':
                               {
                                'variable': 'SICOMO',
                                   'unit': '-',
                                   'lat_name': 'lat',
                                   'lon_name': 'lon',
                                   'model_suffix': '1860-2100.ext',
                                   'file_format': 'nc',
                                   'model_prefix': '',
                                   'scale_factor': 100.,
                                   'valid_mask': 'ocean',
                                   'custom_path': '/home/m300028/shared/dev/svn/pyCMBS/dirk',
                                   'level': 0
                               },
              },

              'budg': {'CMIP5':
                       {
                           'variable': 'budg',
                           'unit': 'mm/d',
                           'lat_name': 'lat',
                           'lon_name': 'lon',
                           'model_suffix': 'ensmean',
                           'file_format': 'nc',
                           'model_prefix': 'Amon',
                           'scale_factor': 86400.,
                           'valid_mask': 'ocean',
                           'custom_path': '/net/nas2/export/eo/workspace/m300036/pycmbs-cmsaf/data'
                       }
                       },


              'sis': {'JSBACH_RAW2':
                      {
                          'variable': 'swdown_acc',
                          'unit': '$W/m^2$',
                          'lat_name': 'lat',
                          'lon_name': 'lon',
                          'file_format': 'nc',
                          'scale_factor': 1.,
                          'valid_mask': 'land'
                      },
                      'CMIP5' :
                      {
                          'valid_mask': 'land'
                      },
                      'CMIP5RAWSINGLE' :
                      {
                          'valid_mask': 'land'
                      }

                      },

              'surface_upward_flux': {'JSBACH_RAW2':
                                      {
                                          'variable': 'swdown_reflect_acc',
                                          'unit': '$W/m^2$',
                                          'lat_name': 'lat',
                                          'lon_name': 'lon',
                                          'file_format': 'nc',
                                          'scale_factor': 1.,
                                          'valid_mask': 'land'
                                      }
                                      },

              'albedo_vis': {'JSBACH_RAW2':
                             {
                                 'variable': 'var14',
                                 'unit': '-',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'file_format': 'nc',
                                 'scale_factor': 1.,
                                 'valid_mask': 'land'
                             }
                             },

              'albedo_nir': {'JSBACH_RAW2':
                             {
                                 'variable': 'var15',
                                 'unit': '-',
                                 'lat_name': 'lat',
                                 'lon_name': 'lon',
                                 'file_format': 'nc',
                                 'scale_factor': 1.,
                                 'valid_mask': 'land'
                             }
                             },

              'temperature': {
                  'JSBACH_RAW2':
                  {
                      'variable': 'temp2',
                      'unit': 'K',
                      'lat_name': 'lat',
                      'lon_name': 'lon',
                      'file_format': 'nc',
                      'scale_factor': 1.,
                      'valid_mask': 'global'
                  }
              }

              }



data_dir = '/home/m300028/shared/private/geo2consult/WPbench/work/testdata/miklip/'
model = 'MPI-M:MPI-ESM-LR#1'
experiment = 'amip'

start_time = datetime.datetime(2000,1,1)
stop_time = datetime.datetime(2007,9,30)

file = '/home/m300028/shared/private/geo2consult/WPbench/work/dwd.cfg'
CF = config.ConfigFile(file)
variables = CF.variables
varmethods = CF.get_methods4variables(variables, model_dict)

themodel = CMIP5RAW_SINGLE(data_dir, model, experiment, varmethods,
                        intervals=CF.intervals, lat_name='lat',
                        lon_name='lon', label=model,
                        start_time=start_time,
                        stop_time=stop_time,
                        shift_lon=True)


PCFG = config.PlotOptions()
PCFG.read(CF)
plot_options = PCFG
themodel.plot_options = plot_options
themodel.get_data()




x=themodel.variables['sis']

f = plt.figure()
ax=f.add_subplot(111)
ax.imshow(x.lat)






