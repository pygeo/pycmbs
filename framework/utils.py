#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.1"
__date__ = "2012/10/29"
__email__ = "alexander.loew@zmaw.de"

'''
# Copyright (C) 2012 Alexander Loew, alexander.loew@zmaw.de
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

from pyCMBS import *

import os

def get_data_pool_directory():
    if 'SEP' in os.environ.keys():
        data_pool_directory = os.environ['SEP'] #get directory of pool/SEP
    else:
        data_pool_directory = '/pool/SEP/'

    return data_pool_directory

def get_temporary_directory():
    """
    @return: path of temporary directory
    @rtype: str
    """
    if 'CDOTEMPDIR' in os.environ.keys():
        tempdir = os.environ['CDOTEMPDIR']
    else:
        tempdir = './'
    if tempdir[-1] != '/':
        tempdir = tempdir + '/'

    return tempdir



def get_T63_landseamask(shift_lon,mask_antarctica=True):
    """
    get JSBACH T63 land sea mask
    the LS mask is read from the JSBACH init file

    @todo: put this to the JSBACH model class
    """
    ls_file = get_data_pool_directory() + 'variables/land/land_sea_mask/jsbach_T63_GR15_4tiles_1992.nc'
    ls_mask = Data(ls_file,'slm',read=True,label='T63 land-sea mask',lat_name='lat',lon_name='lon',shift_lon=shift_lon)
    msk=ls_mask.data>0.; ls_mask.data[~msk] = 0.; ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool') #convert to bool
    if mask_antarctica:
        ls_mask.data[ls_mask.lat < -60.] = False

    return ls_mask

def get_T63_weights(shift_lon):
    """
    get JSBACH T63 cell weights

    @todo: put this to the JSBACH model class
    """
    w_file = get_data_pool_directory() + 'variables/land/land_sea_mask/t63_weights.nc'
    weight = Data(w_file,'cell_weights',read=True,label='T63 cell weights',lat_name='lat',lon_name='lon',shift_lon=shift_lon)

    return weight.data
