#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyCMBS import *

import os

def get_data_pool_directory():
    if 'SEP' in os.environ.keys():
        data_pool_directory = os.environ['SEP'] #get directory of pool/SEP
    else:
        data_pool_directory = '/pool/SEP/'

    return data_pool_directory



def get_T63_landseamask(shift_lon):
    '''
    get JSBACH T63 land sea mask
    the LS mask is read from the JSBACH init file

    @todo: put this to the JSBACH model class
    '''
    ls_file = get_data_pool_directory() + 'variables/land/land_sea_mask/jsbach_T63_GR15_4tiles_1992.nc'
    ls_mask = Data(ls_file,'slm',read=True,label='T63 land-sea mask',lat_name='lat',lon_name='lon',shift_lon=shift_lon)
    msk=ls_mask.data>0.; ls_mask.data[~msk] = 0.; ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool') #convert to bool

    return ls_mask

def get_T63_weights(shift_lon):
    '''
    get JSBACH T63 cell weights

    @todo: put this to the JSBACH model class
    '''
    w_file = get_data_pool_directory() + 'variables/land/land_sea_mask/t63_weights.nc'
    weight = Data(w_file,'cell_weights',read=True,label='T63 cell weights',lat_name='lat',lon_name='lon',shift_lon=shift_lon)

    return weight.data
