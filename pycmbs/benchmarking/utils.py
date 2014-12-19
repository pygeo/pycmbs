# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os
from pycmbs.data import Data
from cdo import Cdo
import numpy as np
import warnings


def get_data_pool_directory():
    """
    get data pool directory for /pool/SEP
    """
    current_dir = os.getcwd() + os.sep

    if 'SEP' in os.environ.keys():
        data_pool_directory = os.environ['SEP']

        # test if the data pool directory variable
        # contains multiple pathes
        found = False
        if os.sep == '/':  # linux
            if ':' in data_pool_directory:
                for d in data_pool_directory.split(':'):
                    if os.path.exists(d) and not found:
                        found = True
                        data_pool_directory = d
                if not found:
                    raise ValueError('The data pool directory variable \
                           SEP does not contain a valid pathname!')
        else:
            raise ValueError('Currently this routine is only tested \
                            for Linux!')
    else:
        data_pool_directory = '/pool/SEP/'

    if not os.path.exists(data_pool_directory):
        warnings.warn('Data pool directory %s not found, \
                using current directory: %s' % (data_pool_directory, current_dir), Warning)
        data_pool_directory = current_dir

    return data_pool_directory


def get_temporary_directory():
    """
    returns temporary directory

    Returns
    -------
    r : str
    """
    if 'CDOTEMPDIR' in os.environ.keys():
        tempdir = os.environ['CDOTEMPDIR']
    else:
        tempdir = '.' + os.sep
    if tempdir[-1] != os.sep:
        tempdir += os.sep
    return tempdir


def get_generic_landseamask(shift_lon, mask_antarctica=True,
                            area='land', interpolation_method='remapnn',
                            target_grid='t63grid', force=False):
    """
    get generic land/sea mask. The routine uses the CDO command 'topo'
    to generate a 0.5 degree land/sea mask and remaps this
    using nearest neighbor
    to the target grid

    NOTE: using inconsistent land/sea masks between datasets can
    result in considerable biases. Note also that
    the application of l/s mask is dependent on the spatial resolution

    This routine implements a VERY simple approach, but assuming
    that all areas >0 m height are land and the rest is ocean.

    Parameters
    ----------
    shift_lon : bool
        specifies if longitudes shall be shifted
    interpolation_method : str
        specifies the interpolation method
        that shall be used for remapping the 0.5degree data
        to the target grid. This can be any of ['remapnn','remapcon',
        'remapbil']
    target_grid : str
        specifies target grid to interpolate to as
        similar to CDO remap functions. This can be either a string or
        a filename which includes valid geometry information
    force : bool
        force calculation (removes previous file) = slower

    area : str
        ['land','ocean']. When 'land', then the mask returned
        is True on land pixels, for ocean it is vice versa.
        in any other case, you get a valid field everywhere
        (globally)

    mask_antarctica : bool
        mask antarctica; if True, then the mask is
        FALSE over Antarctice (<60S)

    Returns
    -------
    returns a Data object
    """

    print ('WARNING: Automatic generation of land/sea mask. \
            Ensure that this is what you want!')

    cdo = Cdo()

    #/// construct output filename.
    #If a filename was given for the grid, replace path separators ///
    target_grid1 = target_grid.replace(os.sep, '_')
    outputfile = get_temporary_directory() + 'land_sea_fractions_' \
        + interpolation_method + '_' + target_grid1 + '.nc'

    print 'outfile: ', outputfile
    print 'cmd: ', '-remapnn,' + target_grid + ' -topo'

    #/// interpolate data to grid using CDO ///
    cdo.monmean(options='-f nc', output=outputfile,
                input='-remapnn,' + target_grid + ' -topo', force=force)

    #/// generate L/S mask from topography (land = height > 0.
    ls_mask = Data(outputfile, 'topo', read=True,
                   label='generic land-sea mask',
                   lat_name='lat', lon_name='lon',
                   shift_lon=shift_lon)
    print('Land/sea mask can be found on file: %s' % outputfile)

    if area == 'land':
        msk = ls_mask.data > 0.  # gives land
    elif area == 'ocean':
        msk = ls_mask.data <= 0.
    else:
        msk = np.ones(ls_mask.data.shape).astype('bool')
    ls_mask.data[~msk] = 0.
    ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool')

    #/// mask Antarctica if desired ///
    if mask_antarctica:
        ls_mask.data[ls_mask.lat < -60.] = False

    # ensure that also the mask attribute is set properly
    ls_mask._apply_mask(~msk)

    return ls_mask


def get_T63_landseamask(shift_lon, mask_antarctica=True, area='land'):
    """
    get JSBACH T63 land sea mask
    the LS mask is read from the JSBACH init file

    area : str
        ['land','ocean']: When 'land', then the mask returned
        is True on land pixels, for ocean it is vice versa.
        In any other case, you get a valid field everywhere (globally)

    mask_antarctica : bool
        if True, then the mask is FALSE over Antarctica (<60S)
    """
    ls_file = get_data_pool_directory() \
        + 'data_sources/LSMASK/jsbach_T63_GR15_4tiles_1992.nc'
    ls_mask = Data(ls_file, 'slm', read=True, label='T63 land-sea mask',
                   lat_name='lat', lon_name='lon', shift_lon=shift_lon)
    if area == 'land':
        msk = ls_mask.data > 0.
    elif area == 'ocean':
        msk = ls_mask.data == 0.
    else:
        msk = np.ones(ls_mask.data.shape).astype('bool')

    ls_mask.data[~msk] = 0.
    ls_mask.data[msk] = 1.
    ls_mask.data = ls_mask.data.astype('bool')
    if mask_antarctica:
        ls_mask.data[ls_mask.lat < -60.] = False

    # ensure that also the mask attribute is set properly
    ls_mask._apply_mask(~msk)

    return ls_mask
