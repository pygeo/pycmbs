# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

import os
from pycmbs.data import Data
import tempfile

def get_example_directory():
    """ returns directory where this file is located """

    testfile = 'xxxxdownload_test.txt'

    # by default use examples directory
    r = os.path.dirname(os.path.realpath(__file__))

    # check if one can write into the directory
    try:
        f = open(r + os.sep + testfile, 'w')
        f.write('test')
        f.close()
        os.remove(r + os.sep + testfile)
    except:
        # if no write access then
        r = tempfile.mkdtemp()

    if r[-1] != os.sep:
        r += os.sep

    return r

def get_example_data_directory():
    """ returns directory where the example data should be """
    return get_example_directory() + 'example_data' + os.sep

def get_sample_file(name='air', return_object=True):
    """
    returns Data object of example file including or the filename
    with the full path. If the file is not existing yet,
    then it will be downloaded.

    Parameters
    ----------
    name : str
        specifies which type of sample file should be returned
        ['air','rain']
    return_object : bool
        return Data object if True, otherwise the filename is returned
    """

    files = {
            'air': {'name' : 'air.mon.mean.nc', 'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc', 'variable': 'air'},
            'rain': {'name' : 'pr_wtr.eatm.mon.mean.nc', 'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/pr_wtr.eatm.mon.mean.nc', 'variable': 'pr_wtr'}
            }

    if name not in files.keys():
        raise ValueError('Invalid sample file')

    fname = get_example_data_directory() + files[name]['name']

    # download data if not existing yet
    if not os.path.exists(fname):
        tdir = get_example_data_directory()
        url = files[name]['url']
        _download_file(url, tdir)
        if not os.path.exists(fname):
            print fname
            raise ValueError('Download failed!')

    # ... here everything should be fine
    if return_object:
        return Data(fname, files[name]['variable'], read=True)
    else:
        return fname

def _download_file(url, tdir):
    """ download URL to target directory tdir """
    if not os.path.exists(tdir):
        os.makedirs(tdir)

    print('Downloading file ... this might take a few minutes')
    print('URL: ' + url)
    curdir=os.getcwd()
    os.chdir(tdir)
    os.system('wget --ftp-user=anonymous --ftp-password=nix ' + url)
    os.chdir(curdir)
    print('Downloading finished!')



