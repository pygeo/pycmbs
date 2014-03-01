# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os

def get_example_directory():
    """ returns directory where this file is located """
    return os.path.dirname(os.path.realpath(__file__))

def get_example_data_directory():
    """ returns directory where the example data should be """
    return get_example_directory() + os.sep + 'example_data'

def get_sample_file(name='air'):
    """
    returns filename of example file including full path
    if the file is not existing yet, then it will be downloaded

    Parameters
    ----------
    name : str
        specifies which type of sample file should be returned
        ['air']
    """

    files = {
            'air': {'name' : 'air.mon.mean.nc', 'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc'},
            'rain': {'name' : 'pr_wtr.eatm.mon.mean.nc', 'url': 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/pr_wtr.eatm.mon.mean.nc'}
            }

    if name not in files.keys():
        raise ValueError('Invalid sample file')

    fname = get_example_data_directory() + os.sep + files[name]['name']

    # download data if not existing yet
    if not os.path.exists(fname):
        tdir = get_example_data_directory()
        url = files[name]['url']
        download_file(url, tdir)
        if not os.path.exists(fname):
            raise ValueError('Download failed!')

    # ... here everything should be fine
    return fname


def download_file(url, tdir):
    if not os.path.exists(tdir):
        os.makedirs(tdir)

    print('Downloading file ... this might take a few minutes')
    curdir=os.getcwd()
    os.chdir(tdir)
    os.system('wget --ftp-user=anonymous --ftp-password=nix ' + url)
    os.chdir(curdir)



