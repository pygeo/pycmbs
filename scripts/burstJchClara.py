#!/usr/bin/env python
import netCDF4 as nc
from pylab import imshow, show, colorbar
import Nio
import numpy as np
import os,sys

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""


def rebin_sum(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).sum(-1).sum(1)

def rebin_mean(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def get_ncdata(jch_filename,cfc_filename):
    
    nfile = nc.Dataset(jch_filename)

    jliq = nfile.variables['jch_liq']
    jice = nfile.variables['jch_ice']
    jall = jice[0,:,:,:,:] + jliq[0,:,:,:,:]

    time = nfile.variables['time']
    lat = nfile.variables['lat']
    lon = nfile.variables['lon']

    cfile = nc.Dataset(cfc_filename)
    cnobs = cfile.variables['cfc_day_nobs'][0,:,:]
    cnobs_rs = rebin_sum(cnobs,shape=(jall.shape[-2],jall.shape[-1]))

    
    return jall,cnobs_rs,time,lat,lon


def sum2d(arr, axis=0):
    ar1 = arr.sum(axis=axis)
    ar2 = ar1.sum(axis=axis)
    del ar1
    return ar2

def set_attrs(newvar,oldvar):
    for att in oldvar.ncattrs():
        newvar.setncattr(att,str(oldvar.getncattr(att)))

def create_ofile(ofilename,time,lat,lon,cli,clisccp,shape=(96,192)):
    ofile = nc.Dataset(ofilename,'w',format='NETCDF3_CLASSIC')
    ofile.createDimension('time',time.shape[0])
    ofile.createDimension('lat',shape[0])
    ofile.createDimension('lon',shape[1])
    otime = ofile.createVariable('time','f8', ('time',))
    olat  = ofile.createVariable('lat','f8',('lat',))
    olon  = ofile.createVariable('lon','f8',('lon',))
    oclisccp = ofile.createVariable('jch','f8',('time','lat','lon'))

    set_attrs(otime,time)
    set_attrs(olat,lat)
    set_attrs(olon,lon)
    #set_attrs(oclisccp,cli)

    print clisccp.shape

    otime[:] = time[:]
    oclisccp[:] = clisccp
    olat[:] = lat[:]
    olon[:] = lon[:]

    ofile.close()


def burst_9_types(jchFilename,cfcFilename,ldict,outputdir="."):
    
    cli,cnobs,time,lat,lon = get_ncdata(jchFilename,cfcFilename)
    print cli.sum(),cnobs.sum()

    for key in ldict.keys():
        name, bnds = key,ldict[key]
        tb = bnds[0]
        pb = bnds[1]
        clouds = cli[tb[0]:tb[1],pb[0]:pb[1],:,:].copy()
        clsum  = sum2d(clouds,axis=0).astype(np.float32) / cnobs * 100
        print clsum
        print clouds.sum(),clsum.sum()

        basename = os.path.basename(jchFilename)
        dirname  = outputdir
        ofilename = dirname + "/" + "%s-%s" % (name,basename)

        create_ofile(ofilename,time,lat,lon,cli,clsum,shape=(lat.shape[0],lon.shape[0]))
        print "*** Created new file: %s" % ofilename


def main():

    tlow = [0,5]; tmid = [5,9]; thigh = [9,13]
    plow = [10,14]; pmid = [6,10]; phigh = [0,6]
    
    ldict = {
        'ci': [ tlow ,  phigh ],
        'cs': [ tmid ,  phigh ],
        'cb': [ thigh,  phigh ],
        'ac': [ tlow ,  pmid  ],
        'as': [ tmid ,  pmid  ],
        'ns': [ thigh,  pmid  ],
        'cu': [ tlow ,  plow  ],
        'sc': [ tmid ,  plow  ],
        'st': [ thigh,  plow  ] }


    jchFile = sys.argv[1]
    cfcFile = sys.argv[2]
    
    burst_9_types(jchFile,cfcFile,ldict,outputdir=sys.argv[3])


if __name__ == "__main__":
    main()
