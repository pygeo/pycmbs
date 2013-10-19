#!/usr/bin/env python
import netCDF4 as nc
from pylab import imshow, show, colorbar
import Nio
import numpy as np
import os,sys


def get_ncdata(ifilename):
    
    nfile = nc.Dataset(ifilename)
    cli = nfile.variables['clisccp']
    
    time = nfile.variables['time']
    lat = nfile.variables['lat']
    lon = nfile.variables['lon']

    return cli,time,lat,lon


def sum2d(arr, axis=0):
    ar1 = arr.sum(axis=axis)
    ar2 = ar1.sum(axis=axis)
    del ar1
    return ar2

def set_attrs(newvar,oldvar):
    for att in oldvar.ncattrs():
        newvar.setncattr(att,oldvar.getncattr(att))

def create_ofile(ofilename,time,lat,lon,cli,clisccp,shape=(96,192)):
    ofile = nc.Dataset(ofilename,'w',format='NETCDF3_CLASSIC')
    ofile.createDimension('time',time.shape[0])
    ofile.createDimension('lat',shape[0])
    ofile.createDimension('lon',shape[1])
    otime = ofile.createVariable('time','f8', ('time',))
    olat  = ofile.createVariable('lat','f8',('lat',))
    olon  = ofile.createVariable('lon','f8',('lon',))
    oclisccp = ofile.createVariable('clisccp','f8',('time','lat','lon'))

    set_attrs(otime,time)
    set_attrs(olat,lat)
    set_attrs(olon,lon)
    set_attrs(oclisccp,cli)
    
    otime[:] = time[:]
    oclisccp[:] = clisccp[:]
    olat[:] = lat[:]
    olon[:] = lon[:]

    ofile.close()


def burst_9_types(filename,ldict,outputdir="."):
    
    cli,time,lat,lon = get_ncdata(filename)

    for key in ldict.keys():
        name, bnds = key,ldict[key]
        tb = bnds[0]
        pb = bnds[1]
        clouds = cli[:,tb[0]:tb[1],pb[0]:pb[1]].copy()
        clsum  = sum2d(clouds,axis=1)

        basename = os.path.basename(filename)
        dirname  = outputdir
        ofilename = dirname + "/" + "%s-%s" % (name,basename)

        create_ofile(ofilename,time,lat,lon,cli,clsum)
        print "*** Created new file: %s" % ofilename


def main():

    plow = [0,2]; pmid = [2,4]; phigh = [4,6]
    tlow = [0,3]; tmid = [3,5]; thigh = [5,7]
    
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


    ifile = sys.argv[1]
    burst_9_types(ifile,ldict,outputdir=sys.argv[2])


if __name__ == "__main__":
    main()






