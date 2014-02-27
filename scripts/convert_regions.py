"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Convert region information to a GRID which can be used then for further analysis

Region information is often provided as coordinates, but needs to be rasterized to work with pyCMBS functions.
This class implements this functionality. It allows basically to

1) construct a vector (polygon) dataset from coordinates
2) store this as a shapefile
3) rasterize it and store the information as netCDF file

Please note that additional external dependencies exists for this tool. You need to install
o pyshp [http://code.google.com/p/pyshp/]
"""

import shapefile
import os
import numpy as np

class RegionConverter(object):
    """
    main class to convert coordinates into a shape file
    """
    def __init__(self,shpfile=None):
        """
        @param shpfile: name of output shapefile
        @type shpfile: str
        """
        if shpfile == None:
            raise ValueError, 'Please provide name of output shape file!'
        if shpfile[-4:] == '.shp':
            shpfile = shpfile[0:-4]

        self.shpfile = shpfile
        self.dbf = shapefile.Writer(shapefile.POLYGON)
        self.dbf.field('REGION','C','40') #field for region
        self.dbf.field('RID','N')

        self.regions={}

    def add_region(self,lons,lats,name,id):
        """
        sequence of lon/lat coordinates are converted

        @param lons: list of longitudes (degree)
        @type lons: list
        @param lats: list of latitudes (degree)
        @type lats: list

        @param name: name of region
        @type name: str
        @param id: unique id for region
        @type id: int

        """
        if len(lons) != len(lats):
            raise ValueError, 'lat/lon need to have same length!'
        #w.poly(parts=[[[-105.,60.],[-168.022,60.],[-168.022,72.554],[-105.0,72.554]]])
        self.dbf.poly(parts=[zip(lons,lats)]) #add geometry
        self.dbf.record(name,id)
        self.regions.update({id:{'name':name,'lon':lons,'lat':lats}})

    def save(self):
        """
        save results as shapefile
        """
        #--- delete previous shapefile ---
        if os.path.exists(self.shpfile + '.shp'):
            os.remove(self.shpfile + '.shp')
        if os.path.exists(self.shpfile + '.dbf'):
            os.remove(self.shpfile + '.dbf')
        if os.path.exists(self.shpfile + '.shx'):
            os.remove(self.shpfile + '.shx')
        self.dbf.save(self.shpfile)

    def rasterize(self,filename):
        """
        save results as netCDF file
        """
        pass


########################################################
# convert IPCC extremes regions to shapefile
########################################################

def _get_coordinate(s):
    s = s.replace('(','').replace(')','')
    if 'S' in s:
        s = '-' + s
        s=s.replace('S','')
    if 'W' in s:
        s = '-' + s
        s=s.replace('W','')
    s = s.replace('N','').replace('E','').replace(',','')

    return float(s)



R=RegionConverter(shpfile='ipcc.shp')
f = open('../framework/regions/ipcc_regions.txt','r')
for l in f.readlines():
    if l[0] == '#':
        continue
    s = l.split(' ')
    name = s[0]; id = int(s[1])
    lons = []; lats=[]
    for i in range(2,len(s),2):
        lats.append(_get_coordinate(s[i]))
        lons.append(_get_coordinate(s[i+1]))

    R.add_region(lons,lats,name,id)
R.save()

