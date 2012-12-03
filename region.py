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

class Region():
    """
    class to specify a Region in pyCMBS. A region defines either
    a rectangle in the data matrix or it can be defined by
    lat/lon coordinates
    """
    def __init__(self,x1,x2,y1,y2,label,type='index',mask=None):
        """
        constructor of class

        @param x1: start position in x-direction (either x index or longitude)
        @type x1: float

        @param x2: stop position in x-direction (either x index or longitude)
        @type x2: float

        @param y1: start position in y-direction (either y index or latitude)
        @type y1: float

        @param y2: stop position in y-direction (either y index or latitude)
        @type y2: float

        @param label: label of the region
        @type label: str

        @param type: type of coordiantes (x1,x2,y1,y2): I{index} or I{latlon}
        @type type: str

        @param mask: mask that will be applied in addition to the coordinates/indices
        @type mask: array(:,:)

        """

        if x2 < x1:
            raise ValueError, 'Invalid X boundaries for region' #, x1, x2
        if y2 < y1:
            raise ValueError, 'Invalid Y boundaries for region' #, y1, y2

        if type == 'index': #coordinates are considered as inidces
            self.x1=x1; self.x2=x2
            self.y1=y1; self.y2=y2
        elif type=='latlon':
            self.latmin = y1; self.latmax=y2
            self.lonmin = x1; self.lonmax=x2

        self.label=label
        self.type = type
        self.mask = mask

#-----------------------------------------------------------------------

    def get_corners(self):
        """
        return a list of corner coordinates (either indices or lat/lon, dependent on the data)

        @return list of coordinates
        """
        if self.type == 'latlon':
            return self._get_corners_latlon()
        else:
            return self._get_corners_index()

#-----------------------------------------------------------------------

    def _get_corners_latlon(self):
        """
        return a list of corner lat/lon
        """
        l = []
        l.append( (self.lonmin,self.latmin)   )
        l.append( (self.lonmin,self.latmax)   )
        l.append( (self.lonmax,self.latmax)   )
        l.append( (self.lonmax,self.latmin)   )
        return l

#-----------------------------------------------------------------------

    def _get_corners_index(self):
        """
        return a list of corner indices
        """
        l = []
        l.append( (self.x1,self.y1)   )
        l.append( (self.x1,self.y2)   )
        l.append( (self.x2,self.y2)   )
        l.append( (self.x2,self.y1)   )
        return l

#-----------------------------------------------------------------------

    def get_subset(self,x):
        """
        extract region subset from data array x

        @param x: array where the data is extracted from
        @type x: array
        """
        if x.ndim == 3:
            return x[:,self.y1:self.y2,self.x1:self.x2]
        elif x.ndim == 2:
            return x[self.y1:self.y2,self.x1:self.x2]
        else:
            raise ValueError, 'Invalid data array for subsetting!', np.shape(x)
