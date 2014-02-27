# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""


class Region(object):
    """
    class to specify a Region in pyCMBS. A region defines either
    a rectangle in the data matrix or it can be defined by
    lat/lon coordinates
    """
    def __init__(self, x1, x2, y1, y2, label, type='index', mask=None):
        """
        constructor of class

        x1 : float
            start position in x-direction (either x index or longitude)
        x2 : float
            stop position in x-direction (either x index or longitude)
        y1 : float
            start position in y-direction (either y index or latitude)
        y2 : float
            stop position in y-direction (either y index or latitude)
        label : str
            label of the region
        type : str
            type of coordiantes (x1,x2,y1,y2): I{index} or I{latlon}
        mask : ndarray
            mask that will be applied in addition to
            the coordinates/indices
        """

        if x2 < x1:
            raise ValueError('Invalid X boundaries for region')
        if y2 < y1:
            raise ValueError('Invalid Y boundaries for region')

        if type == 'index':  # coordinates are considered as inidces
            self.x1 = x1
            self.x2 = x2
            self.y1 = y1
            self.y2 = y2
        elif type == 'latlon':
            self.latmin = y1
            self.latmax = y2
            self.lonmin = x1
            self.lonmax = x2

        self.label = label
        self.type = type
        self.mask = mask

    def _get_label(self):
        return self.label.replace(' ', '')

#-----------------------------------------------------------------------

    def get_corners(self):
        """
        return a list of corner coordinates
        (either indices or lat/lon, dependent on the data)

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
        l = [(self.lonmin, self.latmin), (self.lonmin, self.latmax),
             (self.lonmax, self.latmax), (self.lonmax, self.latmin)]
        return l

#-----------------------------------------------------------------------

    def _get_corners_index(self):
        """
        return a list of corner indices
        """
        l = [(self.x1, self.y1), (self.x1, self.y2),
             (self.x2, self.y2), (self.x2, self.y1)]
        return l

#-----------------------------------------------------------------------

    def get_subset(self, x):
        """
        extract region subset from data array x

        Parameters
        ----------
        x : ndarray
            array where the data is extracted from
        """
        if x.ndim == 3:
            return x[:, self.y1:self.y2, self.x1:self.x2]
        elif x.ndim == 2:
            return x[self.y1:self.y2, self.x1:self.x2]
        else:
            raise ValueError('Invalid data array for subsetting! %s'
                             % np.shape(x))
