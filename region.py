#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

class Region():
    '''
    class to specify a Region in pyCMBS
    '''
    def __init__(self,x1,x2,y1,y2,label,type='index',mask=None):
        
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
    
    def get_corners(self):
        '''
        return a list of corner coordinates
        '''
        if self.type == 'latlon':
            return self._get_corners_latlon()
        else:
            return self._get_corners_index()
    
    def _get_corners_latlon(self):
        l = []
        l.append( (self.lonmin,self.latmin)   )
        l.append( (self.lonmin,self.latmax)   )
        l.append( (self.lonmax,self.latmax)   )
        l.append( (self.lonmax,self.latmin)   )
        return l
        
    def _get_corners_index(self):
        sys.exit('Routine for corner indices not implemented yet!')
        
        
        
        
    def get_subset(self,x):
        '''
        extract region subset from data array x
        '''
        if x.ndim == 3:
            return x[:,self.y1:self.y2,self.x1:self.x2]
        elif x.ndim == 2:
            return x[self.y1:self.y2,self.x1:self.x2]
        else:
            raise ValueError, 'Invalid data array for subsetting!', np.shape(x)
