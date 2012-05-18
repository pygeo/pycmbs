#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

import os

class pyCDO():
    '''
    class to use cdo commands via the shell
    alternative to cdo class provided by MPI-M
    '''
    
    def __init__(self,filename,date1,date2):
        self.filename = filename
        if not os.path.exists(self.filename):
            print 'WARNING: file not existing ', self.filename
        self.date1=date1; self.date2=date2
        self.options = '-f nc'
    
    def seasmean(self,force=False):
        '''
        calculate seasmean
        
        force: True = do calculations; False = do calculations only if data not existing yet
        '''
        
        oname = self.filename[:-3] + '_' + self.date1 + '_' + self.date2 + '_' + 'seasmean' + '.nc'
        cmd2  = self.__seldate(oname)
        cmd1   = 'seasmean'
        cmd = 'cdo ' + self.options + ' ' +  cmd1 + ' -' + cmd2
        self.run(cmd,oname,force)
        
        return oname
        
    def remap(self,method='remapcon',force=False,target_grid='t63grid'):
        
        if method == 'remapcon':
            remap_str = 'remapcon'
        elif method == 'remapnn':
            remap_str = 'remapnn'
        else:
            raise ValueError, 'Unknown remapping method'
            
        oname = self.filename[:-3] + '_'  + remap_str + '.nc'
        cmd = 'cdo ' + self.options + ' ' + remap_str + ',' + target_grid +  ' ' + self.filename + ' ' + oname
        self.run(cmd,oname,force)

        return oname
        
        
    def yearmean(self,force=False):
        oname = self.filename[:-3] + '_'  + 'yearmean' + '.nc'
        cmd = 'cdo ' + self.options + ' ' + 'yearmean' + ' ' + self.filename + ' ' + oname
        self.run(cmd,oname,force)

        return oname

    def yseasmean(self,force=False):
        oname = self.filename[:-3] + '_'  + 'yseasmean' + '.nc'
        cmd = 'cdo ' + self.options + ' ' + 'yseasmean' + ' ' + self.filename + ' ' + oname
        self.run(cmd,oname,force)

        return oname

    def yseasstd(self,force=False):
        oname = self.filename[:-3] + '_'  + 'yseasstd' + '.nc'
        cmd = 'cdo ' + self.options + ' ' + 'yseasstd' + ' ' + self.filename + ' ' + oname
        self.run(cmd,oname,force)

        return oname




    def run(self,cmd,oname,force=False):
        
        #todo: generate output filename such that data is written to temporary directory
        
        f_calc=False
        if os.path.exists(oname):
            if force:
                f_calc = True
        else:
            f_calc = True
        
        if f_calc:
            print cmd; os.system(cmd)
        else:
            print 'INFO - File existing, no calculations will be performed ', oname
        
        
    
    def __seldate(self,oname):
        ''' returns a string for seldate command '''
        return 'seldate,' + str(self.date1) + ',' + str(self.date2) + ' ' + self.filename + ' ' + oname
    
        
        
        



#~ c = pyCDO('/home/m300028/shared/dev/svn/alex/sahel_albedo_jsbach/data/model/tra0074_echam6_BOT_mm_1979-2006_precip.nc','2001-01-01','2001-12-31')
#~ c.seasmean(force=True)
