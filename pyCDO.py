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
    
    def __init__(self,filename,date1,date2,force=False):
        '''
        force: ALWAYS force calculation
        '''
        self.filename = filename
        if not os.path.exists(self.filename):
            print 'WARNING: file not existing ', self.filename
        self.date1=date1; self.date2=date2
        self.options = '-f nc'
        self.force = force
    
    def seasmean(self,force=False):
        '''
        calculate seasmean
        
        force: True = do calculations; False = do calculations only if data not existing yet
        '''
        
        oname = self.filename[:-3] + '_' + self.date1 + '_' + self.date2 + '_' + 'seasmean' + '.nc'
        cmd1   = 'seasmean'
        cmd = 'cdo ' + self.options + ' ' + 'seasmean' + ' '  + self.filename + ' ' + oname
        self.run(cmd,oname,force)
        return oname
        
        
    def selmon(self,months,force=False):
        '''
        select months
        
        input: List with months 
        e.g. [1,2,3,4] for JFMA
        '''
        
        s=''
        if 1 in months:
            s=s+'J'
        if 2 in months:
            s=s+'F'
        if 3 in months:
            s=s+'M'
        if 4 in months:
            s=s+'A'
        if 5 in months:
            s=s+'M'
        if 6 in months:
            s=s+'J'
        if 7 in months:
            s=s+'J'
        if 8 in months:
            s=s+'A'
        if 9 in months:
            s=s+'S'
        if 10 in months:
            s=s+'O'
        if 11 in months:
            s=s+'N'
        if 12 in months:
            s=s+'D'
        
        arg = str(months).replace('[','').replace(']','').replace(' ','')
        
        oname = self.filename[:-3] + '_'  + s + '.nc'
        cmd = 'cdo ' + self.options + ' ' + 'selmon,' + arg + ' ' + self.filename + ' ' + oname
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
            if self.force:
                f_calc = True
        else:
            f_calc = True
        
        if f_calc:
            print cmd; os.system(cmd)
        else:
            print 'INFO - File existing, no calculations will be performed ', oname
        
        
    
    def seldate(self,force=False):
        ''' returns a string for seldate command '''
        oname = self.filename[:-3] + '_' + self.date1 + '_' + self.date2 + '.nc'
        cmd = 'cdo ' + self.options + ' ' + 'seldate,' + str(self.date1) + ',' + str(self.date2) + ' ' + self.filename + ' ' + oname
        self.run(cmd,oname,force)
        return oname
    
        
        
        



#~ c = pyCDO('/home/m300028/shared/dev/svn/alex/sahel_albedo_jsbach/data/model/tra0074_echam6_BOT_mm_1979-2006_precip.nc','2001-01-01','2001-12-31')
#~ c.seasmean(force=True)
