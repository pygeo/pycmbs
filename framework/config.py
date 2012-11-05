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

import os

class ConfigFile():
    '''
    class to read pyCMBS configuration file
    '''
    def __init__(self,file):
        '''
        file: name of configuration file
        '''
        self.file = file
        if not os.path.exists(self.file):
            raise ValueError, 'Configuration file not existing: ', self.file
        else:
            self.f=open(file,'r')

        self.read()

    def __read_header(self):
        '''
        read commented lines until next valid line
        '''
        x = '####'
        while x[0] == '#':
            a = self.f.readline().replace('\n','')
            x = a.lstrip()
            if len(x) == 0: #only whitespaces
                x='#'
        return x

    def __check_bool(self,x):
        s=x.split(',')
        if int(s[1]) == 1:
            return True
        else:
            return False

    def __check_var(self,x):
        s=x.split(',')
        if int(s[1]) == 1:
            return s[0]
        else:
            return None

    def __read_options(self):
        #read header of variable plot
        self.options={};
        l=self.__read_header()
        l=l.lstrip()

        if 'BASEMAP' in l.upper():
            self.options.update({'basemap':self.__check_bool(l)})

        l = self.f.readline().replace('\n','')
        if 'REPORT=' in l.upper():
            s = l[7:]
            self.options.update({'report' : s.replace(' ','') })
        else:
            sys.exit('report missing in configuration file!')

        l = self.f.readline().replace('\n','')
        if 'TEMP_DIR=' in l.upper():
            s = l[9:]
            if s[-1] != '/':
                s = s+'/'
            self.options.update({'tempdir' : s.replace(' ','') })
        else:
            raise ValueError, 'Temporary directory not specified!'

        l = self.f.readline().replace('\n','')
        if 'CLEAN_TEMPDIR' in l.upper():
            self.options.update({'cleandir':self.__check_bool(l)})
        else:
            raise ValueError, 'Invalid option for clean_tempdir!'

        #//// create / remove directories
        if not os.path.exists(self.options['tempdir']):
            print 'Creating temporary output directory: ', self.options['tempdir']
            os.makedirs(self.options['tempdir'])
        else:
            if self.options['cleandir']:
                print 'Cleaning output directory: ', self.options['tempdir']
                import glob
                files = glob.glob(self.options['tempdir'] + '*.nc')
                for f in files:
                    os.remove(f)
            else:
                print 'Temporary output directory already existing: ', self.options['tempdir']

        #update global variable for CDO temporary directory (needed for CDO processing)
        os.environ.update({'CDOTEMPDIR' : self.options['tempdir'] } )




    def __read_var_block(self):
        #read header of variable plot
        vars=[];
        l=self.__read_header()
        r=self.__check_var(l)
        if r != None:
            vars.append(r)
        while l[0] != '#':
            l = self.f.readline().replace('\n','')
            l = l.lstrip()

            if (len(l) > 0):
                if l[0] == '#':
                    pass
                else:
                    r=self.__check_var(l)
                    if r != None:
                        vars.append(r)
            else:
                l=' '
        return vars

    def __get_model_details(self,s):
        return s.split(',')

    def __read_date_block(self):
        date1 = self.__read_header()
        date2 = self.f.readline().replace('\n','')
        return date1,date2

    def __read_model_block(self):
        models=[];types=[]; experiments=[]; ddirs=[]
        l=self.__read_header()
        model,ty,experiment,ddir = self.__get_model_details(l)
        models.append(model); types.append(ty)
        experiments.append(experiment); ddirs.append(ddir.replace('\n',''))

        has_eof = False
        while not has_eof:
            try:
                l=self.f.next(); l=l.lstrip()
                if (len(l) > 0) & (l[0] != '#'):
                    model,ty,experiment,ddir = self.__get_model_details(l)
                    models.append(model); types.append(ty)
                    experiments.append(experiment); ddirs.append(ddir.replace('\n',''))
            except:
                has_eof=True

        return models,experiments,types,ddirs

    def read(self):
        '''
        read configuration files in 3 blocks
        '''

        print '********************************************************'
        print '* BEGIN Config file'
        print '********************************************************'

        self.__read_options()
        self.variables  = self.__read_var_block()
        self.start_date,self.stop_date  = self.__read_date_block()
        self.models,self.experiments,self.dtypes,self.dirs = self.__read_model_block()

        for k in self.dtypes:
            if k.upper() not in ['CMIP5','JSBACH_BOT','JSBACH_RAW']:
                print k
                raise ValueError, 'Unknown model type'

        print '********************************************************'
        print '* END Config file'
        print '********************************************************'
        print ''



