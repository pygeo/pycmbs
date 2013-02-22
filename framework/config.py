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
import sys
from utils import get_data_pool_directory

class ConfigFile():
    """
    class to read pyCMBS configuration file
    """
    def __init__(self,file):
        """
        @param file: name of parameter file to parse
        @type file: str
        """
        self.file = file
        if not os.path.exists(self.file):
            raise ValueError, 'Configuration file not existing: ', self.file
        else:
            self.f=open(file,'r')
        self.read()

    def __read_header(self):
        """
        read commented lines until next valid line
        """
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
            return s[0],s[2] #name,interval
        else:
            return None,None

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
                sys.stdout.write('     Temporary output directory already existing: ' + self.options['tempdir'] + '\n')

        #update global variable for CDO temporary directory (needed for CDO processing)
        os.environ.update({'CDOTEMPDIR' : self.options['tempdir'] } )


    def __read_var_block(self):
        #read header of variable plot
        vars=[]; vars_interval={}
        l=self.__read_header()
        r,interval=self.__check_var(l)
        if r != None:
            vars.append(r); vars_interval.update({r:interval})
        while l[0] != '#':
            l = self.f.readline().replace('\n','')
            l = l.lstrip()
            if (len(l) > 0):
                if l[0] == '#':
                    pass
                else:
                    r,interval=self.__check_var(l)
                    if r != None:
                        vars.append(r); vars_interval.update({r:interval})
            else:
                l=' '
        return vars,vars_interval

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
        """
        read configuration files in 3 blocks
        """

        sys.stdout.write("\n *** Reading config file... \n")

        self.__read_options()
        self.variables,self.intervals  = self.__read_var_block()
        self.start_date,self.stop_date  = self.__read_date_block()
        self.models,self.experiments,self.dtypes,self.dirs = self.__read_model_block()

        for k in self.dtypes:
            if k.upper() not in ['CMIP5','JSBACH_BOT','JSBACH_RAW']:
                print k
                raise ValueError, 'Unknown model type'

        sys.stdout.write(" *** Done reading config file. \n")



#-----------------------------------------------------------------------------------------------------------------------

class PlotOptions():
    """
    class for plot options

    @todo: consistency check for Gleckler positions and certain mandatory fields! --> avoid problems in the plotting, but find errors in INI files already when reading!
    """

    def __init__(self):
        self.options = {}


    def read(self,cfg):
        """
        read plot option files and store results in a dictionary

        @param cfg Instance of ConfigFile class which has been already initialized (config file has been read already)
        @type cfg ConfigFile


        @return:
        """

        from ConfigParser import SafeConfigParser
        parser = SafeConfigParser()

        for var in cfg.variables:

            """ The plot options are assumed to be in a file that has the same name as the variable to look be analyzed """
            file = './configuration/' + var + '.ini'
            if os.path.exists(file):
                parser.read(file)
                sys.stdout.write('\n *** Reading configuration for %s: ' % var + "\n")
            else:
                raise ValueError, 'Plot option file not existing: ' + file

            """
            generate now a dictionary for each variable

            in each file there needs to be an OPTIONS section that specifies
            the options that shall be applied to all plots for all observational datasets

            The other sections specify the details for each observational dataset
            """
            dl = {}
            for section_name in parser.sections():
                if section_name.upper() == 'OPTIONS': #add global plotting options
                    o = {}
                    for name, value in parser.items(section_name):
                        o.update({name:value})
                    dl.update({'OPTIONS':o})
                else:
                    o = {} #observation specific dictionary
                    for name, value in parser.items(section_name):
                        o.update({name:value})
                    dl.update({section_name:o})

            #/// update options dictionary for this variable
            self.options.update({var:dl})

        #convert options to bool/numerical values
        self._convert_options()

        #check options consistency
        self._check()

    def _convert_options(self):
        """
        Convert options, which are only strings in the beginning to numerical values or execute functions to set directory
        names appropriately
        @return:
        """

        for v in self.options.keys():
            var = self.options[v]
            for s in var.keys(): #each section
                sec = var[s]
                for k in sec.keys():
                    sec.update({k:self.__convert(sec[k])}) #update current variable with valid value



    def __convert(self,s):
        """
        convert a single string into a valid value. The routine recognizes if the string is numerical, or boolean
        or if a command needs to be executed to generate a valid value

        @param s: string with value of option
        @type s: str
        @return: value of the option
        """

        s.replace('\n','')

        if len(s) == 0:
            return None

        #1) check for numerical value
        try:
            x = float(s)
            return x
        except:
            pass

        #2) check for boolean value
        h = s.lstrip().rstrip()
        if h.upper() == 'TRUE':
            return True
        if h.upper() == 'FALSE':
            return False

        #3) check if some code should be executed (specified by starting and trailing #)
        if (h[0] == '#') and (h[-1] == '#'):
            cmd = h[1:-1]
            exec('res = ' + cmd)
            return res

        #4) check if a list is provided
        if (h[0] == '[') and (h[-1] == ']'): #is list
            return self.__get_list(s)

        #5) in any other case return s
        return s

    def __get_list(self,s):
        #return a list made out from a string '[1,2,3,4,5]'
        l = s.replace('[','')
        l = l.replace(']','')
        l=l.split(',')

        o = []
        for i in l:
            try:
                o.append(float(i)) #try numerical conversion
            except:
                #else: return the string itself
                o.append(i)
        return o

    def _check(self):
        """
        check consistency of options specified
        @return:
        """
        cerr = 0
        o = self.options

        #specify here the options that need to be given!
        locopt = ['obs_file','obs_var','gleckler_position','scale_data'] #variables that need to be specified (MUST!) for each observational dataset
        globopt = ['cticks','map_difference','map_seasons','preprocess','reichler_plot','gleckler_plot','hovmoeller_plot'] #options for each variable type

        for v in o.keys(): #all variables
            d = o[v] #dictionary for a specific variable
            if not 'OPTIONS' in d.keys():
                sys.stdout.write('Error: missing OPTIONS %s' % v)
                cerr += 1

            #check global options
            for k in globopt:
                if not k in d['OPTIONS'].keys():
                    sys.stdout.write('Error: missing global option: %s (%s)' % (k,v)  )
                    cerr += 1

            #check local options
            for odat in d.keys(): #odat is key for a specific observational dataset
                if odat.upper() == 'OPTIONS':
                    continue
                for k in locopt: #k is not the index for a specific obs. record
                    if not k in d[odat].keys():
                        sys.stdout.write('Error: missing local option: %s (%s,%s)' % (k,odat,v)  )
                        cerr +=1

        if cerr > 0:
            raise ValueError, 'There were errors in the initialization of plotting options!'


















