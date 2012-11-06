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

class ExternalAnalysis():
    def __init__(self,executable,template,tags,output_directory='./',options=''):
        '''
        constructor for ExternalAnalysis class

        @param executable: name of system executable to be called to run analysis. This string need to include all options and have a tag <INPUTFILE>
                           to specify the filename to be called e.g. 'matlab -r <INPUTFILE>'
        @type executable: str

        @param template: name of template file to be modified
        @type template: str

        @param tags: tag that specifies which variables should be replaced with which value.
                     The dictionary must have the list consists of dictonary entries like the following example
                     d = [{'tag':'myfile','value':'myvaluestring'},{'tag':'myfile1','value':'myvaluestring1'}]
        @type tags: list

        @param output_directory: output directory to write output of external analysis to
        @type output_directory: str

        @param options: options that will be appended to the filename in the end
        @type options: str

        '''

        if not '<INPUTFILE>' in executable:
            raise ValueError, 'The commmand needs to include <INPUTFILE> tag'


        self.exe = executable
        self.tags = tags
        self.output_directory = output_directory
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
        self.template = template
        self.options = options


    def _create_script(self):
        '''
        replace tags in a template file

        returns filename
        '''

        #--- check if template file existing
        if not os.path.exists(self.template):
            raise ValueError, '*** Template file not existing! ' + self.template

        #--- copy template file
        filename = self.output_directory + os.path.basename(self.template) #target file
        if os.path.exists(filename):
            os.remove(filename)

        #--- replace template content
        f = open(self.template,'r')
        s = []
        for l in f.readlines():
            d = l
            for tag in self.tags: #replace all tags
                d = d.replace(tag['tag'],tag['value'])
            s.append(d)
        f.close()

        #--- write output
        of = open(filename,'w')
        of.writelines(s); of.close()

        return filename


    def run(self,execute=True,remove_extension=False):
        '''
        run external program
        the execution of the program will be done in the directory where the modified script is located
        this is done, because programs like e.g. matlab have problems when the scriptname is provided with
        a pathname

        @param execute: execute command = run in shell; if false, then the command is simply printed (good for debugging)
        @type execute: bool

        @param remove_extension: remove extension from filename. This is for instance needed for matlab calls from the command line, as matlab does NOT accept script names with a fileextension
        @type remove_extension: bool
        '''

        #/// generate run script (returns filename) ///
        filename = self._create_script()
        if remove_extension:
            filename = os.path.splitext(filename)[0]
            
        #/// split filename into path and 
        thedir = os.path.dirname(filename)
        thefile = os.path.basename(filename)
        curdir = os.getcwd()

        #/// run script
        cmd = self.exe.replace('<INPUTFILE>',thefile)+self.options
        if execute:
            os.chdir(thedir)  #change to directory where script is located
            r = os.system(cmd) #@todo: use subprocess
            os.chdir(curdir) #go back
        else:
            print cmd




#~ template='test.txt'
#~ tags = [{'tag':'ALEX','value':'peakXXXmat'}]
#~
#~ E=ExternalAnalysis('matlab -r <INPUTFILE>',template,tags,output_directory='./tmp/')
#~ E.run()

