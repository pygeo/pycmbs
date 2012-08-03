#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

    def __check_var(self,x):
        s=x.split(',')
        if int(s[1]) == 1:
            return s[0]
        else:
            return None

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
        self.variables  = self.__read_var_block()
        self.start_date,self.stop_date  = self.__read_date_block()
        self.models,self.experiments,self.dtypes,self.dirs = self.__read_model_block()

        for k in self.dtypes:
            if k.upper() not in ['CMIP5','JSBACH']:
                print k
                raise ValueError, 'Unknown model type'
