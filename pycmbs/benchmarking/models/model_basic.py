# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""

from cdo import Cdo
from pycmbs.data import Data
import tempfile as tempfile
import copy
import glob
import os
import sys
import numpy as np

from pycmbs.benchmarking import preprocessor
from pycmbs.benchmarking.utils import get_T63_landseamask, get_temporary_directory


class Model(Data):
    """
    This class is the main class, specifying a climate model or a particular run
    Sub-classes for particular models or experiments are herited from this class
    """
    def __init__(self, data_dir, dic_variables, name='', intervals=None, **kwargs):
        """
        constructor for Model class

        Parameters
        ----------
        intervals : dict
            a dictionary from configuration, that specifies the temporal interval to be used within each analyis

        INPUT
        -----
        filename: name of the file to read data from (single file currently)
        could be improved later by more abstract class definitions

        dic_variables: dictionary specifiying variable names for a model
        e.g. 'rainfall','var4'
        """

        # check
        if intervals is None:
            raise ValueError('Invalid intervals for Model data: needs specification!')

        # set a list with different datasets for different models
        self.dic_vars = dic_variables
        self.intervals = intervals

        # set some metadata
        self.name = name
        self.data_dir = data_dir

        if 'start_time' in kwargs.keys():
            self.start_time = kwargs['start_time']
        else:
            self.start_time = None
        if 'stop_time' in kwargs.keys():
            self.stop_time = kwargs['stop_time']
        else:
            self.stop_time = None

    def save(self, directory, prefix=None):
        """
        save model variables to file

        Parameters
        ----------
        directory : str
            directory where to store the data to
        prefix : str
            file prefix [obligatory]
        """
        if prefix is None:
            raise ValueError('File prefix needs to be given!')
        if not os.path.exists(directory):
            os.makedirs(directory)
        if directory[-1] != os.sep:
            directory += os.sep

        for k in self.variables.keys():
            if isinstance(self.variables[k], tuple):
                pass
            else:
                if self.variables[k] is not None:
                    self.variables[k].save(directory + prefix + '_' + k.strip().upper() + '.nc', varname=k.strip().lower(), delete=True, mean=False, timmean=False)

    def get_data(self):
        """
        central routine to extract data for all variables
        using functions specified in derived class
        """
        self.variables = {}
        for k in self.dic_vars.keys():
            self._actplot_options = self.plot_options.options[k]['OPTIONS']  # set variable specific options (needed for interpolation when reading the data)

            routine = self.dic_vars[k]  # get name of routine to perform data extraction
            interval = self.intervals[k]

            cmd = 'dat = self.' + routine

            if hasattr(self, routine[0:routine.index('(')]):  # check if routine name is there
                print cmd
                exec(cmd)

                # if a tuple is returned, then it is the data + a tuple for the original global mean field
                if 'tuple' in str(type(dat)):
                    self.variables.update({k : dat[0]})  # update field with data
                    self.variables.update({k + '_org': dat[1]})  # (time, meanfield, originalfield)
                else:
                    self.variables.update({k : dat})  # update field with data
            else:
                print k
                print self.dic_vars
                print routine
                print('WARNING: unknown function to read data (skip!), variable: %s ' % k)
                self.variables.update({k: None})
                assert False


class MedianModel(Model):
    """
    implements a class that allows to handle multiple ensembles and calculate
    ensemble statistics
    """
    def __init__(self, dic_variables, intervals=None, **kwargs):
        super(MeanModel, self).__init__(None, dic_variables,
                                        name='median-model',
                                        intervals=intervals, **kwargs)
        self.n = 0  # specifies only the number of models that were used in general, does NOT specify if the data is valid!
        self._unique_name = 'model_median'
        self.ensmean_called = False

    def add_member(self, M):
        """
        add ensemble member
        this routine adds an ensemble member to the already existing
        dataset; note: to save memory, we
        don't store the individual members, but calculate the median
        on the fly
        """

        if not (M, Model):
            if not issubclass(M, Model):
                raise ValueError('Model instance or derived class expected here!')

        if self.n == 0:
            # copy only variables that are valid. This will be the basis
            # for the whole output
            self.variables = {}
            for k in M.variables.keys():
                if M.variables[k] is not None:
                    self.variables.update({k: M.variables[k]})
            self.name = 'median-model'
            self._unique_name = 'model_median'
            self.n_above = {}
            self.n_below = {}

            # init arrays that specify how much datasets are above/below
            # the current value
            for k in self.variables.keys():
                if self.variables[k] is not None:
                    self.n_above.update({k: np.zeros_like(self.variables[k].data)})
                    self.n_below.update({k: np.zeros_like(self.variables[k].data)})
        else:
            # Processing in case that the data arrays already contain
            # data. The actual value is compared to the new value
            # and the information is stored how many values are
            # above/below
            for k in self.variables.keys():
                # the variables[] list contains Data objects!
                # hlp1 = self.variables[k]  # is a Data object or a tuple! The Data object contains already the climatological mean value!
                if k in M.variables.keys():
                    hlp2 = M.variables[k]
                else:
                    print('Warning: variable %s is not available!' % k)
                    hlp2 = None
                if hlp2 is None:
                    print('WARNING: Provided model is missing data field for variable %s' % k.upper())
                    continue

                if isinstance(self.variables[k], tuple):
                    self.variables.update({k: (None, None, None)})
                else:  # mean model!
                    msk_below = M.variables[k].data <= self.variables[k].data
                    raise ValueError('Implementation of median model has not yet been finished!')
                    # could get_percentiles be used for that in some sense???
                    stop

                    up = self.n_above[k]

                    theD = hlp1.copy()
                    theD.add(hlp2, copy=False)  # SUM: by using masked arrays, the resulting field is automatically only valid, when both datasets contain valid information!
                    theD.label = 'Mean-model'
                    self.variables.update({k: theD})
                    nn = self.N[k]
                    self.N.update({k: nn + 1})
                del hlp1, hlp2


class MeanModel(Model):
    """
    implements a class that allows to handle multiple ensembles and calculate
    ensemble statistics
    """
    def __init__(self, dic_variables, intervals=None, **kwargs):
        super(MeanModel, self).__init__(None, dic_variables,
                                        name='mean-model',
                                        intervals=intervals, **kwargs)
        self.n = 0  # specifies only the number of models that were used in general, does NOT specify if the data is valid!
        self._unique_name = 'model_mean'
        self.model_mean = None
        self.ensmean_called = False

    def add_member(self, M):
        """
        add ensemble member
        this routine adds an ensemble member to the already existing
        dataset; note: to save memory, we
        don't store the individual members, but just the sum of them!

        one needs to call ensmean() afterwards to get the ensemble mean!
        """

        if not (M, Model):
            if not issubclass(M, Model):
                raise ValueError('Model instance or derived class expected here!')

        if self.n == 0:
            tmp = M.copy()
            self.variables = tmp.variables
            del tmp
            self.name = 'mean-model'
            self._unique_name = 'model_mean'
            self.N = {}

            for k in self.variables.keys():  # count for each variable the number of valid models
                if self.variables[k] is not None:
                    self.N.update({k: 1})
                else:
                    self.N.update({k: 0})
        else:
            # Do processing for each variable ...
            for k in self.variables.keys():
                # the variables[] list contains Data objects!
                hlp1 = self.variables[k]  # is a Data object or a tuple! The Data object contains already the climatological mean value!
                if k in M.variables.keys():
                    hlp2 = M.variables[k]
                else:
                    print('Warning: variable %s is not available!' % k)
                    hlp2 = None

                if hlp2 is None:
                    print('WARNING: Provided model is missing data field for variable %s' % k.upper())
                    continue

                if isinstance(hlp1, tuple):
                    self.variables.update({k: (None, None, None)})
                else:  # mean model!
                    if hlp1 is None:
                        continue
                    theD = hlp1.copy()
                    theD.add(hlp2, copy=False)  # SUM: by using masked arrays, the resulting field is automatically only valid, when both datasets contain valid information!
                    theD.label = 'Mean-model'
                    self.variables.update({k: theD})
                    nn = self.N[k]
                    self.N.update({k: nn + 1})
                del hlp1, hlp2

        self.n += 1  # specifies only the number of models that were used in general, does NOT specify if the data is valid!

    def ensmean(self):
        """
        calculate ensemble mean statistics
        """
        if self.ensmean_called:
            raise ValueError('Ensemble mean has been already called! MUST NOT be called a second time !')

        # now we have the sum and can calculate the average
        for k in self.variables.keys():
            if isinstance(self.variables[k], tuple):
                pass
            else:
                if self.variables[k] is not None:
                    self.variables[k].mulc(1. / float(self.N[k]), copy=False)  # weight with number of models
                    self.variables[k].label = 'ensmean_' + k
        self.ensmean_called = True

#------------------------------------------------------------------------------
