# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
For COPYRIGHT, LICENSE and AUTHORSHIP please referr to
the pyCMBS licensing details.
"""
import os
import glob


class EnsemblePreprocessor(object):
    """
    class to perform preprocessing for model ensemble data.
    The purpose is in particular to preprocess ensemble members
    which are stored in different files and provide functionality for
    ensemble statistic calculations

    The processor basically provides a logic. Number crunching is
    exclusively done using the CDO's
    """

    def __init__(self, data_dir, outfile):
        """
        Parameters
        ----------
        outfile : str
            name of the output file to be generated
        """
        self.output_dir, self.outfile = os.path.split(outfile)
        if len(self.output_dir) == 0:
            self.output_dir = './'
        self.data_dir = data_dir
        if self.data_dir[-1] != os.sep:
            self.data_dir += os.sep
        if self.output_dir[-1] != os.sep:
            self.output_dir += os.sep
        self.mergetime_files = []

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)


class CMIP5Preprocessor(EnsemblePreprocessor):
    def __init__(self, data_dir, outfile, variable, model, experiment, mip='Amon', realm='atmos', institute=None):
        super(CMIP5Preprocessor, self).__init__(data_dir, outfile)
        if institute is None:
            raise ValueError('An institute name needs to be provided!')
        self.institute = institute
        self.variable = variable
        self.model = model
        self.experiment = experiment
        self.mip = mip
        self.realm = realm

    def _filelist(self, l):
        r = ''
        for x in l:
            r += x + ' '
        return r

    def _log(self, s):
        logfile = self.output_dir + 'cmip5_preprocessor.log'
        o = open(logfile, 'a')
        o.write(s + '\n')
        o.close()

    def mergetime(self, n, start_time=None, stop_time=None, delete=False):
        """
        perform mergetime for a single ensemble member
        In case that only a single file is used, generate only a link
        to the original data

        Parameters
        ----------
        n : int
            number of ensemble member
        start_time : datetime
            start time of period
        stop_time : datetime
            end time of period
        """
        if n not in self.ensemble_files.keys():
            print n
            raise ValueError('Ensemble not existing!')

        if start_time is not None:
            if stop_time is None:
                raise ValueError('Both, start_time and stop_time need to be specified!')
        if stop_time is not None:
            if start_time is None:
                raise ValueError('Both, start_time and stop_time need to be specified!')
        if start_time is not None:
            selstr = 'seldate,' + str(start_time)[0:10] + ',' + str(stop_time)[0:10]
        else:
            selstr = ''

        ens_files = self.ensemble_files[n]
        self._log(self.model + '\t' + str(n) + '\t' + str(len(ens_files)))  # model name, ensemble nr., files per ensemble member

        # create string with filenames
        fstr = self._filelist(ens_files)

        # output file
        ofile = self.output_dir + os.path.basename(ens_files[0]).split('r'+str(n)+'i1p1')[0] + 'r' + str(n) + 'i1p1' + '_mergetime.nc'

        if start_time is not None:
            ofile += '_' + str(start_time)[0:10] + '_' + str(stop_time)[0:10]
        ofile += '.nc'

        self.mergetime_files.append(ofile)

        # cdo

        if selstr == '':
            cmd = 'cdo -f nc mergetime ' + fstr + ' ' + ofile
        else:
            # merge first and the select time (not the most performant way to do it)
            tmpfile = ofile + '.tmp.nc'
            if os.path.exists(tmpfile):
                os.remove(tmpfile)
            cmd1 = 'cdo -f nc mergetime ' + fstr + ' ' + tmpfile
            os.system(cmd1)

            if not os.path.exists(tmpfile):
                tmp_s = 'Error in creating temporary file: ABORT! ' + self.model + ' ' + self.experiment + ' ' + tmpfile
                print tmp_s
                self._log(tmp_s)

            cmd = 'cdo -f nc ' + selstr + ' ' + tmpfile + ' ' + ofile

        # calculate final output file
        if os.path.exists(ofile):
            if delete:
                os.remove(ofile)
            else:
                self.mergetime_files.append(ofile)
                print('File already existing ... no processing is done')
                return

        print('Doing temporal preprocessing for ensemble member ... %s' % n)
        os.system(cmd)
        if os.path.exists(ofile):  # just ensure that everything wen well
            self.mergetime_files.append(ofile)

        if os.path.exists(tmpfile):
            os.remove(tmpfile)

    def get_ensemble_files(self, maxens=50):
        """
        create a dictionary with filenames for the different ensemble
        members

        Parameters
        ----------
        maxens : int
            maximum ensemble member size
        """

        # file wildcard
        #nfiles = len(glob.glob(w + '*.nc'))

        # create filelist for each ensemble
        res = {}
        cnt = 0
        for i in xrange(1, maxens+1):
            w1 = self._get_file_wildcard(i)
            files = glob.glob(w1)
            cnt += len(files)
            if len(files) > 0:
                res.update({i: files})

        # check that all files were used
        #if nfiles != cnt:
        #    print nfiles, cnt
        #    raise ValueError('ERROR: not all files were processed!')
        self.ensemble_files = res

    def _get_file_wildcard(self, ens):
        #return self.data_dir + self.variable + '_Amon_' + self.model + '_' + self.experiment + '_r'
        p = self.data_dir + self.institute + os.sep + self.model + os.sep + self.experiment + os.sep + 'mon' + os.sep + self.realm + os.sep + self.mip + os.sep + 'r' + str(ens) + 'i1p1' + os.sep + self.variable + os.sep + self.variable + '_' + self.mip + '_' + self.model + '_' + self.experiment + '_r' + '*.nc'
        return p

    def mergetime_ensembles(self, delete=False, start_time=None, stop_time=None):
        self.get_ensemble_files()
        for i in self.ensemble_files.keys():
            self.mergetime(i, delete=delete, start_time=start_time, stop_time=stop_time)

    def ensemble_mean(self, delete=False, start_time=None, stop_time=None):
        """
        Parameters
        ----------
        start_time : datetime
            start time of period
        stop_time : datetime
            end time of period
        """
        print 'Doing ensemble mean calculation ...'
        # if temporal merged files are not yet there, do preprocessing
        if len(self.mergetime_files) < 2:
            self.mergetime_ensembles(delete=delete, start_time=start_time, stop_time=stop_time)
        if len(self.mergetime_files) < 2:
            print self.mergetime_files
            print 'No ensemble mean calculation possible as not enough files!'
            self._log('No ensemble mean calculation possible as not enough files! ' + self.institute + ' ' + self.model + ' ' + self.experiment)
        fstr = self._filelist(self.mergetime_files)

        # ensemble mean calculation
        ofile = self.output_dir + self.outfile
        ofile = os.path.splitext(ofile)[0]
        if start_time is not None:
            ofile += '_' + str(start_time)[0:10] + '_' + str(stop_time)[0:10]
        ofile += '.nc'
        self.outfile = ofile
        cmd = 'cdo -f nc ensmean ' + fstr + ' ' + ofile
        if os.path.exists(ofile):
            if delete:
                os.remove(ofile)
                os.system(cmd)
            else:
                print('File already existing ... no processing is done')
        else:
            os.system(cmd)

        # ensemble standard deviation
        ofilestd = self.outfile.replace('_ensmean', '_ensstd')
        cmd = 'cdo -f nc ensstd ' + fstr + ' ' + ofilestd
        if os.path.exists(ofilestd):
            if delete:
                os.remove(ofilestd)
                os.system(cmd)
            else:
                print('File already existing ... no processing is done')
        else:
            os.system(cmd)

        return ofile


class CMIP5ModelParser(object):
    """
    a parser to retrieve model and institute names
    """

    def __init__(self, root_dir):
        self.root_dir = root_dir
        if not os.path.exists(self.root_dir):
            raise ValueError('Path not existing!')

        self.institutes = None
        self.models = None

    def get_institutes(self):
        return [os.path.basename(f) for f in glob.glob(self.root_dir + '*')]

    def get_all_models(self):
        """ get a list with all models and institutes """
        r = {}
        for i in self.get_institutes():
            r.update({i: self._get_models4institute(i)})
        return r

    def _get_models4institute(self, institute):
        return [os.path.basename(f) for f in glob.glob(self.root_dir + institute + os.sep + '*')]

import datetime as dt
# from pyCMBS import *

output_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation/dummy_out/amip_rsus/'

data_dir = '/data/share/mpiles/TRS/PROJECT_RESULTS/EvaClimod/CMIP5_RAWDATA_NEW/radiation/'

the_experiment = 'amip'
the_variable = 'rsds'


# init parser that returns a list of institutes and models
CP = CMIP5ModelParser(data_dir)
model_list = CP.get_all_models()

for institute in model_list.keys():
    for model in model_list[institute]:
        output_file = output_dir + the_variable + '_Amon_' + model + '_' + the_experiment + '_ensmean.nc'
        E = CMIP5Preprocessor(data_dir, output_file, the_variable, model, the_experiment, institute=institute)
        E.get_ensemble_files()
        E.ensemble_mean(delete=False, start_time=dt.datetime(1979,1,1), stop_time=dt.datetime(2012,12,31))

#d=Data(E.output_dir + E.outfile, 'rsds', read=True)


#~ E.mergetime(2, delete=True, start_time=dt.datetime(1999,3,10), stop_time=dt.datetime(2003,5,16))
#~ E.mergetime(2, delete=True)

#~ tmp='/home/m300028/shared/temp/ensmean/rsds/amip/raw/testout/rsds_Amon_GFDL-HIRAM-C180_amip_r2_mergetime.nc'
#~ d = Data(tmp, 'rsds', read=True)
