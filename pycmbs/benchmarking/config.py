# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import os
import sys
import numpy as np
import pylab as pl
import glob
from cdo import Cdo

from pycmbs.data import Data
from pycmbs.region import Region, RegionParser
from pycmbs.polygon import Raster
from pycmbs.polygon import Polygon as pycmbsPolygon

from pycmbs.benchmarking.utils import get_data_pool_directory
from pycmbs.benchmarking.utils import get_generic_landseamask, get_T63_landseamask


class ConfigFile(object):
    """
    class to read pyCMBS configuration file
    """
    def __init__(self, file):
        """
        file : str
            name of parameter file to parse
        """
        self.file = file
        if not os.path.exists(self.file):
            raise ValueError('Configuration file not \
                              existing: % ' % self.file)
        else:
            self.f = open(file, 'r')
        self.read()

    def __read_header(self):
        """
        read commented lines until next valid line
        """
        x = '####'
        while x[0] == '#':
            a = self.f.readline().replace('\n', '')
            x = a.lstrip()
            if len(x) == 0:  # only whitespaces
                x = '#'
        return x

    def __check_bool(self, x):
        s = x.split(',')
        if int(s[1]) == 1:
            return True
        else:
            return False

    def __check_var(self, x):
        s = x.split(',')
        if int(s[1]) == 1:
            return s[0], s[2]  # name,interval
        else:
            return None, None

    def __read_options(self):
        # read header of variable plot
        self.options = {}
        l = self.__read_header()
        l = l.lstrip()

        if 'BASEMAP' in l.upper():
            self.options.update({'basemap': self.__check_bool(l)})

        l = self.f.readline().replace('\n', '')
        if 'REPORT=' in l.upper():
            s = l[7:]
            self.options.update({'report': s.replace(' ', '')})
        else:
            raise ValueError('Report missing in configuration file!')

        l = self.f.readline().replace('\n', '')
        if 'REPORT_FORMAT=' in l.upper():
            s = l[14:].strip()
            s = s.lower()
            if s not in ['png', 'pdf']:
                raise ValueError('Invlid option for report format [png,pdf]: %s' % s)
            else:
                self.options.update({'report_format': s})
        else:
            raise ValueError('report format missing in configuration file!')

        l = self.f.readline().replace('\n', '')
        if 'AUTHOR=' in l.upper():
            s = l[7:]
            self.options.update({'author': s})
        else:
            raise ValueError('author missing in configuration file!')

        l = self.f.readline().replace('\n', '')
        if 'TEMP_DIR=' in l.upper():
            s = l[9:]
            if s[-1] != os.sep:
                s = s + os.sep
            self.options.update({'tempdir': s.replace(' ', '')})
        else:
            raise ValueError('Temporary directory not specified!')

        l = self.f.readline().replace('\n', '')
        if 'CLEAN_TEMPDIR' in l.upper():
            self.options.update({'cleandir': self.__check_bool(l)})
        else:
            raise ValueError('Invalid option for clean_tempdir!')

        l = self.f.readline().replace('\n', '')
        if 'SUMMARY_ONLY' in l.upper():
            self.options.update({'summary': self.__check_bool(l)})
        else:
            raise ValueError('Invalid option for SUMMARY_ONLY!')

        l = self.f.readline().replace('\n', '')
        if 'CONFIG_DIR=' in l.upper():
            s = l[11:]
            if s[-1] != os.sep:
                s = s + os.sep
            if not os.path.exists(s):
                raise ValueError('Configuration path is invalid: %s' % s)
            self.options.update({'configdir': s.replace(' ', '')})
        else:
            raise ValueError('CONFIG directory not specified!')

        l = self.f.readline().replace('\n', '')
        if 'OUTPUT_DIRECTORY=' in l.upper():
            s = l[17:]
            if s[-1] != os.sep:
                s = s + os.sep  # this is the root output directory
            if not os.path.exists(s):
                os.makedirs(s)
            self.options.update({'outputdir': s.replace(' ', '')})
        else:
            raise ValueError('OUTPUT directory not specified!')

        # create / remove directories
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

        # update global variable for CDO temporary directory (needed for CDO processing)
        os.environ.update({'CDOTEMPDIR': self.options['tempdir']})

    def __read_var_block(self):
        # read header of variable plot
        vars = []
        vars_interval = {}
        l = self.__read_header()
        r, interval = self.__check_var(l)
        if r is not None:
            vars.append(r)
            vars_interval.update({r: interval})
        while l[0] != '#':
            l = self.f.readline().replace('\n', '')
            l = l.lstrip()
            if len(l) > 0:
                if l[0] == '#':
                    pass
                else:
                    r, interval = self.__check_var(l)
                    if r is not None:
                        vars.append(r)
                        vars_interval.update({r: interval})
            else:
                l = ' '

        return vars, vars_interval

    def __get_model_details(self, s):
        return s.split(',')

    def __read_date_block(self):
        date1 = self.__read_header()
        date2 = self.f.readline().replace('\n', '')
        tmp = self.f.readline().replace('\n', '')  # same time for observations
        same_for_obs = self.__check_bool(tmp)
        return date1, date2, same_for_obs

    def __read_model_block(self):
        models = []
        types = []
        experiments = []
        ddirs = []
        l = self.__read_header()
        model, ty, experiment, ddir = self.__get_model_details(l)
        ddir = ddir.rstrip()
        if ddir[-1] != os.sep:
            ddir = ddir + os.sep
        models.append(model)
        types.append(ty)
        experiments.append(experiment)
        ddirs.append(ddir.replace('\n', ''))

        has_eof = False
        while not has_eof:
            try:
                l = self.f.next()
                l = l.lstrip()
                if (len(l) > 0) & (l[0] != '#'):
                    model, ty, experiment, ddir = self.__get_model_details(l)
                    ddir = ddir.rstrip()
                    if ddir[-1] != os.sep:
                        ddir = ddir + os.sep
                    models.append(model)
                    types.append(ty)
                    experiments.append(experiment)
                    ddirs.append(ddir.replace('\n', ''))
            except:
                has_eof = True

        return models, experiments, types, ddirs

    def read(self):
        """
        read configuration files in 3 blocks
        """
        sys.stdout.write("\n *** Reading config file... \n")

        self.__read_options()
        self.variables, self.intervals = self.__read_var_block()
        self.start_date, self.stop_date, self.same_time4obs = self.__read_date_block()
        self.models, self.experiments, self.dtypes, self.dirs = self.__read_model_block()

        for k in self.dtypes:
            if k.upper() not in ['CMIP5', 'JSBACH_BOT', 'JSBACH_RAW',
                                 'CMIP3', 'JSBACH_RAW2', 'CMIP5RAW', 'CMIP5RAWSINGLE', 'JSBACH_SPECIAL']:
                raise ValueError('Unknown model type: %s' % k)

        # ensure that up/down fluxes are analyzed in case of albedo
        # in that case the same time sampling is used!
        if 'albedo' in self.variables:
            if 'sis' not in self.variables:
                self.variables.append('sis')
            self.intervals.update({'sis': self.intervals['albedo']})

            if 'surface_upward_flux' not in self.variables:
                self.variables.append('surface_upward_flux')
            self.intervals.update({'surface_upward_flux': self.intervals['albedo']})

        sys.stdout.write(" *** Done reading config file. \n")

    def get_analysis_scripts(self):
        """
        returns names of analysis scripts for all variables as a dictionary
        in general.

        The names of the different analysis routines are taken from file
        <configuration_dir>/analysis_scripts.json

        which has the format

        VARIABLE,NAME OF ANALYSIS ROUTINE
        """
        import json
        jsonfile = self.options['configdir'] + 'analysis_routines.json'
        if not os.path.exists(jsonfile):
            raise ValueError('REQUIRED file analysis_routines.json not existing!')
        d = json.load(open(jsonfile, 'r'))
        return d

    def get_methods4variables(self, variables):
        """
        for a given list of variables, return a dictionary
        with information on methods how to read the data

        The actual information is coming from a json file

        IMPORTANT: all options provided to the routines need to be
        specified here and arguments must be set in calling
        routine get_data()

        Parameters
        ----------
        variables : list
            list of variables to be analzed
        """

        jsonfile = self.options['configdir'] + 'model_data_routines.json'

        import json
        if not os.path.exists(jsonfile):
            raise ValueError('File model_data_analysis.json MISSING!')
        hlp = json.load(open(jsonfile, 'r'))

        res = {}
        for k in hlp.keys():
            # only use the variables that should be analyzed!
            if k in variables:
                res.update({k: hlp[k]})

        # ensure that for albedo processing also the routines
        # for upward and downward shortwave flux are known
        if 'albedo' in variables:
            for k in ['surface_upward_flux', 'sis']:
                if k in hlp.keys():
                    if k not in res.keys():
                        res.update({k: hlp[k]})
                else:
                    err_msg = 'For albedo processing also the ' + k.upper() + ' routines need to be spectified!'
                    raise ValueError(err_msg)

        # implement here also dependencies between variables for analysis
        # e.g. phenology needs faPAR and snow cover fraction. Ensure here that
        # snow cover is also read, even if only phenology option is set
        if ('phenology_faPAR' in variables) and not ('snow' in variables):
            res.update({'snow': hlp['snow']})

        return res


class PlotOptions(object):
    """
    Class for plot options
    """

    def __init__(self):
        self.options = {}

    def _import_regional_file(self, region_file, varname, targetgrid=None, logfile=None):
        """
        check if the regional file can be either imported or if
        regions are provided as vector data. In the latter case
        the regions are rasterized and results are stored in a netCDF
        file

        Parameters
        ----------
        region_file : str
            name of file defining the region. This is either a netCDF
            file which contains the mask as different integer values
            or it is a *.reg file which contains the regions as
            vector data.
        varname : str
            name of variable in netCDF file
        targetgrid : str
            name of targetgrid; either 't63grid' or the name of a file
            with a valid geometry

        Returns
        -------
            region_filename, region_file_varname
        """

        if not os.path.exists(region_file):
            raise ValueError('ERROR: region file is not existing: ' + region_file)

        ext = os.path.splitext(region_file)[1]
        if ext == '.nc':
            # netCDF file was given. Try to read variable
            if varname is None:
                raise ValueError('ERROR: no variable name given!')
            try:
                tmp = Data(region_file, varname, read=True)
            except:
                raise ValueError('ERROR: the regional masking file can not be read!')
            del tmp

            # everything is fine
            return region_file, varname

        elif ext == '.reg':
            # regions were given as vector files. Read it and
            # rasterize the data and store results in a temporary
            # file
            import tempfile

            if targetgrid is None:
                raise ValueError('ERROR: targetgrid needs to be specified for vectorization of regions!')

            if targetgrid == 't63grid':
                ls_mask = get_T63_landseamask(True, area='global', mask_antarctica=False)
            else:
                ls_mask = get_generic_landseamask(True, area='global', target_grid=targetgrid,
                                                  mask_antarctica=False)

            # temporary netCDF filename
            region_file1 = tempfile.mktemp(prefix='region_mask_', suffix='.nc')
            R = RegionParser(region_file)  # read region vector data
            M = Raster(ls_mask.lon, ls_mask.lat)
            polylist = []
            if logfile is not None:
                logf = open(logfile, 'w')
            else:
                logf = None

            id = 1
            for k in R.regions.keys():
                reg = R.regions[k]
                polylist.append(pycmbsPolygon(id, zip(reg.lon, reg.lat)))
                if logf is not None:  # store mapping table
                    logf.write(k + '\t' + str(id) + '\n')
                id += 1

            M.rasterize_polygons(polylist)
            if logf is not None:
                logf.close()

            # generate dummy output file
            O = Data(None, None)
            O.data = M.mask
            O.lat = ls_mask.lat
            O.lon = ls_mask.lon
            varname = 'regions'
            O.save(region_file1, varname=varname, format='nc', delete=True)
            print('Regionfile was store in file: %s' % region_file1)

            # check again that file is readable
            try:
                tmp = Data(region_file1, varname, read=True)
            except:
                print region_file1, varname
                raise ValueError('ERROR: the generated region file is not readable!')
            del tmp

            return region_file1, varname

        else:
            raise ValueError('ERROR: unsupported file type')

    def read(self, cfg):
        """
        read plot option files and store results in a dictionary

        Parameters
        ----------
        cfg : ConfigFile instance
            cfg Instance of ConfigFile class which has been already
            initialized (config file has been read already)
        """
        from ConfigParser import SafeConfigParser

        thevariables = cfg.variables
        # ensure that up/down information also given, when albedo is used
        #~ if 'albedo' in thevariables:
            #~ if 'surface_upward_flux' not in thevariables:
                #~ thevariables.append('surface_upward_flux')
            #~ if 'sis' not in thevariables:
                #~ thevariables.append('sis')

        for var in thevariables:
            parser = SafeConfigParser()

            # The plot options are assumed to be in a file that has the same name as the variable to look be analyzed
            file = cfg.options['configdir'] + var + '.ini'
            if os.path.exists(file):
                sys.stdout.write('\n *** Reading configuration for %s: ' % var + "\n")
                parser.read(file)
            else:
                raise ValueError('Plot option file not existing: %s' % file)

            """
            generate now a dictionary for each variable

            in each file there needs to be an OPTIONS section that specifies
            the options that shall be applied to all plots for all observational datasets

            The other sections specify the details for each observational dataset
            """

            # print('\n*** VARIABLE: %s ***' % var)
            dl = {}
            # print parser.sections()
            for section_name in parser.sections():
                # add global plotting options
                if section_name.upper() == 'OPTIONS':
                    o = {}
                    for name, value in parser.items(section_name):
                        o.update({name: value})
                    dl.update({'OPTIONS': o})
                else:  # observation specific dictionary
                    o = {}
                    for name, value in parser.items(section_name):
                        o.update({name: value})
                    dl.update({section_name: o})

            # update options dictionary for this variable
            self.options.update({var: dl})

            # destroy parser (important, as otherwise problems)
            del parser

        # convert options to bool/numerical values
        self._convert_options()

        # check options consistency
        self._check()

        #reset plotting options if the report should only produce a summary
        if cfg.options['summary']:
            # set the plot options to FALSE which are *not* relevant
            # for summary report
            false_vars = ['map_difference', 'map_seasons',
                          'reichler_plot', 'hovmoeller_plot',
                          'regional_analysis']
            for var in thevariables:
                lopt = self.options[var]
                for vv in false_vars:
                    if vv in lopt['OPTIONS'].keys():
                        print 'Setting variable ', vv, ' to FALSE because of global option for ', var
                        lopt['OPTIONS'].update({vv: False})

        # if the option is set that the observation time shall be
        # the same as the models
        # then overwrite options that were set in the INI files
        if cfg.same_time4obs:
            for var in thevariables:
                lopt = self.options[var]
                lopt['OPTIONS']['start'] = cfg.start_date
                lopt['OPTIONS']['stop'] = cfg.stop_date

        # map interpolation methods
        # the interpolation method is used by the CDOs. It needs to be
        # a value of [bilinear,conservative,nearest]
        for var in thevariables:
            lopt = self.options[var]
            if lopt['OPTIONS']['interpolation'] == 'bilinear':
                lopt['OPTIONS'].update({'interpolation': 'remapbil'})
            elif lopt['OPTIONS']['interpolation'] == 'conservative':
                lopt['OPTIONS'].update({'interpolation': 'remapcon'})
            elif lopt['OPTIONS']['interpolation'] == 'nearest':
                lopt['OPTIONS'].update({'interpolation': 'remapnn'})
            else:
                raise ValueError('ERROR: invalid interpolation method: \
                                  %s' % lopt['OPTIONS']['interpolation'])

    def _convert_options(self):
        """
        Convert options, which are only strings in the beginning
        to numerical values or execute functions to set directory
        names appropriately
        """
        for v in self.options.keys():
            var = self.options[v]
            for s in var.keys():  # each section
                sec = var[s]
                for k in sec.keys():
                    if k == 'start':
                        sec.update({k: pl.num2date(pl.datestr2num(sec[k]))})
                    elif k == 'stop':
                        sec.update({k: pl.num2date(pl.datestr2num(sec[k]))})
                    else:
                        # update current variable with valid value
                        sec.update({k: self.__convert(sec[k])})

    def __convert(self, s):
        """
        convert a single string into a valid value. The routine
        recognizes if the string is numerical, or boolean
        or if a command needs to be executed to generate a valid value

        Parameters
        ----------
        s : str
            string with value of option
        """
        s.replace('\n', '')
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
            cdo = Cdo()
            cmd = h[1:-1]
            exec('res = ' + cmd)
            return res

        #4) check if a list is provided
        if (h[0] == '[') and (h[-1] == ']'):  # is list
            return self.__get_list(s)

        #5) in any other case return s
        return s

    def __get_list(self, s):
        # return a list made out from a string '[1,2,3,4,5]'
        l = s.replace('[', '')
        l = l.replace(']', '')
        l = l.split(',')

        o = []
        for i in l:
            try:
                o.append(float(i))  # try numerical conversion
            except:
                # else: return the string itself
                o.append(i)
        return o

    def _check(self):
        """
        check consistency of options specified
        """
        cerr = 0
        o = self.options

        # specify here the options that need to be given!
        # variables that need to be specified (MUST!) for each
        # observational dataset
        locopt = ['obs_file', 'obs_var', 'gleckler_position',
                  'scale_data']
        globopt = ['cticks', 'map_difference', 'map_seasons',
                   'preprocess', 'reichler_plot', 'gleckler_plot',
                   'hovmoeller_plot', 'regional_analysis',
                   'interpolation', 'targetgrid', 'projection',
                   'global_mean', 'vmin', 'vmax', 'dmin', 'dmax',
                   'cmin', 'cmax', 'pattern_correlation']  # options for each variable type

        # all variables
        for v in o.keys():
            d = o[v]  # dictionary for a specific variable
            if not 'OPTIONS' in d.keys():
                sys.stdout.write('Error: missing OPTIONS %s' % v)
                cerr += 1

            # check global options
            for k in globopt:
                if not k in d['OPTIONS'].keys():
                    sys.stdout.write('Error: missing global option: %s (%s)' % (k, v))
                    cerr += 1

                if k == 'cticks':
                    if isinstance(d['OPTIONS'][k], list):
                        if np.any(np.diff(d['OPTIONS'][k]) < 0):
                            raise ValueError('CTICKS are not in \
                                               increasing order!')
                    else:
                        raise ValueError('CTICKS option needs to \
                                           be a list')

                if k == 'regional_analysis':
                    if d['OPTIONS'][k] is True:
                        if 'region_file' not in d['OPTIONS'].keys():
                            raise ValueError('ERROR: You need to provide a region file name if '
                                             'you want to use regional_analysis!')

            # check local options
            # odat is key for a specific observational dataset
            for odat in d.keys():
                if odat.upper() == 'OPTIONS':
                    continue
                # k is not the index for a specific obs. record
                for k in locopt:
                    if not k in d[odat].keys():
                        sys.stdout.write('Error: missing local option: %s (%s,%s)' % (k, odat, v))
                        cerr += 1
                    if k == 'obs_file':
                        d[odat]['obs_file'] = d[odat]['obs_file'].rstrip()
                        if ((d[odat]['obs_file'][-1] == os.sep) or (d[odat]['obs_file'][-3:] == '.nc') or (d[odat]['obs_file'][-4:] == '.nc4')):
                            pass
                        else:
                            d[odat]['obs_file'] = d[odat]['obs_file'] + os.sep

        # check if region file is given
        if 'region_file' in d['OPTIONS'].keys():
            if d['OPTIONS']['region_file'].lower() == 'none':
                pass
            elif len(d['OPTIONS']['region_file']) == 0:
                pass
            else:
                if not os.path.exists(d['OPTIONS']['region_file']):
                    print d['OPTIONS']['region_file']
                    raise ValueError('Regional masking file not existing: %s' % d['OPTIONS']['region_file'])
                else:
                    region_filename, region_file_varname = self._import_regional_file(d['OPTIONS']['region_file'], d['OPTIONS'].pop('region_file_varname', None), d['OPTIONS']['targetgrid'])
                    self.options[v]['OPTIONS'].update({'region_file': region_filename})
                    self.options[v]['OPTIONS'].update({'region_file_varname': region_file_varname})

        if cerr > 0:
            raise ValueError('There were errors in the initialization \
                              of plotting options!')


class xxxxxxxxxAnalysisRegions():
    """
    Class to handle information about regions to analyze in the framework
    """

    def __init__(self, dir='.' + os.sep + 'regions' + os.sep):
        """
        read regions

        dir : str
            directory where to search for region specification
        """
        self.dir = dir
        self.regions = []
        self.read()

    def read(self):
        files = glob.glob(self.dir + '*.reg')
        for f in files:
            self._read_region_file(f)

        print '--------'
        print 'Regions:'
        print '--------'
        for r in self.regions:
            print r.label

    def _read_region_file(self, file):
        """
        read a single file specifying regions. A file can contain
        more than one region. Results will be stored in self.regions

        Parameters
        ----------
        file : str
            filename
        """
        if not os.path.exists(file):
            return

        f = open(file, 'r')
        for line in f.readlines():
            l = line.lstrip()
            if l[0] == '#':
                continue
            else:
                tmp = l.split(',')
                if len(tmp) != 5:
                    raise ValueError("Error in REGION file: %s" % file)
                else:
                    # interpret data
                    label = tmp[0]
                    lon1 = tmp[1]
                    lon2 = tmp[2]
                    lat1 = tmp[3]
                    lat2 = tmp[4]
                    # set region
                    self.regions.append(Region(lon1, lon2, lat1,
                                               lat2, label,
                                               type='latlon'))


class CFGWriter(object):
    """
    This class is supposed to provide a standardized interface to write
    a pyCMBS configuration file
    """

    def __init__(self, filename, generator='pyCMBS CONFIGURATION WRITER'):
        """
        Parameters
        ---------
        filename : str
            filename of configuration file that shall be generated
        generator : str
            identifier of the caller of the class
        """
        self.generator = generator
        self.filename = filename
        self.output_dir = os.path.dirname(filename)
        if os.path.exists(filename):
            os.remove(self.filename)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def save(self, temp_dir=None, vars=None, start_date=None,
             stop_date=None, models=None, interval='monthly',
             format='png', basemap=False, clean_temp=False,
             summary_only=False):
        """
        save configuration file
        """
        # list which specifies which default variables should be written
        # to the standard configuration file
        supported_vars = ['albedo', 'sis', 'albedo_vis',
                          'albedo_nir', 'surface_upward_flux', 'tree',
                          'temperature', 'rain', 'evap', 'hair', 'wind',
                          'twpa', 'wvpa', 'late', 'budg',
                          'gpp']
        # 'phenology_faPAR'

        if format.lower() not in ['pdf', 'png']:
            raise ValueError('Invalid output format for report: %s' % format)
        if interval not in ['monthly', 'season']:
            raise ValueError('Invalid interval option specified!')
        if temp_dir is None:
            raise ValueError('No temporary output directory specified!')
        if vars is None:
            raise ValueError('No VARIABLES specified!')
        if start_date is None:
            raise ValueError('No start_date specified')
        if stop_date is None:
            raise ValueError('No stop_date specified')
        if models is None:
            raise ValueError('No models specified!')

        import time
        self._write('######################################################')
        self._write('# AUTOMATICALLY GENERATED configuration file for pyCMBS')
        self._write('# https://code.zmaw.de/projects/pycmbs')
        self._write('# generated by ' + self.generator)
        self._write('#')
        self._write('# generated at: ' + time.asctime())
        self._write('######################################################')

        self._write('basemap,' + str(basemap.real))
        self._write('report=reportname_here')
        self._write('report_format=' + format.upper())
        self._write('author=TheAuthorName')
        self._write('temp_dir=' + temp_dir)
        self._write('clean_tempdir,' + str(clean_temp.real))
        self._write('summary_only,' + str(summary_only.real))
        self._write('config_dir=' + self.output_dir + '/configuration/')
        if not os.path.exists(self.output_dir + '/configuration/'):
            os.makedirs(self.output_dir + '/configuration/')
        self._write('output_directory=' + self.output_dir + '/reports/')
        self._write('')

        self._write('################################')
        self._write('# Specify variables to analyze')
        self._write('#')
        self._write("# comments are by '#'")
        self._write('#')
        self._write('# analysis details for each variable are:')
        self._write('# name, [0,1], [monthly,season]')
        self._write('#')
        self._write("# 'name' specifies the variable name to be analyzed; needs to be consistent with routines defined"
                    " in main()")
        self._write('# [0,1] specified if the data shall be used')
        self._write('# [monthly,season] specifies the temporal scale of the analysis')
        self._write('#')
        self._write('################################')

        for v in supported_vars:
            if v in vars:
                self._write(v + ',1,' + interval)
            else:
                self._write(v + ',0,' + interval)
        self._write('')

        self._write('################################')
        self._write('# specify period to analyze')
        self._write('# start-time YYYY-MM-DD')
        self._write('# stop-time  YYYY-MM-DD')
        self._write('################################')
        self._write(start_date)
        self._write(stop_date)
        self._write('use_for_observations,0')
        self._write('')

        self._write('################################')
        self._write('# Register models to analyze')
        self._write('# ID,TYPE,EXPERIMENET,PATH')
        self._write('#')
        self._write('# ID: unique ID to specify model, for CMIP5 ID is also part of the filenames!')
        self._write('# TYPE: Type of model to be anaylzed (JSBACH_BOT, CMIP5, JSBACH_RAW)')
        self._write('# EXPERIMENT: an experiment identifier')
        self._write('# PATH: directory path where data is located')
        self._write('#')
        self._write('# The modes MUST NOT be separated with whitepsaces at the moment!')
        self._write('################################')
        self._write('')
        self._write('#--- MODELS TO ANALYZE ---')

        for i in xrange(len(models)):
            self._write(models[i]['id'] + ',' + models[i]['type']
                        + ',' + models[i]['experiment'] + ','
                        + models[i]['path'])

    def _write(self, s):
        if os.path.exists(self.filename):
            mode = 'a'
        else:
            mode = 'w'
        f = open(self.filename, mode)
        f.write(s + '\n')
        f.close()
