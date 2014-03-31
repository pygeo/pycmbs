# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

import numpy as np
import os
import scipy as sci
from scipy import stats

from matplotlib import pylab as plt

from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.axes as maxes
import matplotlib.cm as cm
import matplotlib.colors as col
import matplotlib as mpl
from pycmbs.plots import pm_bar, add_nice_legend
from pycmbs.mapping import map_plot
from pycmbs.data import Data
from scipy import linalg, dot
import matplotlib.gridspec as gridspec
from pycmbs.anova import *
from pycmbs.taylor import Taylor
### from pylab import *
import pickle

class RegionalAnalysis(object):
    """
    a class to perform comparisons between two datasets on a regional basis
    """
    def __init__(self, x, y, region, f_standard=True, f_correlation=True):
        """

        Parameters
        ----------
        x : Data
            first dataset (is assumed to be the reference dataset!)
        y : Data
            second dataset
        region : Data
            region to analyze in both datasets; needs to
        f_standard : bool
            calculate standard first and second moment statistics
        f_correlation : bool
            calculate correlation statistics between datasets
        """

        self.region = region
        self._preprocessed = False
        self.f_standard = f_standard
        self.f_correlation = f_correlation
        self.statistics = {}

        if x is not None:
            if not isinstance(x, Data):
                raise ValueError('Error: RegionalAnalysis - X \
                                    is not of type Data!')
        if y is not None:
            if not isinstance(y, Data):
                raise ValueError('Error: RegionalAnalysis - Y \
                               is not of type Data!')

        if (x is None) or (y is None):
            # in case of dummy data, do not perform data check
            if x is None:
                self.x = None
            else:
                self.x = x.copy()
            if y is None:
                self.y = None
            else:
                self.y = y.copy()
        else:
            self.x = x.copy()
            self.y = y.copy()
            self._check()  # check input data

#---

    def _check(self):
        """
        check consistency of data
        """
        if not isinstance(self.region, Data):
            raise ValueError('Region needs to be of instance C{Data}!')
        if self.x is not None:
            if self.x.ndim == 2:
                if self.x.shape != self.region.shape:
                    print self.x.shape, self.region.shape
                    raise ValueError('Inconsistent shape of \
                                       X-data with Region')
            elif self.x.ndim == 3:
                if self.x.data[0, :, :].shape != self.region.shape:
                    print self.x.data[0, :, :].shape, self.region.shape
                    raise ValueError('Inconsistent shape of \
                                       X-data with Region')
            else:
                raise ValueError('Unknown geometry!')

        if self.y is not None:
            if self.y.ndim == 2:
                if self.y.shape != self.region.shape:
                    print self.y.shape, self.region.shape
                    raise ValueError('Inconsistent shape of \
                                       Y-data with Region')
            elif self.y.ndim == 3:
                if self.y.data[0, :, :].shape != self.region.shape:
                    print self.y.data[0, :, :].shape, self.region.shape
                    raise ValueError('Inconsistent shape of \
                                      Y-data with Region')
            else:
                raise ValueError('Unknown geometry!')

        #--- check datatypes
        #if not isinstance(self.region,Region):
        #    raise ValueError, 'Error: RegionalAnalysis
        #- region is not of type Region!'

        #--- check geometries
        if self.x.shape != self.y.shape:
            print self.x.shape, self.y.shape
            raise ValueError('ERROR: RegionalAnalyis - \
                               inconsistent geometries!')

#---

    def xxxxx_prepare(self):

        #--- apply regional mask
        if self.region.type == 'latlon':
            self.x.get_aoi_lat_lon(self.region)
            self.y.get_aoi_lat_lon(self.region)
        else:
            raise ValueError('RegionalAnalysis does not work with \
                              regions that are not specified by lat/lon!')
            self.x = self.x.get_aoi(self.region)
            self.y = self.y.get_aoi(self.region)

        #--- reduce data volume by cutting unnecessary data
        self.x = self.x.cut_bounding_box(return_object=True)
        self.y = self.y.cut_bounding_box(return_object=True)

        #--- temporal mean if desired
        if self.use_mean:
            self.x = self.x.timmean(return_object=True)
            self.y = self.y.timmean(return_object=True)
        self._preprocessed = True

#---

    def _get_correlation(self, pthres=1.01):
        """
        calculate correlation between two fields

        in general one can think of three ways to calculate the correlation between X and Y
        a) correlate per pixel --> map of correlation measures; calculate then regional mean of this skill score (e.g. R)
        b) use multidimensional array [time,ny,nx] and correlate whole dataset for a region, without averaging
        c) calculate spatial mean field for each region and then correlate the means

        Here we implement all three approaches to allow for comparability

        Parameters
        ----------
        pthres : float
            threshold to mask insignificant correlations (default>1. ==> no masking)
        """

        if (self.x is None) or (self.y is None):
            return {'analysis_A': None, 'analysis_B': None, 'analysis_C': None}

        #x.get_valid_data

        #=======================================================================
        # A) calculate once correlation and then calculate regional statistics
        RO, PO = self.x.correlate(self.y, pthres=pthres,
                                  spearman=False, detrend=False)
        corrstat1 = RO.condstat(self.region)  # gives a dictionary already

        correlations = []
        correlations1 = []
        slopes = []
        slopes1 = []
        pvalues = []
        pvalues1 = []
        intercepts = []
        intercepts1 = []
        ids = []
        stdx = []
        stdy = []
        vals = np.unique(self.region.data.flatten())
        for v in vals:
            # print('Regional analysis - correlation for ID: %s ' % str(v).zfill(3))
            msk = self.region.data == v  # generate mask
            x = self.x.copy()
            y = self.y.copy()
            x._apply_mask(msk)
            y._apply_mask(msk)
            del msk

            #=======================================================================
            # B) calculate regional statistics based on entire dataset for a region
            xvec = x.data.flatten()
            yvec = y.data.flatten()
            slope, intercept, r_value, p_value, std_err = stats.mstats.linregress(xvec, yvec)
            ids.append(v)
            slopes.append(slope)
            correlations.append(r_value)
            pvalues.append(p_value)
            intercepts.append(intercept)
            stdx.append(xvec.std())
            stdy.append(yvec.std())
            del xvec, yvec, slope, intercept, r_value, p_value, std_err

            #=======================================================================
            # C) fldmean() for each region and then correlate
            xm = x.fldmean(return_data=True)
            ym = y.fldmean(return_data=True)
            sh = xm.shape
            if xm.ndim != 3:
                raise ValueError('Invalid shape: %s' % sh)
            if sh[0] != xm.nt:
                raise ValueError('Timeseries should be of dimension [nt,1,1], %s' % str(sh))
            if sh[1] != 1:
                raise ValueError('Timeseries should be of dimension [nt,1,1], %s' % str(sh))
            if sh[2] != 1:
                raise ValueError('Timeseries should be of dimension [nt,1,1], %s' % str(sh))
            slope1, intercept1, r_value1, p_value1, std_err1 = stats.mstats.linregress(xm.data[:, 0, 0],
                                                                                       ym.data[:, 0, 0])

            slopes1.append(slope1)
            correlations1.append(r_value1)
            pvalues1.append(p_value1)
            intercepts1.append(intercept1)

            print('TODO: how to deal with area weighting in global correlation analyis ????')

            del x, y, xm, ym
        ids = np.asarray(ids)
        slopes = np.asarray(slopes)
        slopes1 = np.asarray(slopes1)
        correlations = np.asarray(correlations)
        correlations1 = np.asarray(correlations1)
        pvalues = np.asarray(pvalues)
        pvalues1 = np.asarray(pvalues1)
        intercepts = np.asarray(intercepts)
        intercepts1 = np.asarray(intercepts1)
        stdx = np.asarray(stdx)
        stdy = np.asarray(stdy)

        def _reshuffle(d):
            """reshuffle structure of output dictionary"""
            r = {}
            for i in xrange(len(d['id'])):
                id = d['id'][i]
                r.update({id: {'slope': d['slope'][i], 'intercept': d['intercept'][i],
                               'correlation': d['correlation'][i], 'pvalue': d['pvalue'][i]}})
                #stdx, stdy are not remapped at the moment
            return r

        corrstat2 = {'id': vals, 'slope': slopes,
                     'correlation': correlations,
                     'pvalue': pvalues,
                     'intercept': intercepts,
                     'stdx': stdx,
                     'stdy': stdy}

        corrstat3 = {'id': vals, 'slope': slopes1,
                     'correlation': correlations1,
                     'pvalue': pvalues1,
                     'intercept': intercepts1}

        #--- return result ---
        return {'analysis_A': corrstat1, 'analysis_B': _reshuffle(corrstat2), 'analysis_C': _reshuffle(corrstat3)}

    def calculate(self, pthres=1.01):
        """
        perform calculation of regional statistics

        The routine calculates the following statistics for EACH regions
        a) timeseries of first and second order moments (mean, std)
        b)

        Returns
        -------
        returns a dictionary with details on regional statistics
        """

        #--- 1) standard statistic for X and Y datasets ---
        xstat = None
        ystat = None
        if self.f_standard:
            if self.x is not None:
                # returns a dictionary with statistics for each region (could be for all timesteps!)
                xstat = self.x.condstat(self.region)
            if self.y is not None:
                ystat = self.y.condstat(self.region)

        self.statistics.update({'xstat': xstat})
        self.statistics.update({'ystat': ystat})

        #--- 2) correlation statistics ---
        corrstat = None
        if self.f_correlation:
            self.statistics.update({'corrstat': self._get_correlation(pthres=pthres)})

        #--- 3) weighted squared difference --> Reichler index for different regions !
        #todo: how to do the weighting ????
        #stop

#---

    def save(self, prefix='', format='txt', dir=None):
        """
        save statistic results to ASCII files

        Parameters
        ----------
        prefix : str
            prefix for output filenames
        format : str
            output format ['pkl','txt']
        dir : str
            directory where results shall be written to
        """

        if format not in ['txt', 'pkl']:
            raise ValueError('Invalid format for output [txt,pkl]: %s' % format)
        if dir is None:
            raise ValueError('You need to specify an output directory!')
        if not os.path.exists(dir):
            os.makedirs(dir)  # try to generate directory
        if not os.path.exists(dir):
            raise ValueError('ERROR: output directory not existing and it also can not be created!')
        if dir[-1] != os.sep:
            dir += os.sep

        if format == 'pkl':
            oname = dir + prefix + '_regional_statistics.pkl'
            if os.path.exists(oname):
                os.remove(oname)
            pickle.dump(self.statistics, open(oname, 'w'))
        elif format == 'txt':
            self._save_standard_statistics(dir + prefix + '_regional_statistics_standard.txt')
            self._save_correlation_statistics_A(dir + prefix + '_regional_statistics_correlation_A.txt')
            self._save_correlation_statistics_B(dir + prefix + '_regional_statistics_correlation_B.txt')
            self._save_correlation_statistics_C(dir + prefix + '_regional_statistics_correlation_C.txt')
        else:
            raise ValueError('Unsupported output format!')

    def _save_correlation_statistics_B(self, fname, tok='B'):
        """
        save correlation B results to ASCII file

        | id | slope | intercept | correlation | pvalue |

        """
        sep = '\t'
        if os.path.exists(fname):
            os.remove(fname)
        if os.path.splitext(fname)[1] != '.txt':
            fname += '.txt'
        o = open(fname, 'w')
        # header
        o.write('id' + sep + 'slope' + sep + 'intercept' + sep + 'correlation' + sep + 'pvalue' + '\n')
        # data
        corrstat = self.statistics['corrstat']['analysis_' + tok]
        for k in corrstat.keys():
            s = str(k) \
                + sep + str(corrstat[k]['slope']) \
                + sep + str(corrstat[k]['intercept']) \
                + sep + str(corrstat[k]['correlation']) \
                + sep + str(corrstat[k]['pvalue']) + '\n'
            o.write(s)
        o.close()

    def _save_correlation_statistics_C(self, fname):
        """
        save correlation B results to ASCII file

        | id | slope | intercept | correlation | pvalue |

        This routine just acts as a wrapper
        """
        self._save_correlation_statistics_B(fname, tok='C')

    def _save_correlation_statistics_A(self, fname):
        """
        save correlation A results to ASCII file
        The results correspond to e.g the *mean* correlation coefficient for a particular region

        | id | rmean | rstd | rsum | rmin | rmax |

        """
        sep = '\t'
        if os.path.exists(fname):
            os.remove(fname)
        if os.path.splitext(fname)[1] != '.txt':
            fname += '.txt'
        o = open(fname, 'w')
        # header
        o.write('id' + sep + 'r-mean' + sep + 'r-std' + sep + 'r-sum' + sep + 'r-min' + sep + 'r-max' + '\n')

        # data
        corrstat = self.statistics['corrstat']['analysis_A']
        for k in corrstat.keys():
            s = str(k) \
                + sep + str(corrstat[k]['mean'][0]) \
                + sep + str(corrstat[k]['std'][0]) \
                + sep + str(corrstat[k]['sum'][0]) \
                + sep + str(corrstat[k]['min'][0]) \
                + sep + str(corrstat[k]['max'][0]) + '\n'
            o.write(s)
        o.close()

    def _save_standard_statistics(self, fname):
        """
        save standard (first and second order) statistics to ASCII files as follows

        time | xmean | ymean | xstd | ystd | xmin | ymin | xmax | ymax

        This information is stored separatley for each region ID in a separate file!

        """
        root = os.path.splitext(fname)[0]
        sep = '\t'

        ids = self.statistics['xstat'].keys()
        for id in ids:
            fname = root + '_' + str(id).zfill(16) + '.txt'
            if os.path.exists(fname):
                os.remove(fname)
            o = open(fname, 'w')

            # header
            o.write('time' + sep + 'xmean' + sep + 'ymean' + sep + 'xstd' + sep + 'ystd' + '\n')

            # data
            for i in xrange(len(self.statistics['xstat'][id]['mean'])):
                s = str(i) + sep \
                    + str(self.statistics['xstat'][id]['mean'][i]) + sep \
                    + str(self.statistics['ystat'][id]['mean'][i]) + sep \
                    + str(self.statistics['xstat'][id]['std'][i]) + sep \
                    + str(self.statistics['ystat'][id]['std'][i]) + sep \
                    + str(self.statistics['xstat'][id]['min'][i]) + sep \
                    + str(self.statistics['ystat'][id]['min'][i]) + sep \
                    + str(self.statistics['xstat'][id]['max'][i]) + sep \
                    + str(self.statistics['ystat'][id]['max'][i]) + '\n'
                o.write(s)
            o.close()

    def xxxxxxxxxxxxxprint_table(self, format='txt', filename=None):
        """
        print table of regional statistics
        """

        def _get_string(d, id):
            #d: dictionary (self.statistics)
            #id region id

            #xmean,xstd,ymean,ystd

            sep = '\t'
            s = str(id) + sep

            #xstat and ystat needs to be printed separately, as it can be a function of time !!!
            # xstat = d['xstat']
            # if xstat == None:
            #     s += '' + sep + '' + sep
            # else:
            #     m = d['xstat']['id'] == id
            #     s += str(d['xstat']['mean'][m]) + sep + str(d['xstat']['std'][m]) + sep
            #
            # ystat = d['ystat']
            # if ystat == None:
            #     s += '' + sep + '' + sep
            # else:
            #     m = d['ystat']['id'] == id
            #     s += str(d['ystat']['mean'][m]) + sep + str(d['ystat']['std'][m]) + sep

            # id, mean correlation, std of correlation
            if d['corrstat']['corrstat1'] is None:
                s += '' + sep + '' + sep
            else:
                stat = d['corrstat']['corrstat1']
                m = stat['id'] == id
                if len(m) > 0:
                    if sum(m) == 1:
                        s += str(stat['mean'][0][m][0]) + sep + str(stat['std'][0][m][0]) + sep
                    else:
                        print id
                        print m
                        raise ValueError('The ID seems not to be unique here!')
                else:
                    s += '' + sep + '' + sep

            #--- corrstat2 ---
            if d['corrstat']['corrstat2'] is None:
                s += '' + sep + '' + sep + '' + sep + '' + sep + '' + sep + '' + sep
            else:
                stat = d['corrstat']['corrstat2']
                m = stat['id'] == id

                if len(m) > 0:
                    if sum(m) == 1:
                        s += str(stat['correlation'][m][0]) + sep + str(stat['pvalue'][m][0]) \
                            + sep + str(stat['slope'][m][0]) + sep + str(stat['intercept'][m][0]) \
                            + sep + str(stat['stdx'][m][0]) + sep + str(stat['stdy'][m][0]) + sep
                    else:
                        print id
                        print m
                        raise ValueError('The ID seems not to be unique here!')
                else:
                    s += '' + sep + '' + sep + '' + sep + '' + sep + '' + sep + '' + sep
            return s

        if filename is not None:
            if os.path.exists(filename):
                os.remove(filename)
            o = open(filename, 'w')

        #--- loop over all regions ---
        keys = np.unique(self.region.data.flatten())
        keys.sort()
        sep = '\t'
        header = 'id' + sep + 'r1' + sep + 'sig_r1' + sep + 'r2' + sep + 'pval2' + sep + 'slope2'\
                 + sep + 'intercept2' + sep + 'stdx' + sep + 'stdy' + sep
        print header
        if filename is not None:
            o.write(header + '\n')
        for k in keys:
            s = _get_string(self.statistics, k)  # get formatted string for a particular region
            print s
            if filename is not None:
                if format == 'txt':
                    o.write(s + '\n')
                else:
                    raise ValueError('Unsupported output format!')

        if filename is not None:
            o.close()

#---

    def plot_taylor(self, dia=None, color='red'):
        """
        Taylor plot of statistics

        The plot produced contains the IDs of each region as default.
        todo: as an alternative one should be able to provide a dictionary that specifies how to plot results
        (e.hg. labels, markerstyles etc.)

        requires that statistics have alredy been calculated

        Parameters
        ----------
        color : str
            color for specifiying the current plot
        dia : Taylor
            Taylor plot instance
        """

        # check
        keys = np.unique(self.region.data.flatten())
        keys.sort()

        r = self.statistics['corrstat']['corrstat2']['correlation']
        sx = self.statistics['corrstat']['corrstat2']['stdx']
        sy = self.statistics['corrstat']['corrstat2']['stdy']

        ratio = sy / sx

        if dia is None:
            tay = Taylor(stdmax=max(ratio.max() * 1.2, 2.))
        else:
            if not isinstance(dia, Taylor):
                print type(dia)
                raise ValueError('Provided argument is no taylor class! ')
            else:
                tay = dia
                tay.stdmax = max(max(ratio.max() * 1.2, 2.), dia.stdmax)  # preserve stdmax information when possible

        sid = map(str, self.statistics['corrstat']['corrstat2']['id'])
        tay.plot(r, ratio, labels=sid, color=color)
        return tay
