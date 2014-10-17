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
from scipy import linalg, dot
import matplotlib.gridspec as gridspec
import pickle

from pycmbs.plots import pm_bar, add_nice_legend
from pycmbs.mapping import map_plot
from pycmbs.data import Data
from pycmbs.anova import *
from pycmbs.taylor import Taylor
from pycmbs.benchmarking.report import Report

from pycmbs.plots.violin import ViolinPlot


class RegionalAnalysis(object):
    """
    a class to perform comparisons between two datasets on a regional basis
    """
    def __init__(self, x, y, region, f_correlation=True, f_statistic=True,  f_aggregated_violin=False, report=None):
        """

        Parameters
        ----------
        x : Data
            first dataset (is assumed to be the reference dataset!)
        y : Data
            second dataset
        region : Data
            region to analyze in both datasets; needs to have same geometry
            as the datasets
        f_statistic : bool
            calculate standard first and second moment statistics
        f_correlation : bool
            calculate correlation statistics between datasets
        f_aggregated_violin : bool
            generate a violin plot for the whole statistics over
            all timesteps
        report : Report
            Report class. When provided, figures are integrated
            into the report on the fly
        """

        self.region = region
        self._preprocessed = False
        self.f_correlation = f_correlation
        self.f_statistic = f_statistic
        self.f_aggregated_violin = f_aggregated_violin
        self.statistics = {}
        self.report=report

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
        if self.report is not None:
            if not isinstance(self.report, Report):
                raise ValueError('Report is of invalid object type')
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

        # check geometries
        if self.x.shape != self.y.shape:
            print self.x.shape, self.y.shape
            raise ValueError('ERROR: RegionalAnalyis - \
                               inconsistent geometries!')


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
        vals = self._get_unique_region_ids()
        for v in vals:
            x = self._get_masked_data(self.x, v)
            y = self._get_masked_data(self.y, v)

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

            print('TODO: how to deal with area weighting in global correlation analyis ????')  # TODO

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

        Parameters
        ----------
        pthres : float
            threshold for significance level
            e.g. 0.05 corresponds to the 95% level
            values > 1. mean that no masking is applied to
            insignificant results

        Returns
        -------
        returns a dictionary with details on regional statistics
        """

        if self.f_statistic:  # calculate general statistic
            self._calc_general_statistic()
        if self.f_correlation:
            self._calc_correlation(pthres=pthres)
        if self.f_aggregated_violin:
            self._calc_global_violin()

    def _get_unique_region_ids(self):
        return np.unique(self.region.data.flatten())

    def _get_masked_data(self, x, id):
        """
        mask dataobject for a particular region
        and return masked object. The geometry is not changed.

        Parameters
        ----------
        x : Data
            data to be masked
        id : int
            region ID
        """

        # TODO cut bounding box for better performance when
        # bug in cut_bouding_box has been fixed

        if id not in self._get_unique_region_ids():
            print('WARNING: ID %s not found in data!')
            return x

        msk = self.region.data == id
        d = x.copy()
        d._apply_mask(msk)
        del msk
        return d

    def _calc_pdf_analysis(self):
        """
        do PDF analysis for regions
        """
        vals = self._get_unique_region_ids()
        for v in vals:
            x = self._get_masked_data(self.x, v)
            y = self._get_masked_data(self.y, v)
            x.label = x._get_label() + '-' + ' Region ' + str(v).zfill(5)
            y.label = y._get_label() + '-' + ' Region ' + str(v).zfill(5)
            pdf = PDFAnalysis(x, Y=y, report=self.report, qq=False, labels=None)
            pdf.plot()
            del pdf


    def _calc_global_violin(self):
        """
        generate single violin plot that summarizes information
        for all regions
        """

        assert False # not finally implemented


        xdata = []
        labels=[]
        if self.y is not None:
            ydata = []
        else:
            ydata = None

        vals = self._get_unique_region_ids()
        for v in vals:
            x = self._get_masked_data(self.x, v)
            xdata.append(x.data.flatten())
            labels.append(str(v).zfill(5))  # todo replace by original region labels
            del x

            if self.y is not None:
                y = self._get_masked_data(self.y, v)
                ydata.append(y.data.flatten())
                del y

        print 'Generating Violin plot ...'
        V = ViolinPlot(xdata, ydata, labels=labels)
        V.plot()


        print 'Done !'
        if self.report is not None:
            self.report.figure(V.ax.figure, caption='Violin plot for ' + self.variable.upper())
            plt.close(V.ax.figure.number)
        del V

#~ [ngroups, nsamp]
#~
    #~ def __init__(self, data, data2=None, labels=None, ax=None,
                 #~ boxplot=True, figsize=(10, 6)):
#~
#~
#~
#~
#~ Calculate Violin plot for each region
#~ generate a list with all regions as ID/label and generate a single plot
#~ that plots a vilon plot for each region
#~
#~ In best case the CTRL and EXPERIMENT data is plotted together in the Violin plot!



    def _calc_general_statistic(self):
        """
        calculate standard statistics for each region
        """
        xstat = None
        ystat = None
        if self.x is not None:
            xstat = self.x.condstat(self.region)
        if self.y is not None:
            ystat = self.y.condstat(self.region)
        self.statistics.update({'xstat': xstat})
        self.statistics.update({'ystat': ystat})

    def _calc_correlation(self, pthres=1.01):
        """
        calculate correlation between X and Y in different ways
        """
        ### TODO weighting when calculating correlation
        self.statistics.update({'corrstat': self._get_correlation(pthres=pthres)})

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

        if format not in ['txt', 'pkl', 'tex']:
            raise ValueError('Invalid format for output [txt,pkl]: %s' % format)
        else:
            self.format = format
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
        elif format in ['txt', 'tex']:
            self._save_standard_statistics(dir + prefix + '_regional_statistics_standard.' + format)
            self._save_correlation_statistics_A(dir + prefix + '_regional_statistics_correlation_A.' + format)
            self._save_correlation_statistics_B(dir + prefix + '_regional_statistics_correlation_B.' + format)
            self._save_correlation_statistics_C(dir + prefix + '_regional_statistics_correlation_C.' + format)
        else:
            raise ValueError('Unsupported output format!')

    def _save_correlation_statistics_B(self, fname, tok='B'):
        """
        save correlation B results to ASCII file

        | id | slope | intercept | correlation | pvalue |

        """
        if self.format == 'txt':
            sep = '\t'
            eol = '\n'
        elif self.format == 'tex':
            sep = ' & '
            eol = ' \\\\  \n'

        if os.path.exists(fname):
            os.remove(fname)
        if os.path.splitext(fname)[1] != '.' + self.format:
            fname += '.' + self.format
        o = open(fname, 'w')
        # header
        if self.format == 'tex':
            if self.report is not None:
                self.report.open_table()
            o.write('\\begin{tabular}{lcccc}' + eol)
        o.write('id' + sep + 'slope' + sep + 'intercept' + sep + 'correlation' + sep + 'pvalue' + eol)
        # data
        corrstat = self.statistics['corrstat']['analysis_' + tok]
        for k in corrstat.keys():
            s = str(k) \
                + sep + str(corrstat[k]['slope']) \
                + sep + str(corrstat[k]['intercept']) \
                + sep + str(corrstat[k]['correlation']) \
                + sep + str(corrstat[k]['pvalue']) + eol
            o.write(s)
        if self.format == 'tex':
            if self.report is not None:
                self.report.write('    \\input{' + fname + '}')
                self.report.close_table(caption='Regional statistics, correlation ' + tok)
            o.write('\end{tabular}')
            if self.report is not None:
                self.report.barrier()

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

        if self.format == 'txt':
            sep = '\t'
            eol = '\n'
        elif self.format == 'tex':
            sep = ' & '
            eol = ' \\\\  \n'

        if os.path.exists(fname):
            os.remove(fname)
        if os.path.splitext(fname)[1] != '.' + self.format:
            fname += '.' + self.format
        o = open(fname, 'w')
        # header
        if self.format == 'tex':
            if self.report is not None:
                self.report.open_table()
            o.write('\\begin{tabular}{lccccc}' + eol)
        o.write('id' + sep + 'r-mean' + sep + 'r-std' + sep + 'r-sum' + sep + 'r-min' + sep + 'r-max' + eol)

        # data
        corrstat = self.statistics['corrstat']['analysis_A']
        for k in corrstat.keys():
            s = str(k) \
                + sep + str(corrstat[k]['mean'][0]) \
                + sep + str(corrstat[k]['std'][0]) \
                + sep + str(corrstat[k]['sum'][0]) \
                + sep + str(corrstat[k]['min'][0]) \
                + sep + str(corrstat[k]['max'][0]) + eol
            o.write(s)
        if self.format == 'tex':
            if self.report is not None:
                self.report.write('    \\input{' + fname + '}')
                self.report.close_table(caption='Regional statistics, correlation A')
            o.write('\end{tabular}')
            if self.report is not None:
                self.report.barrier()


        o.close()

    def _save_standard_statistics(self, fname):
        """
        save standard (first and second order) statistics to ASCII files as follows

        time | xmean | ymean | xstd | ystd | xmin | ymin | xmax | ymax

        This information is stored separatley for each region ID in a separate file!

        The output is written dependent on the format.

        Parameters
        ----------
        fname : str
            filename for storing results
        """
        root = os.path.splitext(fname)[0]
        if self.format == 'txt':
            sep = '\t'
            eol = '\n'
        elif self.format == 'tex':
            sep = ' & '
            eol = ' \\\\  \n'
        else:
            raise ValueError('Invalid format')

        ids = self.statistics['xstat'].keys()
        for id in ids:
            fname = root + '_' + str(id).zfill(16) + '.' + self.format
            if os.path.exists(fname):
                os.remove(fname)
            o = open(fname, 'w')

            # header
            if self.format == 'tex':
                if self.report is not None:
                    self.report.open_table()
                o.write('\\begin{tabular}{lcccccccc}' + eol)
            o.write('time' + sep + 'xmean' + sep + 'ymean' + sep + 'xstd' + sep + 'ystd' + sep + 'xmin' + sep + 'ymin' + sep + 'xmax' + sep +'ymax' + eol)

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
                    + str(self.statistics['ystat'][id]['max'][i]) + eol
                o.write(s)
            if self.format == 'tex':
                if self.report is not None:
                    self.report.write('    \\input{' + fname + '}')
                    self.report.close_table(caption='General regional statistics')
                o.write('\end{tabular}')
                if self.report is not None:
                    self.report.barrier()

            o.close()

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
