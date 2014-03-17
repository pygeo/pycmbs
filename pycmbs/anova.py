# -*- coding: utf-8 -*-

"""
This file is part of pyCMBS. (c) 2012-2014
For COPYING and LICENSE details, please refer to the file
COPYRIGHT.md
"""

"""
Module for ANOVA analysis

REFERENCES
==========

1. von Storch and Zwiers, 1999, p. 171ff
"""

from scipy import stats

import numpy as np


class Anova1():
    """
    Analysis of variance class for one-way anova analyis
    """
    def __init__(self, x):
        """

        an alternative code for ANOVA analysis in python can be found
        http://adorio-research.org/wordpress/?p=1102

         x : ndarray
            data array which contains multiple ensemble members
            of time series [nrens,ntime]. Each single row contains
            an entire timeseries. Equal number of samples is assumed
        """

        if x.ndim != 2:
            raise ValueError('Only 2D supported for one-way ANOVA')

        self.x = x.copy()
        self.xmean = x.mean(axis=0)  # ensemble mean
        self.mean = x.mean()  # overall mean

        self.n, self.nt = x.shape

    def one_way_anova(self, verbose=False):
        self._calc_sst()  # calculate different variance components
        self.p = self.get_significance(self._get_f())

        if verbose:
            self.print_results()

    def _get_f(self):
        n = self.ssa / (self.nt - 1)
        d = self.sse / (self.nt * (self.n - 1))
        return n / d

    def print_results(self):
        print

        print 'p  :', self.p
        print 'F  :', self._get_f()
        print 'R2 :', self.get_fractional_variance_explained()
        print 'R2a: ', self.get_fractional_variance_explained(adjust=False)
        print 'SSA: ', self.ssa
        print 'SSE: ', self.sse
        print 'SST: ', self.sst

    def _calc_sst(self):
        if not hasattr(self, 'ssa'):
            self._calc_ssa()
        if not hasattr(self, 'sse'):
            self._calc_sse()

        self.sst = self.ssa + self.sse

    def _calc_ssa(self):
        self.ssa = self.n * sum((self.xmean - self.mean) ** 2)

    def get_significance(self, f):
        p = 1. - stats.f.cdf(f, (self.nt - 1), self.nt * (self.n - 1))
        return p

    def _calc_sse(self):
        tmp = []
        for i in range(self.n):
            tmp.append(sum((self.x[i, :] - self.xmean) ** 2))
        tmp = np.asarray(tmp)
        self.sse = sum(tmp)

    def get_fractional_variance_explained(self, adjust=True):
        if adjust:
            return (self.ssa - self.sse * (self.nt - 1) / (self.nt * (self.n - 1))) / self.sst
        else:
            return self.ssa / self.sst

#-----------------------------------------------------------------------


class Anova2():
    def __init__(self, x):
        '''
        x [blocks,treatmens,nrens]
        '''
        self.x = x.copy()
        self.mean = x.mean()  # overall mean

        # mean value for each block (Ybar_0j0)
        self.bmean = x.mean(axis=1).mean(axis=1)
        # mean value for all treatments (ybar_i00)
        self.tmean = x.mean(axis=0).mean(axis=1)

        # nr of experiments, timesteps, ensemble members
        self.J, self.I, self.n = self.x.shape

    def two_way_anova_with_replication(self, verbose=False):
        self._calc_sst()
        self._calc_f()
        if verbose:
            self.print_results()

    def get_fractional_variance_explained(self, s, adjust=True):
        """
        s: a,b,i,e

        @todo: not absolutely sure about these formulas, need some reference!
        """
        if adjust:
            print 'WARNING: AJUSTED FRACTIONAL VARIANCE NOT VALIDATED WITH SOME REFERENCE DATA'
            if s == 'a':
                adj = (self.sse / self.sst) * (self.df_ssa / self.df_sse)
            elif s == 'b':
                adj = (self.sse / self.sst) * (self.df_ssb / self.df_sse)
            elif s == 'i':
                adj = (self.sse / self.sst) * (self.df_ssi / self.df_sse)
            elif s == 'e':
                adj = (self.sse / self.sst) * (self.df_sse / self.df_sse)
            else:
                raise ValueError('unkown type')

        else:
            adj = 0.

        if s == 'a':
            return self.ssa / self.sst - adj
        elif s == 'b':
            return self.ssb / self.sst - adj
        elif s == 'i':
            return self.ssi / self.sst - adj
        elif s == 'e':
            return self.sse / self.sst - adj
        else:
            raise ValueError('Invalid identified for variance!')

    def get_significance(self, f, df1, df2):
        p = 1. - stats.f.cdf(f, df1, df2)
        return p

    def _calc_sst(self):
        if not hasattr(self, 'ssa'):
            self._calc_ssa()
        if not hasattr(self, 'ssb'):
            self._calc_ssb()
        if not hasattr(self, 'sse'):
            self._calc_sse()
        if not hasattr(self, 'ssi'):
            self._calc_ssi()

        self.sst = self.ssa + self.sse + self.ssb + self.ssi  # eq. 9.19 in von Storch

    def _calc_ssa(self):
        self.ssa = self.n * self.J * sum((self.tmean - self.mean) ** 2)  # eq. 9.21 in von Storch

    def _calc_ssb(self):
        self.ssb = self.n * self.I * sum((self.bmean - self.mean) ** 2)  # eq. 9.22 in von Storch

    def _calc_ssi(self):
        ensmean = self.x.mean(axis=2)  # mean over all replications [J,I]
        s = 0.
        for i in np.arange(self.I):  # eq.9.23
            for j in np.arange(self.J):
                s += (ensmean[j, i] - self.tmean[i] - self.bmean[j] + self.mean) ** 2
        self.ssi = self.n * s

    def _calc_sse(self):
        ensmean = self.x.mean(axis=2)
        #~ print ensmean.shape
        s = []
        for l in np.arange(self.n):
            tmp = (self.x[:, :, l] - ensmean) ** 2
            #~ print tmp.shape, tmp.sum()
            s.append(tmp.sum())
        s = np.asarray(s)
        if len(s) != self.n:
            raise ValueError('invalid lengths!')

        self.sse = sum(s)

    def _calc_f(self):
        #- calculate degree of freedoms
        self.df_ssa = (self.I - 1)
        self.df_ssb = (self.J - 1)
        self.df_ssi = (self.I - 1) * (self.J - 1)
        self.df_sse = self.I * self.J * (self.n - 1)

        #- calculate MS
        mssa = self.ssa / self.df_ssa
        mssb = self.ssb / self.df_ssb
        mssi = self.ssi / self.df_ssi
        msse = self.sse / self.df_sse

        #- calculate F-value (caution: not derived from von Storch & Zwiers, but from http://people.richland.edu/james/lecture/m170/ch13-2wy.html)
        #@todo: clarify !
        self.f_ssa = mssa / msse
        self.f_ssb = mssb / msse
        self.f_ssi = mssi / msse

        self.p_ssa = self.get_significance(self.f_ssa, self.df_ssa, self.df_sse)
        self.p_ssb = self.get_significance(self.f_ssb, self.df_ssb, self.df_sse)
        self.p_ssi = self.get_significance(self.f_ssi, self.df_ssi, self.df_sse)

    def print_results(self):
        print 'Source of variation', 'SS', 'df', 'F', 'p'
        print 'A          : ', self.ssa, self.df_ssa, self.f_ssa, self.p_ssa
        print 'B          : ', self.ssb, self.df_ssb, self.f_ssb, self.p_ssb
        print 'interaction: ', self.ssi, self.df_ssi, self.f_ssi, self.p_ssi
        print 'within     : ', self.sse, self.df_sse
        print 'total      : ', self.sst

#-----------------------------------------------------------------------


def __example_one_way():
    '''
    http://adorio-research.org/wordpress/?p=1102
    '''

    groups = [[48, 49, 50, 49],
              [47, 49, 48, 48],
              [49, 51, 50, 50]]
    groups = np.asarray(groups).T  # transpose necessary, as reference routine works different

    A = Anova1(groups)
    A.one_way_anova()

    print self.p
    print self.get_fractional_variance_explained()
    print self.get_fractional_variance_explained(adjust=False)
    print self.ssa
    print self.sse
    print self.sst


#http://adorio-research.org/wordpress/?p=11857
#~ groups = [
           #~ [[64, 72, 74],
            #~ [66, 81, 51],
            #~ [70, 64, 65]
           #~ ],
           #~ [[65, 57, 47],
            #~ [63, 43, 58],
            #~ [58, 52, 67]
           #~ ],
           #~ [[59, 66, 58],
            #~ [68, 71, 39],
            #~ [65, 59, 42]
           #~ ],
           #~ [[58, 57, 53],
            #~ [41, 61, 59],
            #~ [46, 53, 38]
           #~ ]
         #~ ]
#~
#~ groups = [
          #~ [[6,4,5,5,4],
           #~ [5,7,4,6,8]
          #~ ],
          #~ [[10,8,7,7,9],
           #~ [7,9,12,8,8]
          #~ ],
          #~ [[7,5,6,5,9],
           #~ [9,7,5,4,6]
          #~ ],
          #~ [[8,4,6,5,5],
           #~ [5,7,9,7,10]
          #~ ]
         #~ ]
#~
#~
#~
#~ groups = [
          #~ [[4,7,10],
          #~ [6,13,12]],
#~
          #~ [[5,9,12],
          #~ [6,15,13]],
#~
          #~ [[6,8,11],
          #~ [4,12,10]],
#~
          #~ [[5,12,9],
          #~ [4,12,13]]
#~
          #~ ]
#~
#~
#~ #http://people.richland.edu/james/lecture/m170/ch13-2wy.html
#~ groups = [
          #~ [[106,95,94,103,100],
          #~ [110,98,100,108,105],
          #~ [94,86,98,99,94]],
#~
          #~ [[110,100,107,104,102],
          #~ [112,99,101,112,107],
          #~ [97,87,99,101,98]]
#~
          #~ ]
#~
#~
#~
#~ ny,nx = np.shape(groups[0])
#~ nr = len(groups)
#~
#~ trans = False #do transpose
#~
#~ if trans:
    #~ g1=np.zeros((nx,ny,nr))*np.nan
#~ else:
    #~ g1=np.zeros((ny,nx,nr))*np.nan
#~ for i in range(len(groups)):
    #~ if trans:
        #~ g1[:,:,i] = np.asarray(np.asarray(groups[i])).T
    #~ else:
        #~ g1[:,:,i] = np.asarray(np.asarray(groups[i]))
#~
#~
#~ A = Anova2(g1)
#~ A.two_way_anova_with_replication()
#~ A._calc_f()
#~
#~ print 'SSA: ', A.ssa, A.df_ssa
#~ print 'SSB: ', A.ssb, A.df_ssb
#~ print 'SSI: ', A.ssi, A.df_ssi
#~ print 'SST: ', A.sst
#~ print 'SSE: ', A.sse, A.df_sse
#~ print 'n  : ', A.n
#~ print 'I  : ', A.I
#~ print 'J  : ', A.J
#~
#~
#~ print A.p_ssa
#~ print A.p_ssb
#~ print A.p_ssi
