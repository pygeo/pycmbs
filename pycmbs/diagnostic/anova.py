# -*- coding: utf-8 -*-
"""
This file is part of pyCMBS.
(c) 2012- Alexander Loew
For COPYING and LICENSE details, please refer to the LICENSE file
"""


class ANOVA(object):
    """
    main class to perform an ANOVA analysis
    using C{Data} objects

    treatments are derived from data timeseries
    blocks are obtained by different experiments
    """
    def __init__(self):
        self.experiments = []
        self.data = {}

    def add_experiment(self, label):
        self.experiments.append(self.__trim(label))
        self.data.update({self.__trim(label): []})

    def __trim(self, s):
        return s.replace(' ', '_')

    def add_data(self, e, d):
        """
        adds data for each experiment to a list
        """
        k = self.__trim(e)  # experiment key
        if k in self.experiments:
            #experiment was registered
            #~ has_mask = True
            #~ try: #this is a hack; how to do it better ??? #TODO !!!!
                #~ d.data.shape == d.mask.shape
            #~ except:
                #~ has_mask = False
            #~ if not has_mask:
                #~ print d.data
                #~ print e
                #~ raise ValueError, 'All data objects need to me masked arrays!' #this is needed,as
                # otherwise common mask generation is not possible
            self.data[k].append(d)
        else:
            raise ValueError('Experiment was not yet registered!')

    def analysis(self, analysis_type=None):
        """
        perform ANOVA analysis based on the data given
        """

        if analysis_type is None:
            raise ValueError('Needs to specify type of ANOVA to be performed')

        #1) check if all data has the same length
        # this is not a strong prerequesite for ANOVA as such, but for
        # the implementation here
        # generate a common mask with points that are valid throughout
        # all timesteps and are not masked
        self._get_valid_mask()  # --> self.mask

        #2) rearange data to work only with valid (not masked data)
        #   check also that only data which is valid in all cases is used
        #   it is realized by getting the valid indices of the data
        idx = np.argwhere(self.mask)  # returns list of indices of valid pixels

        #3) perform ANOVA analysis on all pixels individually
        resA = np.zeros(self.refshape) * np.nan
        resB = np.zeros(self.refshape) * np.nan
        resI = np.zeros(self.refshape) * np.nan
        resE = np.zeros(self.refshape) * np.nan
        resPA = np.zeros(self.refshape) * np.nan
        resPB = np.zeros(self.refshape) * np.nan
        resPI = np.zeros(self.refshape) * np.nan

        resAa = np.zeros(self.refshape) * np.nan
        resBa = np.zeros(self.refshape) * np.nan
        resIa = np.zeros(self.refshape) * np.nan

        for p in idx:
            if analysis_type == 'one':  # one way ANOVA
                m = self._data2anova1(p)  # [nrens,ntime]
                A = Anova1(m)
                A.one_way_anova(verbose=False)

                resA[p[0], p[1]] = A.get_fractional_variance_explained(adjust=False)  # ssa
                resPA[p[0], p[1]] = A.p
                resB[p[0], p[1]] = A.sse / A.sst  # todo: adjustment here ???

            elif analysis_type == 'two':  # two way ANOVA
                m = self._data2anova2(p)
                #- perform 2-way anova
                A = Anova2(m)
                A.two_way_anova_with_replication()

                resA[p[0], p[1]] = A.get_fractional_variance_explained('a', adjust=False)  # todo: adjust variance
                resB[p[0], p[1]] = A.get_fractional_variance_explained('b', adjust=False)  # todo: adjust variance
                resI[p[0], p[1]] = A.get_fractional_variance_explained('i', adjust=False)  # todo: adjust variance
                resE[p[0], p[1]] = A.get_fractional_variance_explained('e', adjust=False)  # todo: adjust variance

                resPA[p[0], p[1]] = A.p_ssa
                resPB[p[0], p[1]] = A.p_ssb
                resPI[p[0], p[1]] = A.p_ssi

                #~ resAa[p[0],p[1]] = A.get_fractional_variance_explained('a',adjust=True) #todo: adjust variance
                #~ resBa[p[0],p[1]] = A.get_fractional_variance_explained('b',adjust=True) #todo: adjust variance
                #~ resIa[p[0],p[1]] = A.get_fractional_variance_explained('i',adjust=True) #todo: adjust variance

                #@todo significance

            else:
                raise ValueError('Invalid ANOVA type')

        self.resA = resA
        self.resB = resB
        self.resI = resI
        self.resE = resE

        self.resPA = resPA
        self.resPB = resPB
        self.resPI = resPI

        self.resAa = resAa
        self.resBa = resBa
        self.resIa = resIa

    def _data2anova1(self, p):
        """
        extract from the database all the data
        relevant for a single location, given by the indices in p

        p = indices

        ---> time (nt)
        |
        |
        v experiment (nexp)
        """

        nexp = len(self.experiments)

        if nexp != 1:
            raise ValueError('one-way anova only valid for signle experiments!')

        x = np.zeros((self.n, self.nt)) * np.nan  # [nrens,nt]

        e = self.experiments[0]
        for i in range(self.n):
            d = self.data[e][i]  # data object for experiment 'e' and ensemble nr [i]
            x[i, :] = d.data[:, p[0], p[1]]

        if np.any(np.isnan(x) > 0):
            raise ValueError('Something is wrong: not all data valid!')

        return x

    def _data2anova2(self, p):
        """
        extract from the database all the data
        relevant for a single location, given by the indices in p

        p = indices

        ---> time (nt)
        |
        |
        v experiment (nexp)
        """
        nexp = len(self.experiments)
        x = np.zeros((nexp, self.nt, self.n)) * np.nan

        for j in range(len(self.experiments)):
            e = self.experiments[j]
            for i in range(self.n):
                d = self.data[e][i]  # data object for experiment 'e' and ensemble nr [i]
                x[j, :, i] = d.data[:, p[0], p[1]]

        if np.any(np.isnan(x) > 0):
            raise ValueError('Something is wrong: not all data valid!')

        return x

    def _get_valid_mask(self):
        """
        generate a mask where all datasets are valid
        """
        #- check if geometry in general o.k.
        self.__same_geometry()

    def __same_geometry(self):
        """
        check if all data has the same geometry
        """
        cnt = 0
        for k in self.experiments:
            d = self.data[k]
            if cnt == 0:
                refshape = d[0].data.shape  # first dataset geometry as reference
                self.refshape = (refshape[1], refshape[2])
                nrens = len(d)
                self.n = nrens
                self.nt = refshape[0]
                refmsk = np.ones((refshape[1], refshape[2])).astype('bool')
                cnt += 1

            if len(d) != nrens:
                raise ValueError('Invalid Number of ensemble members found!')
            for i in range(len(d)):

                if d[i].data.shape != refshape:
                    print k, i
                    raise ValueError('Invalid data shape found!')

                #- generate common mask
                msk = d[i].get_valid_mask()
                #~ print sum(msk), k, i
                refmsk = refmsk & msk

        self.mask = refmsk

        return True
