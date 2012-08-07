#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"
__email__ = "alexander.loew@zmaw.de"

'''
main module for diagnostic routines
@todo: check routine performances e.g. Reichler based on reference data solutions
'''

import numpy as np

import scipy as sci

from matplotlib import pylab as plt

from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes
import matplotlib as  mpl

from pyCMBS.plots import map_plot
from pyCMBS.data import Data

from scipy import linalg, dot;

import matplotlib.gridspec as gridspec

from anova import *

#-----------------------------------------------------------------------

class EOF():
    '''
    main class to perform an EOF analyis
    '''

    def __init__(self,x,allow_gaps=False):
        '''
        constructor for EOF analysis

        @param x: C{Data} object with a 3D data. The data is assumed to have structire [time,ny,nx]
        @type x: C{Data} object

        @param allow_gaps: specifies if data gaps are allowed. If True, then temporal gaps are allowed
                           and are considered approprately when calculating the covariance matrix
                           if FALSE, then only timeseries without any gaps will be used.
        @type allow_gaps: bool

        @todo: how to deal with negative eigenvalues, which sometimes occur?
        '''

        #/// check geometries
        if x.data.ndim != 3:
            raise ValueError, 'EOF analysis currently only supported for 3D data matrices of type [time,ny,nx]'

        self._x0 = x #preserve information on original data

        #/// reshape data [time,npoints] ///
        self._shape0 = x.data[0,:,:].shape #original data shape
        n = len(x.data) #number of timestamps
        self.n = n

        #-estimate only valid data, discard any masked values
        if allow_gaps:
            lon,lat,vdata,msk = x.get_valid_data(return_mask=True,mode='one')
        else:
            lon,lat,vdata,msk = x.get_valid_data(return_mask=True,mode='all')

        self._x0mask = msk.copy() #store mask applied to original data

        self.x = vdata.copy()
        self.x.shape = (n,-1)
        self.x = self.x.T #[npoints,time]

        print 'Shape is now ...', self.x.shape

        #/// calculate covariance matrix ///
        if allow_gaps:
            self.C = np.ma.cov(self.x,rowvar=0)
        else:
            self.C = np.cov(self.x,rowvar=0)

        #/// solve eigenvalue problem ///
        eigval,eigvec = np.linalg.eig(self.C)

        self.eigval = eigval
        self.eigvec = eigvec

        #/// calculate expansion coefficients == PC (projection of original data to new parameter space)
        if allow_gaps:
            A = np.ma.dot(self.x,eigvec)
        else:
            A = np.dot(self.x,eigvec)
        self.EOF = A

        #/// explained variance
        evar = eigval/sum(eigval)
        self._var = evar #explained variance

    def get_explained_variance(self):
        '''
        returns vector with explained variance
        '''
        return self._var

    def plot_eof_coefficients(self,k,all=False,norm=True,ax=None,label=None,show_legend=True):
        '''
        plot EOF coefficients = time series

        @param k: list of eof coefficients to be plotted
        @type k: list

        @param all: plot all principle components (overwrites k)
        @type all: bool

        @param norm: normalize coefficients by stdv. to allow better plotting (default=True)
        @type norm: bool

        '''
        if all:
            k=range(self.n)
        else:
            if np.isscalar(k):
                k=[k]

        if ax == None:
            f = plt.figure()
            ax = f.add_subplot(111)
        else:
            f=ax.figure

        if label == None:
            label=''
        else:
            label = label + ' '

        for i in k:
            y = self.eigvec[:,i].copy()
            if norm:
                y = (y-y.mean() ) / y.std()
            ax.plot(plt.num2date(self._x0.time),y,label=label + 'EOF'+str(i+1).zfill(3)) #caution: labeling is k+1

        if show_legend:
            ax.legend()

        return ax

    def plot_EOF(self,k,all=False,use_basemap=False,logplot=False,ax=None,label=None,region=None,vmin=None,vmax=None,show_coef=False,cmap=None,title=None):
        '''
        plot multiple eof patterns

        @param k: scalar or list with principal component indices
        @type k: scalar or list

        @param all: plot all principle components (overwrites k)
        @type all: bool

        @param logplot: take log of data for plotting
        @param logplot: bool

        @param show_coef: show coefficients in a separate plot
        @param show_coef: bool
        '''

        if all:
            k = range(self.n)
            ax = None
        else:
            if np.isscalar(k):
                k=[k]

        if show_coef:
            f = plt.figure()
            gs = gridspec.GridSpec(2, 1, wspace=0.05,hspace=0.05,bottom=0.2,height_ratios = [5,1])
            ax  = f.add_subplot(gs[0])
            ax2 = f.add_subplot(gs[1])

        for i in k:
            self._plot_single_EOF(i,use_basemap=use_basemap,logplot=logplot,ax=ax,label=label,region=region,vmin=vmin,vmax=vmax,cmap=cmap,title=title)
            if show_coef:
                self.plot_eof_coefficients(k,ax=ax2,show_legend=False)
                ax2.grid()
                ti = ax2.get_yticks(); n=len(ti) / 2
                ax2.set_yticks([ti[0],ti[n],ti[-1]])

        if show_coef:
            return f
        else:
            return None



    def _plot_single_EOF(self,k,use_basemap=False,logplot=False,ax=None,label=None,region=None,vmin=None,vmax=None,cmap=None,title=None):
        '''
        plot principal component k

        @param k: number of principal component to plot
        @type k: int

        @param use_basemap: do plot using Basemap
        @type use_basemap: bool

        @param logplot: take log of data for plotting
        @param logplot: bool
        '''
        if k<0:
            raise ValueError, 'k<0'
        #~ d = self.EOF[:,k]; d.shape = self._shape0

        if label == None:
            label = ''
        else:
            label = label + ' '

        if ax == None:
            f = plt.figure()
            ax = f.add_subplot(111)



        #remap data back to original shape
        #1) valid data --> all data
        hlp = np.zeros(len(self._x0mask))*np.nan
        hlp[self._x0mask] = self.EOF[:,k].copy()
        #2) vector --> matrix
        hlp.shape = self._shape0
        hlp = np.ma.array(hlp,mask=np.isnan(hlp))


        #pyCMBS data object
        D = self._x0.copy()
        D.data = hlp
        D.unit = None #reset units as EOF have no physical units
        D.label = label + 'EOF ' + str(k+1).zfill(3) + ' (' + str(round(self._var[k]*100.,2)) + '%)' #caution: labeling is always k+1!
        map_plot(D,use_basemap=use_basemap,logplot=logplot,ax=ax,region=region,vmin=vmin,vmax=vmax,cmap_data=cmap,title=title)

    def reconstruct_data(self,maxn=None,input=None):
        '''
        reconstruct data from EOFs

        @param maxn: specifies the truncation number for EOF reconstruction
        @type maxn: int

        @param input: if this argument is given, then the reconstruction is based on the modes specified
                      in this list. It can be an arbitrary list of mode valid mode indices
        @type input: list of int
        '''

        sh = (self.n,np.prod(self._shape0))
        F = np.zeros(sh)



        #- reconsturction list
        if input == None:
            #use all data up to maxn
            if maxn == None:
                maxn = self.n
            thelist = range(maxn)
        else:
            #use user defined list
            thelist = input


        #- reconstruct data matrix
        for i in thelist:
            #~ a = np.asarray([self.EOF[:,i]]).T
            #remap to original geometry first
            hlp = np.zeros(len(self._x0mask))*np.nan
            hlp[self._x0mask] = self.EOF[:,i].copy()

            a = np.asarray([hlp]).T
            c = np.asarray([self.eigvec[:,i]])
            F += np.dot(a,c).T

        #- generate data object to be returned
        D = self._x0.copy()
        F.shape = self._x0.data.shape
        D.data = np.ma.array(F,mask=np.isnan(F))

        return D

#-----------------------------------------------------------------------

    def get_correlation_matrix(self):
        '''
        correlation matrix of original data [ntimes,ntimes]
        '''
        return np.corrcoef(self.x,rowvar=0)

    def get_eof_data_correlation(self,plot=True):
        '''
        get correlation between original data and PCs
        '''

        c=np.corrcoef(self.x,self.EOF,rowvar=0) #correate PCS and original data
        c1 = c[self.n:,0:self.n]
        if plot:
            f = plt.figure()
            ax = f.add_subplot(111)

            im=ax.imshow(c1**2,interpolation='nearest',vmin=0,vmax=1.,origin='lower')
            plt.colorbar(im)
            ax.set_title('$R^2$ of EOFs with original data')
            ax.set_xlabel('original data channel #')
            ax.set_ylabel('PC #')

    def plot_channnel_correlations(self,samp):
        '''
        generate a scatterplot of correlations of call channles vs. each other

        @param samp: stepsize for subsampling of data for faster plotting
        @type samp: int
        '''

        f = plt.figure()
        cnt = 1
        for i in range(self.n):
            x = self.x[::samp,i]
            for j in range(self.n):
                print i,j
                y = self.x[::samp,j]
                if j >=i:
                    ax = f.add_subplot(self.n,self.n,cnt)
                    #~ ax.set_aspect('equal')
                    ax.hexbin(x,y,mincnt=1,bins='log')
                    ax.set_ylim(ax.get_xlim())
                    ax.set_xticks([])
                    ax.set_yticks([])
                cnt += 1

        f.subplots_adjust(wspace=0.,hspace=0.,right=1.,left=0.,bottom=0.,top=1.)

        return f







#-----------------------------------------------------------------------

class ANOVA():
    '''
    main class to perform an ANOVA analysis
    using C{Data} objects

    treatments are derived from data timeseries
    blocks are obtained by different experiments
    '''
    def __init__(self):
        self.experiments = []
        self.data = {}

    def add_experiment(self,label):
        self.experiments.append(self.__trim(label))
        self.data.update({self.__trim(label):[]})

    def __trim(self,s):
        return s.replace(' ','_')

    def add_data(self,e,d):
        '''
        adds data for each experiment to a list
        '''
        k = self.__trim(e) #experiment key
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
                #~ raise ValueError, 'All data objects need to me masked arrays!' #this is needed,as otherwise common mask generation is not possible
            self.data[k].append(d)
        else:
            raise ValueError, 'Experiment was not yet registered!'

    def analysis(self,analysis_type = None):
        '''
        perform ANOVA analysis based on the data given
        '''

        if analysis_type == None:
            raise ValueError, 'Needs to specify type of ANOVA to be performed'

        #1) check if all data has the same length
        # this is not a strong prerequesite for ANOVA as such, but for
        # the implementation here
        # generate a common mask with points that are valid throughout
        # all timesteps and are not masked
        self._get_valid_mask() #--> self.mask

        #2) rearange data to work only with valid (not masked data)
        #   check also that only data which is valid in all cases is used
        #   it is realized by getting the valid indices of the data
        idx = np.argwhere(self.mask) #returns list of indices of valid pixels

        #3) perform ANOVA analysis on all pixels individually
        resA = np.zeros(self.refshape)*np.nan
        resB = np.zeros(self.refshape)*np.nan
        resI = np.zeros(self.refshape)*np.nan
        resE = np.zeros(self.refshape)*np.nan
        resPA = np.zeros(self.refshape)*np.nan
        resPB = np.zeros(self.refshape)*np.nan
        resPI = np.zeros(self.refshape)*np.nan

        resAa = np.zeros(self.refshape)*np.nan
        resBa = np.zeros(self.refshape)*np.nan
        resIa = np.zeros(self.refshape)*np.nan

        for p in idx:
            if analysis_type == 'one': #one way ANOVA
                m = self._data2anova1(p) #[nrens,ntime]
                A = Anova1(m)
                A.one_way_anova(verbose=False)

                resA[p[0],p[1]] = A.get_fractional_variance_explained(adjust=False) #ssa
                resPA[p[0],p[1]] = A.p
                resB[p[0],p[1]] = A.sse / A.sst #todo: adjustment here ???



            elif analysis_type == 'two': #two way ANOVA
                m = self._data2anova2(p)
                #- perform 2-way anova
                A = Anova2(m)
                A.two_way_anova_with_replication()

                resA[p[0],p[1]] = A.get_fractional_variance_explained('a',adjust=False) #todo: adjust variance
                resB[p[0],p[1]] = A.get_fractional_variance_explained('b',adjust=False) #todo: adjust variance
                resI[p[0],p[1]] = A.get_fractional_variance_explained('i',adjust=False) #todo: adjust variance
                resE[p[0],p[1]] = A.get_fractional_variance_explained('e',adjust=False) #todo: adjust variance

                resPA[p[0],p[1]] = A.p_ssa
                resPB[p[0],p[1]] = A.p_ssb
                resPI[p[0],p[1]] = A.p_ssi

                resAa[p[0],p[1]] = A.get_fractional_variance_explained('a',adjust=True) #todo: adjust variance
                resBa[p[0],p[1]] = A.get_fractional_variance_explained('b',adjust=True) #todo: adjust variance
                resIa[p[0],p[1]] = A.get_fractional_variance_explained('i',adjust=True) #todo: adjust variance

                #todo significance


            else:
                raise ValueError, 'Invalid ANOVA type'

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


    def _data2anova1(self,p):
        '''
        extract from the database all the data
        relevant for a single location, given by the indices in p

        p = indices

        ---> time (nt)
        |
        |
        v experiment (nexp)
        '''

        nexp = len(self.experiments)

        if nexp != 1:
            raise ValueError, 'one-way anova only valid for signle experiments!'

        x = np.zeros((self.n,self.nt))*np.nan #[nrens,nt]

        e = self.experiments[0]
        for i in range(self.n):
            d = self.data[e][i] #data object for experiment 'e' and ensemble nr [i]
            x [i,:] = d.data[:,p[0],p[1]]

        if np.any(np.isnan(x) > 0):
            raise ValueError, 'Something is wrong: not all data valid!'

        return x



    def _data2anova2(self,p):
        '''
        extract from the database all the data
        relevant for a single location, given by the indices in p

        p = indices

        ---> time (nt)
        |
        |
        v experiment (nexp)

        '''
        nexp = len(self.experiments)
        x = np.zeros((nexp,self.nt,self.n))*np.nan

        for j in range(len(self.experiments)):
            e = self.experiments[j]
            for i in range(self.n):
                d = self.data[e][i] #data object for experiment 'e' and ensemble nr [i]
                x [j,:,i] = d.data[:,p[0],p[1]]

        if np.any(np.isnan(x) > 0):
            raise ValueError, 'Something is wrong: not all data valid!'

        return x

    def _get_valid_mask(self):
        '''
        generate a mask where all datasets are valid
        '''
        #- check if geometry in general o.k.
        self.__same_geometry()

    def __same_geometry(self):
        '''
        check if all data has the same geometry
        '''
        cnt = 0
        for k in self.experiments:
            d = self.data[k]
            if cnt == 0:
                refshape = d[0].data.shape #first dataset geometry as reference
                self.refshape=(refshape[1],refshape[2])
                nrens = len(d)
                self.n = nrens
                self.nt = refshape[0]
                refmsk = np.ones((refshape[1],refshape[2])).astype('bool')
                cnt += 1

            if len(d) != nrens:
                raise ValueError, 'Invalid Number of ensemble members found!'
            for i in range(len(d)):

                if d[i].data.shape != refshape:
                    print k, i
                    raise ValueError, 'Invalid data shape found!'

                #- generate common mask
                msk = d[i].get_valid_mask()
                #~ print sum(msk), k, i
                refmsk = refmsk & msk

        self.mask = refmsk

        return True

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class SVD():
    '''
    class to perform singular value decomposition analysis
    also known as Maximum covariance analysis (MCA)

    REFERENCES
    ==========
    1. Bjoernsson and Venegas: A manual for EOF and SVD analyses of Climate Data,
                           McGill University, available online
    '''

    def __init__(self,X,Y,scf_threshold = 0.01,label=''):
        '''
        constructor for SVD class

        @param X: x-variable field
        @type X: C{Data} object

        @param Y: y-variable field
        @type Y: C{Data} object

        @param scf_threshold: threshold for explained variance until which result maps are plotted
        @type scf_threshold: float

        @param label: label for labeling figures
        @type label: str
        '''

        x = X.data.copy(); y = Y.data.copy()

        n = len(x)
        if n != len(y):
            raise ValueError, 'Datasets need to have same timelength!'

        x.shape = (n,-1) #[time,position]
        y.shape = (n,-1)

        self.x = x; self.y = y
        self.X = X; self.Y = Y

        self.time = X.time.copy()

        self.label = label

        self.scf_threshold = scf_threshold #threshold for explained variance until which result maps are plotted
        self.dpi  = 150 #output dpi for plotting
        self.ext = 'pdf' #file extension for plotting
        self.use_basemap = False

#-----------------------------------------------------------------------

    def _get_valid_timeseries(self,x):
        '''
        get only points where all
        timeseries are valid

        @param x: array [time,position]
        @type x: numpy masked array

        @return: masked array and mask

        @todo: handling of invalid, gappy data
        '''
        sx = x.data.sum(axis=0) #calculate sum. If all values are valid, then no a float should be there, else Nan
        mx = ~np.isnan(sx) #masked array
        m1 = mx
        #~ m1[x.mask]=False #apply also mask of masked array

        r = np.ones( (len(x.data) , sum(mx)  ))*np.nan
        for i in range(len(x.data)):
            tmp = x.data[i,:]
            if np.any(np.isnan(tmp[m1])):
                print 'nans are not allowed here!'
                print tmp
                stop

            r[i,:] = tmp[m1]*1.
            del tmp



        return np.ma.array(r),m1

#-----------------------------------------------------------------------

    def __detrend_time(self,x):
        '''
        given a variable x[time,position]
        the data is linear detrended individually for each position

        @param x: data array [time,position]
        @type x: numpy array
        @return: return detrended array [time,position]
        '''

        if x.ndim !=2:
            raise ValueError, 'Invalid shape for detrending'
        n,m = x.shape
        for i in range(m):
            h = x[:,i].copy(); h = plt.detrend_linear(h)
            x[:,i] = h.copy()
        return x

#-----------------------------------------------------------------------

    def __time_normalization(self,x):
        '''
        normalize timeseries x [time,position]
        by dividiing by the standard deviation

        @param x: data array [time,position]
        @type x: numpy array

        @return: normalized timeseries numpy array
        '''
        nt,nx = np.shape(x)
        s = x.std(axis=0) #temporal standard deviation
        S = np.repeat(s,nt).reshape(nx,nt).T #generate array with all same std
        x = x / S
        return x

#-----------------------------------------------------------------------

    def svd_analysis(self,detrend=True,varnorm=False):
        '''
        perform SVD analysis

        @param detrend: detrend data
        @type detrend: bool

        @param varnorm: normalize variance of time series
        @type varnorm: bool
        '''

        #/// perform SVN only for data points which are valid throughout entire time series ///
        x,mskx = self._get_valid_timeseries(self.x)
        y,msky = self._get_valid_timeseries(self.y)
        self.mskx = mskx; self.msky = msky

        #/// detrend the data for each grid point ///
        print 'Detrending ...'
        if detrend:
            x = self.__detrend_time(x)
            y = self.__detrend_time(y)
        print 'Detrended!'

        #/// normalized each timeseries by its variance
        if varnorm:
            x = self.__time_normalization(x)
            y = self.__time_normalization(y)

        #/// calculate covariance matrix
        print 'Construct covariance matrix ...'
        C = dot(x.T,y) #this covariance matrix does NOT contain the variances of the individual grid points, but only the covariance terms!
        print 'Done!'
        self.C = C
        self.x_used = x.copy() #store vectors like they are used for SVD calculations
        self.y_used = y.copy()

        #/// singular value decomposition
        print '   Doing singular value decomposition ...'
        U, s, V = linalg.svd( C )
        print 'Done!'
        L = linalg.diagsvd(s, len(C), len(V) ) #construct diagonal maxtrix such that U L V.T = C; this is somewhat python specific

        #/// expansion coefficients (time series)
        A = dot(x,U); B = dot(y,V)

        #/// store results
        self.U = U; self.V = V
        self.L = L; self.A = A; self.B = B
        self.scf = (s*s) / sum(s*s) #fractions of variance explained CAUTION: not properly described in manual if squared or not!
        self.__get_mode_correlation() #calculate correlation between modes

#-----------------------------------------------------------------------

    def __get_mode_correlation(self):
        '''
        calculate correlations between expansion modes
        of the two fields
        '''
        self.mcorr = []
        for i in range(len(self.scf)):
            c = np.corrcoef(self.A[:,i],self.B[:,i])[0][1]
            self.mcorr.append(c)
        self.mcorr = np.asarray(self.mcorr)

#-----------------------------------------------------------------------

    def get_singular_vectors(self,mode):
        '''
        return the singular vectors of both fields for a specific mode
        as a spatial (2D) field

        @param mode: mode to be shown
        @type mode: int

        @return: numpy arrays for U and V
        '''

        #x_used is a vector that only contains the valid values that were used for caluclation of covariance matrix C
        #mskx is the corresponding mask that maps x_used to the original geometry (both are estimated with _get_valid_timeseries()  )
        u = self.U[:,mode]; v = self.V[:,mode] #get singular vectors

        #map singular vectors to 2D
        udat = self._map_valid2org(u,self.mskx,self.X.data[0,:,:].shape)
        vdat = self._map_valid2org(v,self.msky,self.Y.data[0,:,:].shape)

        U = self.X.copy(); U.label = 'U(' + str(mode) + ')'
        U.data = udat.copy()

        V = self.Y.copy(); V.label = 'V(' + str(mode) + ')'
        V.data = vdat.copy()

        return U,V

#-----------------------------------------------------------------------

    def _map_valid2org(self,data,mask,target_shape):
        '''
        map valid data vector back to
        original data shape

        @param data: data vector that was used for SVD calculations (1D)
        @type data: numpy array vector

        @param mask: 1D data mask
        @type mask: numpy data (1D)

        @param target_shape: shape to map to
        @type target_shape: geometry tuple
        '''
        sz = np.shape(data)
        if sz[0] != mask.sum():
            print sz[1]; print mask.sum()
            raise ValueError, 'Inconsistent mask and data'

        res = np.ones(target_shape)*np.nan; res.shape = (-1)
        res[mask] = data.copy()

        res = np.ma.array(res,mask=np.isnan(res))
        res = np.reshape(res,target_shape)

        return res

#-----------------------------------------------------------------------

    def plot_var(self,ax=None,filename=None,maxvar=1.):
        '''
        plot explained variance

        @param ax: axis to put plot in. If None, then a new figure is generated
        @type ax: matplotlib axis

        @param filename: name of the file to store figure to (if None, then no file is saved)
        @type filename: str

        @param maxvar: upper limit of variance plot
        @type maxvar: float
        '''

        def make_patch_spines_invisible(ax):
            #http://matplotlib.sourceforge.net/examples/pylab_examples/multiple_yaxis_with_spines.html
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.itervalues():
                sp.set_visible(False)

        if ax == None:
            fig=plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = ax
            fig = ax.figure

        fig.subplots_adjust(right=0.75)

        ax1 = ax.twinx() #axis for cumulated variance
        ax2 = ax.twinx()

        ax2.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(ax2)
        ax2.spines["right"].set_visible(True)

        n = len(self.scf)
        ax.step(np.arange(n),self.scf,where='post')
        ax.set_ylabel('fraction of variance explained',color='blue')
        ax.set_xlabel('mode')
        ax.set_ylim(0.,maxvar)
        ax.grid()

        ax1.plot(np.cumsum(self.scf),color='red')
        ax1.set_ylabel('cumulated variance [-]',color='red')
        ax1.set_ylim(0.,1.)

        ax2.plot(np.arange(n),self.mcorr,color='green')
        ax2.set_ylabel('mode correlation [-]',color='green')
        ax2.set_ylim(-1,1.)
        ax2.grid(color='green')

        ax .tick_params(axis='y', colors='blue')
        ax1.tick_params(axis='y', colors='red')
        ax2.tick_params(axis='y', colors='green')

        if filename != None:
            oname = filename + '_mode_var.' + self.ext
            ax.figure.savefig(oname,dpi=self.dpi)

#-----------------------------------------------------------------------

    def plot_singular_vectors(self,mode,use_basemap=False,logplot=False,filename=None):
        '''
        generate maps of singular vectors U and V

        mode (list) : list of modes to be plotted
        '''

        #--- mode list
        if mode == None: #plot all modes with variance contained
            mode_list = []
            for i in range(len(self.scf)):
                if self.scf[i] > self.scf_threshold:
                    mode_list.append(i)
        else:
            mode_list = [mode]

        #--- generate plots
        for i in mode_list:
            fig=plt.figure()
            ax1 = fig.add_subplot(121); ax2 = fig.add_subplot(122)

            U,V = self.get_singular_vectors(i) #get singular vector fields

            #--- determine min/max values
            mu = U.data.mean(); su = U.data.std()
            mv = V.data.mean(); sv = V.data.std()

            su = 1.96*su; sv = 1.96*sv #not used at the moment

            if logplot:
                umin=None;umax=None
                vmin=None;vmax=None
            else:
                umin=mu-su;umax=mu+su
                vmin=mv-sv;vmax=mv+sv

            map_plot(U,use_basemap=use_basemap,ax=ax1,logplot=logplot,vmin=umin,vmax=umax)
            map_plot(V,use_basemap=use_basemap,ax=ax2,logplot=logplot,vmin=vmin,vmax=vmax)

            fig.suptitle('Mode: #' + str(i) + ' scf: ' + str(round(self.scf[i]*100.,1) ) + '%', size=14)

            if filename != None:
                fig.savefig(filename + '_singular_vectors_mode_' + str(i).zfill(5) + '.pdf' )

#-----------------------------------------------------------------------

    def plot_correlation_map(self,mode,ax1in=None,ax2in=None,pthres=1.01,plot_var=False,filename=None,region1=None,region2=None,regions_to_plot=None):
        '''
        plot correlation map of an SVN mode
        with original data

        mode specifies the number of the mode that should be correlated

        ctype specifies if homogeneous (homo) or heterogeneous (hetero)
        correlations shall be calculated. Homogeneous, means correlation
        of expansion coefficients with the same geophysical field, while
        heterogeneous means correlation with the other geophysical field

        pthres specifies the significance level. Values > pthres will be masked
        if you want to plot e.g only points with significant correlation at p < 0.05, then
        set pthres = 0.05

        plot_var: plot variance instead of correlation
        '''


        n1,m1 = self.A.shape
        n2,m2 = self.B.shape

        if mode != None:
            if mode > m1-1:
                raise ValueError, 'Mode > A'
            if mode > m2-1:
                raise ValueError, 'Mode > B'

        if mode == None: #plot all modes with variance contained
            mode_list = []
            for i in range(len(self.scf)):
                if self.scf[i] > self.scf_threshold:
                    mode_list.append(i)
        else:
            mode_list = [mode]

        def plot_cmap(R,ax,title,vmin=-1.,vmax=1.,plot_var=False,use_basemap=False,region=None,cmap='RdBu_r',cticks=None,regions_to_plot=None):
            '''
            R data object
            '''

            if plot_var:
                O = R.copy()
                O.data = O.data*O.data
                O.label = 'exp.frac.var.'
            else:
                O = R.copy()
                O.label = 'correlation'

            #calculate mean and stdv
            O.label = O.label + ' (' + str(round(O.data.mean(),2)) + ' ' + str(round(O.data.std(),2)) + ')'

            map_plot(O,use_basemap=use_basemap,ax=ax,region=region,cmap_data=cmap,vmin=vmin,vmax=vmax,cticks=cticks,title=title,regions_to_plot=regions_to_plot,show_stat=True)








        #/// calculate correlations and do plotting



        for i in mode_list:
            fig=plt.figure(figsize=(6,8))
            ax1a = fig.add_subplot(421) #homogeneous plots
            ax1b = fig.add_subplot(423)
            ax1c = fig.add_subplot(425)
            #~ ax1d = fig.add_subplot(524)

            ax2a = fig.add_subplot(422) #heterogeneous plots
            ax2b = fig.add_subplot(424)
            ax2c = fig.add_subplot(426)
            #~ ax2d = fig.add_subplot(528)

            ax3  = fig.add_subplot(515) #expansion coefficients




            #homogeneous correlations
            Rout1_ho,Sout1_ho,Iout1_ho,Pout1_ho,Cout1_ho = self.X.corr_single(self.A[:,i],pthres=pthres)
            Rout2_ho,Sout2_ho,Iout2_ho,Pout2_ho,Cout2_ho = self.Y.corr_single(self.B[:,i],pthres=pthres)

            #heterogeneous correlations
            Rout1_he,Sout1_he,Iout1_he,Pout1_he,Cout1_he = self.X.corr_single(self.B[:,i],pthres=pthres)
            Rout2_he,Sout2_he,Iout2_he,Pout2_he,Cout2_he = self.Y.corr_single(self.A[:,i],pthres=pthres)

            #--- plot maps
            print 'Starting map plotting'
            #homogeneous
            plot_cmap(Rout1_ho,ax1a,'correlation (homo)',plot_var=False,use_basemap=self.use_basemap,region=region1,vmin=-0.8,vmax=0.8,cmap='RdBu_r',cticks=[-1.,-0.5,0.,0.5,1.],regions_to_plot=regions_to_plot) #correlation field 1
            plot_cmap(Rout2_ho,ax1b,'correlation (homo)',plot_var=False,use_basemap=self.use_basemap,region=region2,vmin=-0.8,vmax=0.8,cmap='RdBu_r',cticks=[-1.,-0.5,0.,0.5,1.],regions_to_plot=regions_to_plot) #correlation field 2
            plot_cmap(Rout2_ho,ax1c,'exp.frac.var (homo)'   ,plot_var=True,use_basemap=self.use_basemap,region=region2,vmin=0.,vmax=0.6,cmap='YlOrRd',cticks=[0.,0.25,0.5],regions_to_plot=regions_to_plot)  #explained variance field 2
            #~ plot_cmap(Cout2_ho,ax1d,'covariance (homo)',plot_var=False,use_basemap=self.use_basemap,region=region2,cmap='jet',regions_to_plot=regions_to_plot,vmin=None,vmax=None) #explained covariance

            #heterogeneous
            plot_cmap(Rout1_he,ax2a,'correlation (hetero)',plot_var=False,use_basemap=self.use_basemap,region=region1,vmin=-0.8,vmax=0.8,cmap='RdBu_r',cticks=[-1.,-0.5,0.,0.5,1.],regions_to_plot=regions_to_plot) #correlation field 1
            plot_cmap(Rout2_he,ax2b,'correlation (hetero)',plot_var=False,use_basemap=self.use_basemap,region=region2,vmin=-0.8,vmax=0.8,cmap='RdBu_r',cticks=[-1.,-0.5,0.,0.5,1.],regions_to_plot=regions_to_plot) #correlation field 2
            plot_cmap(Rout2_he,ax2c,'exp.frac.var (hetero)'   ,plot_var=True,use_basemap=self.use_basemap,region=region2,vmin=0.,vmax=0.6,cmap='YlOrRd',cticks=[0.,0.25,0.5],regions_to_plot=regions_to_plot)  #explained variance field 2
            #~ plot_cmap(Cout2_he,ax2d,'covariance (hetero)',plot_var=False,use_basemap=self.use_basemap,region=region2,cmap='jet',regions_to_plot=regions_to_plot,vmin=None,vmax=None) #explained covariance

            #expansion coefficients
            self.plot_expansion_correlation(i,ax=ax3)

            #figure title
            fig.suptitle(self.label + ': Mode: #' + str(i) + ' (scf: ' + str(round(self.scf[i],2)) + ')',size=14)
            fig.subplots_adjust(wspace=0.0,hspace=0.5)
            fig.set_figheight(10.)

            #--- save figure
            if filename != None:
                oname = filename + '_mode_' + str(i) + '.' + self.ext
                ax1a.figure.savefig(oname,dpi=self.dpi)

#-----------------------------------------------------------------------

    def _get_variance_field(self,X,E,mode,pthres=1.01):
        '''
        calculate variance field for a particular mode
        (explained variance by a particular expansion mode)
        This is obtained by correlating an expansion mode to
        a particular data field

        @param X: data field that should be explained by expansion coefficient
        @type X: C{Data} object

        @param E: expansion cofficient to be used for correlation calculation
        @type E: numpy array

        @return: squared correlation as C{Data} object
        '''
        Rout,Sout,Iout,Pout = X.corr_single(E[:,mode],pthres=pthres)
        Rout.data = Rout.data*Rout.data
        return Rout #return squared correlation to get variance

#-----------------------------------------------------------------------

    def reconstruct_variance_fraction(self,X,E,mode_list,pthres=1.01):
        '''
        reconstruct variance of data based on
        a list of modes that should be used
        for that purpose.

        The variances of the different modes are added
        assuming that they are indpendent of each other

        mode_list : list with up to N modes

        X Data object
        E expansion coefficients (data object)

        returns:
        array with variance
        '''

        O = None
        for mode in mode_list:
            V = self._get_variance_field(X,E,mode,pthres=pthres)
            if O == None:
                O = V.data
            else:
                O = O + V.data #add variances

        O = np.ma.array(O,mask=np.isnan(O))

        return O

#-----------------------------------------------------------------------

    def print_mode_statistic(self,filename=None):
        '''
        print statistic of modes

        @param filename: filename to save table to
        @type filename: str
        '''
        sep = ' & '; rnd = 2
        self.__get_mode_correlation() #calculate mode correlations

        if filename != None:
            o = open(filename,'w')
            o.write('mode' + sep + 'scf' + sep + 'r' + ' \\\ ' + '\n')

        for i in np.arange(len(self.scf)):
            if self.scf[i] > self.scf_threshold:
                print i, self.scf[i], self.mcorr[i]
                if filename != None:
                    s = str(i) + sep + str(np.round(self.scf[i],rnd)) + sep +  str(np.round(self.mcorr[i],rnd)) + ' \\\ ' +  '\n'
                    o.write(s)

        if filename != None:
            o.close()

#-----------------------------------------------------------------------

    def plot_expansion_correlation(self,mode,ax=None):
        '''
        plot correlation and time series of expansion coeffcients

        @param mode: mode to plot
        @type mode: int

        @param ax: axis to plot data
        @type ax: matplotlib axis
        '''
        if ax == None:
            fig=plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = ax

        ax.plot(plt.num2date(self.time),self.A[:,mode]/np.std(self.A[:,mode]),label='A',color='red')
        ax.plot(plt.num2date(self.time),self.B[:,mode]/np.std(self.B[:,mode]),label='B',color='blue',linestyle='--')
        c = np.corrcoef(self.A[:,mode],self.B[:,mode])[0][1]
        plt.legend()
        ax.set_title('normalized expansion coefficient #' + str(mode) + ' (r=' + str(round(c,2)) + ')',size=10)
        ax.set_xlabel('time'); ax.set_ylim(-3.,3.)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

class Diagnostic():
    def __init__(self,x,y=None):
        '''
        constructor for diagnostic class
        diagnostic for one or multiple data sets

        @param x: x data to be analyzed
        @type x: C{Data} object

        @param y: y data to be analyzed
        @type y: C{Data} object
        '''
        self.x = x
        if y != None:
            self.y = y

#-----------------------------------------------------------------------

    def get_n(self):
        '''
        return the number of valid samples
        '''
        xm = ~self.xvec.mask #vector with valid sample
        ym = ~self.yvec.mask
        m = xm & ym
        return sum(m)

#-----------------------------------------------------------------------

    def get_rmse_value(self):
        '''
        calculate root-mean-squared error

        @param self: Diagnostic object
        @type self : Diagnostic object
        '''
        return np.sqrt(np.mean((self.xvec - self.yvec)**2))

#-----------------------------------------------------------------------

    def lagged_correlation_vec(self,lags,pthres=1.01,detrend_linear=False,detrend_mean=False):
        '''
        lagged correlation for two vectors

        x,y Data objects, where data needs to have been pre-processed
        to be a single vector

        lags: list of lags
        '''

        if self.x.data.shape != self.y.data.shape:
            raise ValueError, 'Invalid geometries!'

        if not plt.isvector(self.x.data):
            raise ValueError, 'Routine works only with vectorized data!'

        if any(lags < 0.):
            raise ValueError, 'Negative lags currently not supported yet!'

        CO = []
        for lag in lags:
            hlpx = self.x.data.copy(); hlpy = self.y.data.copy()

            if detrend_linear:
                hlpx = plt.detrend_linear(hlpx);
                hlpy = plt.detrend_linear(hlpy);
            if detrend_mean:
                hlpx = plt.detrend_mean(hlpx);
                hlpy = plt.detrend_mean(hlpy);

            #1) temporal subsetting of data
            if lag > 0:
                #temporal subset data
                hlpx = hlpx[lag:]; hlpy = hlpy[:-lag]

            #2) calculation of lagged correlation
            if len(hlpx) > 1:
                slope, intercept, r_value, p_value, std_err = sci.stats.linregress(hlpx,hlpy)
            else:
                r_value = np.nan; p_value =2.

            if p_value < pthres:
                CO.append(r_value)
            else:
                CO.append(np.nan)

        #~ #-- find best lag
        #~ best_lag=abs(np.asarray(CO)).argmax(axis=0).astype('float')
        #~ print 'BEST LAG: ', best_lag

        CO = np.asarray(CO)

        return CO

#-----------------------------------------------------------------------

    def get_correlation_value(self):
        '''
        get correlation between two vectors
        '''
        c = np.ma.corrcoef(self.xvec,self.yvec)[0][1]
        return c #correlation coefficient

#-----------------------------------------------------------------------

    def _mat2vec(self,mask = None):
        '''
        concatenate all information into
        a vector and apply a given
        mask if desired

        @param mask: mask to be applied
        @type mask: numpy array
        '''

        #--- generated copies and mask data if desired
        X = self.x.copy()
        if mask != None:
            X._apply_mask(mask,keep_mask=False)
        if self.y != None:
            Y = self.y.copy()
            if mask != None:
                Y._apply_mask(mask,keep_mask=False)
        else:
            Y = None

        #--- vectorize the data (concatenate in space and time)
        xvec = X.data.copy()
        xvec.shape = (-1)
        if self.y != None:
            yvec = Y.data.copy(); yvec.shape = (-1)

        self.xvec = xvec; self.yvec = yvec

#-----------------------------------------------------------------------

    def calc_reichler_index(self,weights):
        '''
        calculate index after Reichler & Kim (2008)
        for a single model

        it is assumed that the field has a time component

        variable x is assumed to be the reference dataset


        returns E**2
        '''

        if not hasattr(self,'y'):
            raise ValueError, 'Can not calculate Reichler & Kim index without a second variable!'

        weights = weights.copy()
        x = self.x.data.copy(); y = self.y.data.copy()
        std_x = self.x.std.copy()

        if np.shape(x) != np.shape(y):
            print np.shape(x), np.shape(y)
            raise ValueError, 'Invalid shapes of arrays!'

        if x.ndim == 1: #only timeseries
            e2 = sum(weights * (x-y)**2 / std_x)
        else:
            n = len(x)
            x.shape = (n,-1) #[time,index]
            y.shape = (n,-1)
            std_x.shape = (n,-1); weights.shape = (-1)
            if np.shape(x[0,:]) != np.shape(weights):
                raise ValueError, 'Invalid shape for weights!'
            e2 = []
            for i in range(n):
                d = weights * ( (x[i,:]-y[i,:])**2)   / std_x[i,:]
                e2.append(np.sum(d)) #sum at end to avoid nan's   #it is important to use np.sum() !!
            e2 = np.asarray(e2)

        if np.any(np.isnan(e2)):
            raise ValueError, 'Reichler: e2 contains NAN, this happens most likely if STDV == 0'
        return e2

#-----------------------------------------------------------------------

    def get_correlationxxxxx(self,lag=0,nlags=None):
        '''
        calculate correlation between two data sets
        '''

        print 'Calculating correlation with lag=', lag, nlags

        def NormCrossCorrSlow(x1, x2,
                  nlags=400):
            res=[]
            lags=[]
            for i in range(-(nlags/2),nlags/2,1):
                lags.append(i)
                if i<0:
                    xx1=x1[:i]
                    xx2=x2[-i:]
                elif i==0:
                    xx1=x1
                    xx2=x2
                else:
                    xx1=x1[i:]
                    xx2=x2[:-i]

                xx1 = xx1 - xx1.mean()
                xx2 = xx2 - xx2.mean()

                res.append( (xx1*xx2).sum() /( (xx1**2).sum() *(xx2**2).sum() )**0.5)

            return np.array(res), np.array(lags)




        if not hasattr(self,'y'):
            raise ValueError, 'No y variable existing!'

        x = self.x.data.copy(); y = self.y.data.copy()

        s1 = np.shape(x); s2 = np.shape(y)
        if s1 != s2:
            print s1,s2
            raise ValueError, 'Invalid shapes!'


        n = np.shape(x)[0]
        #~ print n
        x.shape = (n,-1); y.shape = (n,-1)




        #~ print s1
        print x.shape
        print s1
        print x.ndim
        print len(s1)
        if len(s1) ==2:
            ndata = s1[1]*s1[2]
        elif len(s1) == 1:
            ndata = 1 #vector
        else:
            raise ValueError, 'Invalid shape!'

        R=np.zeros(ndata)*np.nan

        for i in range(ndata):
            xx=x[:,i]
            yy=y[:,i]
            msk = (~np.isnan(xx) & ~np.isnan(yy))
            if sum(msk) > 3:

                if lag == 0:
                    slope, intercept, r_value, p_value, std_err = sci.stats.linregress(xx[msk],yy[msk])
                else:

                    print nlags
                    if nlags == None:
                        nlags = len(xx[msk])
                        if np.mod(nlags,2) == 0:
                            nlags = nlags + 1
                    print nlags

                    r1,lags = NormCrossCorrSlow(xx[msk],yy[msk],nlags=nlags)
                    idx = nlags/2+lag #todo something does not work with the indices !!!
                    print idx, nlags, len(r1)
                    r_value = r1[idx]
                    #~ print r_value, r1
                    #~ print lags[nlags/2 + lag], r1[nlags/2+lag]
                    #~ pl.close('all')
                    #~ pl.plot(lags,r1)
                    #~
                R[i] = r_value

        #~ print np.shape(self.x.data[0,:])
        #~ print np.shape(R)
        print x.ndim

        print self.x.data.ndim
        if len(s1) == 2:
            R = np.reshape(R,np.shape(self.x.data[0,:]))
        else:
            R = R

        return R

#-----------------------------------------------------------------------

    def slice_corr(self,timmean=True,spearman=False,partial=False,z=None):
        '''
        perform correlation analysis for
        different starting times and length
        of the correlation period

        if timmean=True then the correlation is caluclated
        on basis of the mean spatial fields, thus
        the timeseries of each pixels is averaged over time
        before the correlation calculation

        partial: do partial correlation
        z: condition in case of partial correlation
        '''

        if partial:
            if spearman:
                raise ValueError, 'Spearman and partial correlation not supported'
            if z == None:
                raise ValueError, 'no z-value given for partial correlation!'
            if self.x.data.shape != z.data.shape:
                print self.x.data.shape
                print z.data.shape
                raise ValueError, 'Invalid geometries for partial correlation!'

        x=self.x.data.copy()

        if not hasattr(self,'y'):
            #if no y value is given, then time is used as independent variable
            print 'No y-value specified. Use time as indpendent variable!'

            y = x.copy()
            x = np.ma.array(self.x.time.copy(),mask = self.x.time < 0. )

        else:
            y = self.y.data.copy()

        if np.shape(x) != np.shape(y):
            raise ValueError, 'slice_corr: shapes not matching!'



        #--- reshape data
        n = len(x) #timesteps
        x.shape = (n,-1) #size [time,ngridcells]
        y.shape = (n,-1)

        if partial:
            z = z.data.copy()
            z.shape = (n,-1)




        R=np.ones((n,n))*np.nan; P=np.ones((n,n))*np.nan
        L=np.ones((n,n))*np.nan; S=np.ones((n,n))*np.nan

        #--- perform correlation analysis
        print '   Doing slice correlation analysis ...'
        i1 = 0
        while i1 < n-1: #loop over starting year
            i2 = i1 + 2
            #- loop over different lengths
            while i2 < len(x)-1:
                length = i2-i1

                if timmean:
                    ''' temporal mean -> all grid cells only (temporal mean) '''
                    #print 'drin'
                    xdata = x[i1:i2,:].mean(axis=0)
                    ydata = y[i1:i2,:].mean(axis=0)

                    xmsk  = xdata.mask; ymsk = ydata.mask
                    msk = xmsk | ymsk

                    if partial:
                        raise ValueError, 'No timmean supported yet for partial correlation!'

                else:
                    ''' all grid cells at all times '''
                    xdata = x.data[i1:i2,:]; ydata = y.data[i1:i2,:]
                    xmsk  = x.mask[i1:i2,:]
                    ymsk  = y.mask[i1:i2,:]
                    msk   = xmsk | ymsk

                    if partial:
                        zdata = z.data[i1:i2,:]
                        zmsk  = z.mask[i1:i2,:]
                        msk = msk | zmsk
                        zdata = zdata[~msk].flatten()

                xdata = xdata[~msk].flatten()
                ydata = ydata[~msk].flatten()

                #use spearman correlation
                if spearman:
                    tmpx = xdata.argsort()
                    tmpy = ydata.argsort()
                    xdata = tmpx
                    ydata = tmpy

                if partial:
                    #calculate residuals for individual correlations
                    slope, intercept, r, p, stderr = sci.stats.linregress(zdata,xdata)
                    xdata = (xdata - intercept) / slope

                    slope, intercept, r, p, stderr = sci.stats.linregress(zdata,ydata)
                    ydata = (ydata - intercept) / slope



                slope, intercept, r, p, stderr = sci.stats.linregress(xdata,ydata)
                R[length,i1] = r
                P[length,i1] = p
                L[length,i1] = length
                S[length,i1] = slope

                i2 += 1
            i1 += 1

        self.slice_r = R
        self.slice_p = P
        self.slice_length = L
        self.slice_slope = S

#-----------------------------------------------------------------------

    def slice_corr_gap(self,timmean=True,spearman=False,pthres=None):
        '''
        perform correlation analysis for
        different starting times and gap sizes

        if timmean=True then the correlation is caluclated
        on basis of the mean spatial fields, thus
        the timeseries of each pixels is averaged over time
        before the correlation calculation

        '''

        x=self.x.data.copy()

        if not hasattr(self,'y'):
            #if no y value is given, then time is used as independent variable
            print 'No y-value specified. Use time as indpendent variable!'

            y = x.copy()
            x = np.ma.array(self.x.time.copy(),mask = self.x.time < 0. )

        else:
            y = self.y.data.copy()

        if np.shape(x) != np.shape(y):
            raise ValueError, 'slice_corr: shapes not matching!'

        #--- reshape data
        n = len(x) #timesteps

        gaps = np.arange(n)

        x.shape = (n,-1) #size [time,ngridcells]
        y.shape = (n,-1)
        maxgap = n


        #~ print 'maxgap: ', maxgap


        R=np.ones((maxgap,n))*np.nan; P=np.ones((maxgap,n))*np.nan
        L=np.ones((maxgap,n))*np.nan; S=np.ones((maxgap,n))*np.nan

        #--- perform correlation analysis
        print '   Doing slice correlation analysis ...'
        i1 = 0
        while i1 < n-1: #loop over starting year
            i2 = n #always entire time period
            #- loop over different lengths
            for gap in gaps:

                if gap >=  i2 - i1:
                    continue

                if timmean:
                    ''' temporal mean -> all grid cells only (temporal mean) '''

                    raise ValueError, 'TIMMEAN not supported yet for gap analysis'

                    #print 'drin'
                    xdata = x[i1:i2,:].mean(axis=0)
                    ydata = y[i1:i2,:].mean(axis=0)

                    xmsk  = xdata.mask; ymsk = ydata.mask
                    msk = xmsk | ymsk

                else:
                    ''' all grid cells at all times '''
                    #~ xdata = x.data[i1:i2,:]; ydata = y.data[i1:i2,:]
                    xdata = x.data.copy(); ydata = y.data.copy()
                    xmsk  = x.mask.copy() #  [i1:i2,:]
                    ymsk  = y.mask.copy() #  [i1:i2,:]

                    #mask data which has gaps and use whole period elsewhere
                    xmsk[i1:i1+gap,:] = True
                    ymsk[i1:i1+gap,:] = True

                    msk   = xmsk | ymsk

                xdata = xdata[~msk].flatten()
                ydata = ydata[~msk].flatten()

                #use spearman correlation
                if spearman:
                    tmpx = xdata.argsort()
                    tmpy = ydata.argsort()
                    xdata = tmpx
                    ydata = tmpy

                slope, intercept, r, p, stderr = sci.stats.linregress(xdata,ydata)
                R[gap,i1] = r
                P[gap,i1] = p
                L[gap,i1] = gap-1
                S[gap,i1] = slope


            i1 += 1


        if pthres != None: #mask all insignificant values
            R = np.ma.array(R,mask=P>pthres)
            S = np.ma.array(S,mask=P>pthres)

        self.slice_r_gap = R
        self.slice_p_gap = P
        self.slice_length_gap = L
        self.slice_slope_gap = S


#-----------------------------------------------------------------------

    def _set_year_ticks(self,years,ax,axis='x',size=10):
        '''
        set ticks of timeline with
        yearly ticks

        @param years: list of years
        @type years: list

        @param ax: axis to handle
        @type ax: matplotlib axis

        @param axis: specify which axis to handle 'x' or 'y'
        @type axis: str
        '''
        ticks = ax.get_xticks()

        #- calculate ticks from year
        oticks=[]
        for t in ticks:
            if t < 0:
                oticks.append('')
            elif t > len(years)-1:
                oticks.append('')
            else:
                oticks.append(years[int(t)])
        #- set ticks of axis
        if   axis == 'x':
            ax.set_xticklabels(oticks,size=size)
        elif axis == 'y':
            ax.set_yticklabels(oticks,size=size)
        else:
            raise ValueError, 'Invalid axis (set_year_ticks)'

#-----------------------------------------------------------------------

    def plot_slice_correlation(self,pthres = 1.01):
        '''
        plot slice correlation results

        @param pthres: significance threshold. All results with p-values
                       below this threshold will be plotted
        @type pthres: float
        '''
        cmap1 = plt.cm.get_cmap('RdBu_r', 10)
        cmap2 = plt.cm.get_cmap('jet', 10)

        if not hasattr(self,'slice_r'):
            raise ValueError, 'perform slice_corr() before plotting!'

        #- get years of data for ticks
        years = self.x._get_years()

        #- generate plots
        fig=plt.figure(figsize=(12,6))
        fig.subplots_adjust(hspace=0.5)
        self.slice_fig = fig
        ax1=fig.add_subplot(221); ax2=fig.add_subplot(222)
        ax3=fig.add_subplot(223); ax4=fig.add_subplot(224)

        r_data      = self.slice_r.copy()
        p_data      = self.slice_p.copy()
        length_data = self.slice_length.copy()
        slope_data  = self.slice_slope.copy()

        msk = p_data > pthres
        r_data[msk]      = np.nan; p_data[msk]      = np.nan
        length_data[msk] = np.nan; slope_data[msk]  = np.nan

        #- correlation
        imr=ax1.imshow(r_data,interpolation='nearest',cmap=cmap1)
        ax1.set_title('correlation')
        plt.colorbar(imr,ax=ax1,shrink=0.8)
        ax1.set_xlabel('start year')
        ax1.set_ylabel('correlation period [years]')

        #- significance
        imp=ax2.imshow(p_data,interpolation='nearest',cmap=cmap2)
        ax2.set_title('p-value')
        plt.colorbar(imp,ax=ax2,shrink=0.8)
        ax2.set_xlabel('start year')
        ax2.set_ylabel('correlation period [years]')

        #- length of period
        iml=ax3.imshow(length_data,interpolation='nearest',cmap='RdBu_r')
        ax3.set_title('length')
        plt.colorbar(iml,ax=ax3,shrink=0.8)
        ax3.set_xlabel('start year')
        ax3.set_ylabel('correlation period [years]')

        #- slope
        ims=ax4.imshow(slope_data,interpolation='nearest',cmap=cmap2)
        ax4.set_title('slope')
        plt.colorbar(ims,ax=ax4,shrink=0.8)
        ax4.set_xlabel('start year')
        ax4.set_ylabel('correlation period [years]')

        #/// set tick labels ///
        self._set_year_ticks(years,ax1,axis='x')
        self._set_year_ticks(years,ax2,axis='x')
        self._set_year_ticks(years,ax3,axis='x')
        self._set_year_ticks(years,ax4,axis='x')

        #- contour plots
        CP1 = ax1.contour(p_data,[0.01,0.05,0.1],linewidths=2)
        CP2 = ax2.contour(p_data,[0.01,0.05,0.1],linewidths=2)
        CP3 = ax3.contour(p_data,[0.01,0.05,0.1],linewidths=2)
        CP4 = ax4.contour(p_data,[0.01,0.05,0.1],linewidths=2)

        ax1.clabel(CP1, inline=1, fontsize=10)
        ax2.clabel(CP2, inline=1, fontsize=10)
        ax3.clabel(CP3, inline=1, fontsize=10)
        ax4.clabel(CP4, inline=1, fontsize=10)
