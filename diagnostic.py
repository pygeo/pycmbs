#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

import numpy as np

import scipy as sci

from matplotlib import pylab as plt

from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes 
import matplotlib as  mpl

from pyCMBS.plots import map_plot


from scipy import linalg, dot;


class SVD():
    '''
    class to perform singular value decomposition analysis
    also known as Maximum covariance analysis (MCA)
    
    REFERENCES
    Björnsson and Venegas: A manual for EOF and SVD analyses of Climate Data, McGill University, available online
    '''
    
    def __init__(self,X,Y,scf_threshold = 0.01,label=''):
        '''
        x,y Data objects
        '''
        
        x = X.data.copy(); y = Y.data.copy()
        
        n = len(x)
        if n != len(y):
            raise ValueError, 'Datasets need to have same timelength!'


        x.shape = (n,-1) #[time,position]
        y.shape = (n,-1) 
        
        self.x = x
        self.y = y
        self.X = X
        self.Y = Y
        
        self.label = label
        
        self.scf_threshold = scf_threshold #threshold for explained variance until which result maps are plotted
        self.dpi  = 150 #output dpi for plotting
        self.ext = 'pdf' #file extension for plotting
        self.use_basemap = False
        
    def __get_valid_timeseries(self,x):
        '''
        get only points where all
        timeseries are valid
        
        x: masked array [time,position]
        
        '''
        sx = x.sum(axis=0) #calculate sum. If all values are valid, then no a float should be there, else Nan
        mx = ~np.isnan(sx) #masked array
        m1 = mx.data; m1[mx.mask]=False #apply also mask of masked array
        
        return x[:,m1],m1
        
    
    def __detrend_time(self,x):
        '''
        given a variable x[time,position]
        the data is linear detrended
        '''
        
        if x.ndim !=2:
            raise ValueError, 'Invalid shape for detrending'
        
        n,m = x.shape
        
        for i in range(m):
            h = x[:,i].copy()
            h = plt.detrend_linear(h)
            x[:,i] = h.copy()
        
        return x
        
        
    def svd_analysis(self):
        
        #/// perform SVN only for data points which are valid throughout entire time series ///
        x,mskx = self.__get_valid_timeseries(self.x)
        y,mxky = self.__get_valid_timeseries(self.y)
        
        print ' Dimensions for x in SVN: ', x.shape
        print ' Dimensions for y in SVN: ', y.shape
        
        #/// detrend the data for each grid point ///
        x = self.__detrend_time(x)
        y = self.__detrend_time(y)
        
        #calculate covariance matrix
        C = dot(x.T,y) #this covariance matrix does NOT contain the variances of the individual grid points, but only the covariance terms!
        self.C = C
        
        #singular value decomposition
        print '   Doing singular value decomposition ...'
        U, s, V = linalg.svd( C )
        L = linalg.diagsvd(s, len(C), len(V) ) #construct diagonal maxtrix such that U L V.T = C; this is somewhat python specific

        #expansion coefficients (time series)
        A = dot(x,U); B = dot(y,V)
        
        #/// store results
        self.U = U
        self.V = V
        self.L = L
        self.A = A
        self.B = B
        self.scf = (s*s) / sum(s*s) #fractions of variance explained CAUTION: not properly described in manual if squared or not!
        
        self.__get_mode_correlation() #calculate correlation between modes
        
        
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
        
        

        
        
    def plot_var(self,ax=None,filename=None,maxvar=1.):
        '''
        plot explained variance
        
        ax (optional) specifies axis to plot to
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
        
        
        #~ ax.yaxis.label.set_color('blue')
        #~ ax1.yaxis.label.set_color('red')
        #~ ax2.yaxis.label.set_color('green')
        
        
        ax .tick_params(axis='y', colors='blue')
        ax1.tick_params(axis='y', colors='red')
        ax2.tick_params(axis='y', colors='green')
        
        if filename != None:
            oname = filename + '_mode_var.' + self.ext
            ax.figure.savefig(oname,dpi=self.dpi)
            
            

    def plot_correlation_map(self,mode,ctype='hetero',ax1in=None,ax2in=None,pthres=1.01,plot_var=False,filename=None,region1=None,region2=None):
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
        
        

        
        types = ['homo','hetero']
        if ctype in types:
            pass
        else:
            raise ValueError, 'Invalid ctype'
            
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
        
        def plot_cmap(R,ax,title,vmin=-1.,vmax=1.,plot_var=False,use_basemap=False,region=None,cmap='RdBu_r'):
            '''
            R data object
            '''
            
            if plot_var:
                O = R.copy()
                O.data = O.data*O.data
                O.label = 'variance'
            else:
                O = R.copy()
                O.label = 'correlation'
            
            map_plot(O,use_basemap=use_basemap,ax=ax,region=region,cmap_data=cmap,shift=True,vmin=vmin,vmax=vmax)
            
        
        
        
        
        print 'Mode list for correlatio maps: ', mode_list
        #/// calculate correlations and do plotting
        for i in mode_list:
            
            if ax1in == None:
                fig=plt.figure()
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                ax4 = fig.add_subplot(212)
            else:
                ax1 = ax1in
                ax2 = ax2in
                ax4 = ax1in.figure.add_subplot(313)
            
            
            if ctype == 'homo':
                Rout1,Sout1,Iout1,Pout1 = self.X.corr_single(self.A[:,i],pthres=pthres)
                Rout2,Sout2,Iout2,Pout2 = self.Y.corr_single(self.B[:,i],pthres=pthres)
            
            if ctype == 'hetero':
                Rout1,Sout1,Iout1,Pout1 = self.X.corr_single(self.B[:,i],pthres=pthres)
                Rout2,Sout2,Iout2,Pout2 = self.Y.corr_single(self.A[:,i],pthres=pthres)
            

            
            plot_cmap(Rout1,ax1,ctype,plot_var=False,use_basemap=self.use_basemap,region=region1,vmin=-0.8,vmax=0.8,cmap='RdBu_r') #correlation field 1
            plot_cmap(Rout2,ax2,ctype,plot_var=False,use_basemap=self.use_basemap,region=region2,vmin=-0.8,vmax=0.8,cmap='RdBu_r') #correlation field 2
            plot_cmap(Rout2,ax3,ctype,plot_var=True,use_basemap=self.use_basemap,region=region2,vmin=0.,vmax=1.,cmap='YlOrRd')  #explained variance field 2
            
            ax1.figure.suptitle(self.label + ': Mode: #' + str(i) + ' (scf: ' + str(round(self.scf[i],2)) + ')' +   ' - ' + ctype)
            
            self.plot_expansion_correlation(i,ax=ax4)
            
            #--- save figure
            if filename != None:
                oname = filename + '_' + ctype + '_mode_' + str(i) + '.' + self.ext
                ax1.figure.savefig(oname,dpi=self.dpi)
        
    def _get_variance_field(self,X,E,mode,pthres=1.01):
        '''
        calculate variance field for a particular mode
        
        X: Data object of field X
        E: expansion cofficient to be used for correlation calculation
        '''
    
        Rout,Sout,Iout,Pout = X.corr_single(E[:,mode],pthres=pthres)
        
        Rout.data = Rout.data*Rout.data
        
        
        return Rout #return squared correlation to get variance
        
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
        
        
    def print_mode_statistic(self,filename=None):
        '''
        print statistic of modes
        '''
        
        sep = ' & '
        rnd = 2
        
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
                    
    
    def plot_expansion_correlation(self,mode,ax=None):
        '''
        plot correlation and time series of expansion coeffcients
        '''
        if ax == None:
            fig=plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = ax
        
        ax.plot(self.A[:,mode]/np.std(self.A[:,mode]),label='A')
        ax.plot(self.B[:,mode]/np.std(self.B[:,mode]),label='B')
        c = np.corrcoef(self.A[:,mode],self.B[:,mode])[0][1]
        plt.legend()
        ax.set_title('expansion coefficient #' + str(mode) + ' (r=' + str(round(c,2)) + ')')
        ax.set_xlabel('time')
        ax.set_ylabel('normalized expansion coefficient')
        
        

        

            
        
        
        








class Diagnostic():
    def __init__(self,x,y=None):
        '''
        diagnostic for one or multiple data sets
        
        x,y of type Data
        '''
        
        self.x = x
        if y != None:
            self.y = y
    
    def get_n(self):
        '''
        return the number of valid samples
        '''
        
        xm = ~self.xvec.mask #vector with valid sample
        ym = ~self.yvec.mask
        m = xm & ym
        return sum(m)
        
    
    
    def get_rmse_value(self):
        '''
        calculate root-mean-squared error
        
        @param self: Diagnostic object
        @type self : Diagnostic object        
        '''
        return np.sqrt(np.mean((self.xvec - self.yvec)**2))
        
    def lagged_correlation_vec(self,lags,pthres=1.01):
        '''
        lagged correlation for two vectors
        
        x,y Data objects, where data needs to have been pre-processed
        to be a single vector
        
        lags: list of lags
        '''
        
        if self.x.data.shape != self.y.data.shape:
            raise ValueError, 'Invalid geometries!'
        
        if any(lags < 0.):
            raise ValueError, 'Negative lags currently not supported yet!'
        
        CO = []
        for lag in lags:
            hlpx = self.x.copy(); hlpy = self.y.copy()
            
            if lag > 0:
                #temporal subset data
                hlpx.data = hlpx.data[lag:]; hlpy.data = hlpy.data[:-lag]
            slope, intercept, r_value, p_value, std_err = sci.stats.linregress(hlpx.data,hlpy.data)
            
            if p_value < pthres:
                CO.append(r_value)
            else:
                CO.append(np.nan)

        #~ #-- find best lag
        #~ best_lag=abs(np.asarray(CO)).argmax(axis=0).astype('float')
        #~ print 'BEST LAG: ', best_lag
        
        CO = np.asarray(CO)
        
        return CO
        
        
        
        
        
        
        
    def get_correlation_value(self):
        c = np.ma.corrcoef(self.xvec,self.yvec)[0][1]
        return c #correlation coefficient 

    def _mat2vec(self,mask = None):
        '''
        concatenate all information into
        a vector and apply a given
        mask if desired
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
            yvec = Y.data.copy()
            yvec.shape = (-1)
            
        self.xvec = xvec
        self.yvec = yvec
        
        
        
        


        
        
            
            
    
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
        x = self.x.data.copy()
        y = self.y.data.copy()
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
            std_x.shape = (n,-1)
            weights.shape = (-1)

            if np.shape(x[0,:]) != np.shape(weights):
                raise ValueError, 'Invalid shape for weights!'
        
            e2 = []
            for i in range(n):
                d = weights * ( (x[i,:]-y[i,:])**2)   / std_x[i,:]
                #~ print std_x[i,:]
                #~ print '*** ', i
                #~ print min(x[i,:]),max(x[i,:])
                #~ print min(y[i,:]),max(y[i,:])
                #~ print ''
                e2.append(np.sum(d)) #sum at end to avoid nan's   #it is important to use np.sum() !!

            e2 = np.asarray(e2)

        print 'CALCULATED REICHLER INDEX: ', e2
        if np.any(np.isnan(e2)):
            raise ValueError, 'Reichler: e2 contains NAN, this happens most likely if STDV == 0'

        return e2
    
    
    def get_correlation(self,lag=0,nlags=None):
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
        
    
    
    
    
    def slice_corr(self,timmean=True):
        '''
        perform correlation analysis for
        different starting times and length
        of the correlation period
        
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
        x.shape = (n,-1) #size [time,ngridcells]
        y.shape = (n,-1)

        R=np.ones((n,n))*np.nan
        P=np.ones((n,n))*np.nan
        L=np.ones((n,n))*np.nan
        S=np.ones((n,n))*np.nan
        
        #--- perform correlation analysis
        print '   Doing slice correlation analysis ...'
        i1 = 0
        while i1 < n-1: #loop over starting year
            i2 = i1 + 2
            #print i1
            #- loop over different lengths
            while i2 < len(x)-1:
                length = i2-i1
                #print length
                #print 'timmean',timmean
                #~ print i1,i2
                #~ print x[i1:i2,:]
                #~ print y[i1:i2,:]
                
                #- reshape data to valid vector
                
                #todo
                #~ it makes a considerable differencec if
                #~ timmean is applied or not !!!
                
                
                if timmean:
                    ''' temporal mean -> all grid cells only (temporal mean) '''
                    #print 'drin'
                    xdata = x[i1:i2,:].mean(axis=0)
                    ydata = y[i1:i2,:].mean(axis=0)
                    
                    xmsk  = xdata.mask; ymsk = ydata.mask
                    msk = xmsk | ymsk
                    
                    #~ print 'shapex: ', np.shape(xdata)
                    #~ print 'shapey: ', np.shape(ydata)
                    
                    #print 'MASK: ', sum(msk), sum(xmsk),sum(ymsk)
                else:
                    ''' all grid cells at all times '''
                    #~ print i1,i2
                    #~ print x.data.shape, y.data.shape, type(x.data)
                    xdata = x.data[i1:i2,:]; ydata = y.data[i1:i2,:]
                    xmsk  = x.mask[i1:i2,:]
                    ymsk  = y.mask[i1:i2,:]
                    msk   = xmsk | ymsk
                    
                    #~ print 'xdata: ', xdata
                    #~ print x.data
                    #~ print i1,i2
                    #~ print x.data[i1:i2,:]
                    #~ stop
                
                
                #~ print x
                #~ print x.data
                
                #~ print 'x1', xdata, i1, i2, x.data
                #~ print x.data.shape
                #~ print 'y1', ydata
                #~ print xdata.flatten()
                #~ stop
            
                xdata = xdata[~msk].flatten()
                ydata = ydata[~msk].flatten()
                
                #print 'Final data',len(xdata),len(ydata)
                
                slope, intercept, r, p, stderr = sci.stats.linregress(xdata,ydata)
                R[length,i1] = r
                P[length,i1] = p   
                L[length,i1] = length
                S[length,i1] = slope
                
                #~ print 'r: ', r
                #~ print 'x', xdata
                #~ print 'y', ydata
               
                i2 += 1
            i1 += 1
    
        self.slice_r = R
        self.slice_p = P
        self.slice_length = L
        self.slice_slope = S
    


    def _set_year_ticks(self,years,ax,axis='x'):
        '''
        set ticks of timeline with
        yearly ticks
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
            ax.set_xticklabels(oticks)
        elif axis == 'y':
            ax.set_yticklabels(oticks)
        else:
            raise ValueError, 'Invalid axis (set_year_ticks)'


    def plot_slice_correlation(self,pthres = 1.01):
        '''
        plot slice correlation results
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
        ax1=fig.add_subplot(221)
        ax2=fig.add_subplot(222)
        ax3=fig.add_subplot(223)
        ax4=fig.add_subplot(224)
        #ax1 = fig.add_axes([0.07, 0.07, 0.4, 0.8])
        #ax2 = fig.add_axes([0.55, 0.07, 0.4, 0.8])
        
        r_data      = self.slice_r.copy()
        p_data      = self.slice_p.copy()
        length_data = self.slice_length.copy()
        slope_data  = self.slice_slope.copy()
        
        msk = p_data > pthres
        r_data[msk]      = np.nan
        p_data[msk]      = np.nan
        length_data[msk] = np.nan
        slope_data[msk]  = np.nan
        
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
        
        #pl.figure()
        #pl.imshow(P,interpolation='nearest')
        #~ CP = ax.contour(P,[0.01,0.05,0.1],linewidths=2) #,colors='black')
        #~ ax.clabel(CP, inline=1, fontsize=10) #,colors='black')
#~ 
        #~ 
        #~ ax2.plot(x,label='x')
        #~ ax2.plot(y,label='y')
        #~ ax2.set_title('Region: ' + reg.label)
        #~ ax.set_xlabel('years')
        #~ ax2.grid()
        #~ ax2.legend()

    
    

