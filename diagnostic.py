#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Alexander Loew"
__version__ = "0.0"
__date__ = "0000/00/00"

import numpy as np


class Diagnostic():
    def __init__(self,x,y=None):
        '''
        diagnostic for one or multiple data sets
        
        x,y of type Data
        '''
        
        self.x = x
        if y != None:
            self.y = y
    
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
            #~ if x.ndim !=3:
                #~ sys.exit('Reichler diagnostic: not 3D data!')
            n = len(x)
            #~ print 'Timesteps in Reichler: ', n
            x.shape = (n,-1) #[time,index]
            y.shape = (n,-1)
            std_x.shape = (n,-1)
            weights.shape = (-1)

            if np.shape(x[0,:]) != np.shape(weights):
                raise ValueError, 'Invalid shape for weights!'
        
            e2 = []
            for i in range(n):
                d = weights * ( (x[i,:]-y[i,:])**2)   / std_x[i,:]
                e2.append(np.sum(d)) #sum at end to avoid nan's   #it is important to use np.sum() !!

            e2 = np.asarray(e2)

        return e2
    
    
    def get_correlation(self,lag=0,nlags=None):
        '''
        calculate correlation between two data sets
        '''
        
        print 'Calculating correlation with lag=', lag
        
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
        x.shape = (n,-1)
        y.shape = (n,-1)
    

        

        #~ print s1
        print x.shape
        R=np.zeros(s1[1]*s1[2])*np.nan
        for i in range(s1[1]*s1[2]):
            xx=x[:,i]
            yy=y[:,i]
            msk = (~np.isnan(xx) & ~np.isnan(yy))
            if sum(msk) > 3:
                
                if lag == 0:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(xx[msk],yy[msk])
                else:
                    
                    
                    if nlags == None:
                        nlags = len(xx[msk])
                        if mod(nlags,2) == 0:
                            nlags = nlags + 1
                    
                    r1,lags = NormCrossCorrSlow(xx[msk],yy[msk],nlags=nlags)
                    idx = nlags/2+lag
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
        R = np.reshape(R,np.shape(self.x.data[0,:]))
        
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
            raise KeyError, 'No y-value specified. Processing not possible!'
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
                    xdata = x.data[i1:i2,:]; ydata = y.data[i1:i2,:]
                    xmsk  = x.mask[i1:i2,:]; ymsk  = y.mask[i1:i2,:]
                    msk   = xmsk | ymsk
                    
                xdata = xdata[~msk].flatten()
                ydata = ydata[~msk].flatten()
                
                #print 'Final data',len(xdata),len(ydata)
                
                slope, intercept, r, p, stderr = stats.linregress(xdata,ydata)
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
    


    def __set_year_ticks(self,years,ax,axis='x'):
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
        imr=ax1.imshow(r_data,interpolation='nearest',cmap='RdBu_r')
        ax1.set_title('correlation')
        plt.colorbar(imr,ax=ax1,shrink=0.8)
        ax1.set_xlabel('start year')
        ax1.set_ylabel('correlation period [years]')
        
        #- significance
        imp=ax2.imshow(p_data,interpolation='nearest',cmap='RdBu_r')
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
        ims=ax4.imshow(slope_data,interpolation='nearest',cmap='RdBu_r')
        ax4.set_title('slope')
        plt.colorbar(ims,ax=ax4,shrink=0.8)
        ax4.set_xlabel('start year')
        ax4.set_ylabel('correlation period [years]')        
        
        #/// set tick labels ///
        self.__set_year_ticks(years,ax1,axis='x')
        self.__set_year_ticks(years,ax2,axis='x')
        self.__set_year_ticks(years,ax3,axis='x')
        self.__set_year_ticks(years,ax4,axis='x')
        
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

    
    

