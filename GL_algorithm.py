# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:13:40 2015

@author: Felix Darvas

compute the Gregory-Laredo algorithm on arival times

This function computes the likelihood of a set of arrival times originating 
from a periodic system rather than constant rate (poisson) background noise

based on

Gregory, P. C. and Thomas. J. Loredo, 1992, 
"A New Method For The Detection Of A Periodic Signal Of Unknown Shape And Period" 
in 
The Astrophysical Journal, Astrophysical J., 398, p.146

inputs:

Tlist    -  list of arrival times, numpy int array
m_max    -  max number of bins typically 12-15, we use 12 as default
w_range  -  frequency range to scan numpy float array of frequency values 
           default is  w_lo=20*pi/T at delta w = pi/T to w_hi=pi*N/T
           where N=#arrival times, T=observation time
ni       - number of integration steps, default ni=10
parallel - use paralel execution - default is off

outut:
O_period - Odds ratio for a periodic process vs. constant rate process
p_period - probability of a periodic process 0<=p_period<=1
m_opt    - optimal bin size 1<= m_opt <=m_max
S        - The probability spectrum 
w        - The frequency range for S
w_peak   - the expected frequency of the process
w_conf   - 95% confidence interval of w_peak

"""

import numpy as np 
import multiprocessing as mp 
import functools

def compute_bin(Tlist,m,w,p):
    n=np.zeros(m,'int')
    j=np.floor(m*np.mod(w*Tlist+p,2*np.pi)/(2*np.pi))
    j.astype(int)
    for u in range(0,m):
        n[u]=np.size(np.extract(j==u,j))
    return n

def compute_factor(N,m,v): 
    # compute m^N *(N+m-1) over N /(2pi v)
    # whic is used to scale the multiplicity function (see eqn.5.25 of the paper )
    # we return the log of the integral scale here to avoid numerical overflow
    f1=N*np.log(m) # log of m^N
    f2=np.sum(np.log(np.arange(1,N+m))) # log of (N+m-1)!
    f3=np.sum(np.log(np.arange(1,m+1))) # log of m!
    f=f1+f3-f2-np.log(2*np.pi*v)
    return f
    
def compute_W_scaled(Tlist,m,w,phase,factor):
    # compute the scaled multiplicity 1/Wm(w,phase) (see eqn. 5.25)
    # note that for large arrival time numbers the multiplicity becomes
    # excessively large. Since W_m(phase,w) is never needed explicitly,
    # we use the scaled version. 
    # input factor is the log of the actual factor
    n=compute_bin(Tlist,m,w,phase) # find bin histogram for a given bin number m, frequency w and phase 
    f=0    
    for i in range(0,m):
        f=f+np.sum(np.log(np.arange(2,n[i]+1)))
    y=np.exp(f+factor)
    return y

def compute_Om1(w,Tlist,m,factor,ni):
    #compute  specific odds-ratios (eqn. 5.25)
    p=np.arange(0,ni)/float(ni)*2*np.pi/m # intgration range, only integrate over a single bin, as values of the integral repeat over bins       
    y=np.zeros(np.size(p),'float') # array to hold values of W_scaled over the integration range
    for i in range(0,np.size(y)):
        y[i]=compute_W_scaled(Tlist,m,w,p[i],factor)
    return np.trapz(y,p)*m # return intregrated W_Scaled


def compute_Om1wPar(Tlist,m_max,w,fa,ni): # compute odds-ratios for bins and frequencies
	# parallel version 
    Om1w=np.zeros((m_max,np.size(w)),'float') # odds ratio matrix
    pool=mp.Pool() # use all workers
    
    for m in range(0,m_max):
		Om1w[m,:]=pool.map(functools.partial(compute_Om1,Tlist=Tlist,m=(m+1),factor=fa[m],ni=ni),w)
    return Om1w
    
def compute_Om1w(Tlist,m_max,w,fa,ni): # compute odds-ratios for bins and frequencies
    Om1w=np.zeros((m_max,np.size(w)),'float') # odds ratio matrix
    for m in range(0,m_max):
        for wi in range(0,np.size(w)):
            Om1w[m,wi]=compute_Om1(w[wi],Tlist,m+1,fa[m],ni)
    return Om1w
    
                    

def compute_GL(Tlist,m_max=12,w_range=None,ni=10,parallel=False):
    # initialize output values
    O_period=None
    p_period=None
    m_opt=None
    S=None
    w=None
    w_peak=None
    w_mean=None
    w_conf=None
    N=float(np.size(Tlist)) # need float value to avoid int/int
    if N>0:
        # compute GL algorithm
        v=m_max-1
        T=float(np.max(Tlist)) # duration of the observation
        if w_range is None: # use default frequencies
            w_hi=np.pi*N/T # max default frequency
            
            w_lo=np.minimum(20,N/10)*np.pi/T # min default frequency
            dw=np.pi/T/10 # step size
            w=np.arange(w_lo,w_hi,dw)
            if np.size(w)<2:
                print "error "
                raise ValueError('bad arrival time list')
        else: # use user supplied frequency vector
            w=w_range 
            w_hi=np.max(w_range)
            w_lo=np.min(w_range)
        if w_lo==0:
            # cannot have w_lo =0
            print ("minimum frequency cannot be 0!\n")
            return
       
        fa=np.zeros(m_max)
        for m in range(0,m_max): # precompute factors for each m
            fa[m]=compute_factor(N,m+1,v)    
        if parallel:    
			Om1w=compute_Om1wPar(Tlist,m_max,w,fa,ni)
        else:
            Om1w=compute_Om1w(Tlist,m_max,w,fa,ni)
            
        pw=1/w/np.log(w_hi/w_lo)
        O1m=np.zeros(m_max)
        for i in range(0,m_max): # intgreate odd ratios across frequencies
            O1m[i]=np.trapz(pw*Om1w[i],w)
            
        m_opt=np.argmax(O1m) # find optimum bin number, i.e largest odds-ratio    
        S=Om1w[m_opt]/w # compute Spectral probability
        m_opt=m_opt+1 # start bin index with 1
        C=np.trapz(S,w) # compute normalization
        S=S/C # normalized probability
        O_period=np.sum(O1m) # integrated odds ratio
        p_period=O_period/(1+O_period) # likelihood of periodic event
        cdf=np.array(S)
        for i in range(0,np.size(S)):
            cdf[i]=np.trapz(S[0:i],w[0:i])
        wr=np.extract(np.logical_and(cdf>.025, cdf<.975),w)
        w_peak=w[np.argmax(S)]
        w_mean=np.trapz(S*w,w)
        if np.size(wr)>0:
            w_conf=[np.min(wr),np.max(wr)]
        else:
            w_conf=[w_peak,w_peak]
        return O_period,p_period,m_opt,S,w,w_peak,w_mean,w_conf
    else:
        # throw an error
        print ("No valid arrival time array provided!\n")
        return O_period,p_period,m_opt,S,w,w_peak,w_mean,w_conf
        
