# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 16:45:44 2015

@author: Felix Darvas

script to demonstrate the application of my implementation of the Gregory-Laredo algorithm 
"""

from simulate_arrival_times import simulate_arrival_times
from GL_algorithm import compute_GL,compute_bin
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    psx=False # parallel execution
    
    w0=0.0004  # simulate slow frequency signal
    phase0=1.0 # with phase offset 1.0
     
    # simulate constant rate process as test of Null - hypothesis
    T1=simulate_arrival_times() 
    
    # simulate periodic rate process as positive test
    T2=simulate_arrival_times(w=w0,phase=phase0,noise_rate=0.0000)
    
    #make bin histograms
    n1=compute_bin(T1,m=7,w=w0,p=phase0)
    n2=compute_bin(T2,m=7,w=w0,p=phase0)
    
    # plot the bin histograms 
    ind=np.arange(7) # x-axis for the bin histogram
    width = 0.35 # the width of the bars
    
    # constant rate process
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, n1, width, color='b')
    ax.set_ylabel('bin count')
    ax.set_title('phase bin histogramm for constant rate process')
    
    # periodic process
    fig, ax = plt.subplots()
    rects2 = ax.bar(ind, n2, width, color='b')
    ax.set_ylabel('bin count')
    ax.set_title('phase bin histogramm for periodic rate process')
    
    # test GL algorithm for constant rate process
    O_period1,p_period1,m_opt1,S1,w1,w_peak1,w_mean1,w_conf1=compute_GL(T1,parallel=psx)
    
    fig, ax = plt.subplots()
    ax.plot(w1,S1)
    ax.set_ylabel('probability')
    ax.set_xlabel('f [rad/s]')
    ax.set_title('spectrum for constant rate process')
    
    print ('Likelihood of periodic process =%3.2f %% most likely frequency %e mean frequency %e 95 %% confidence interval = [%e  %e]\n') % (p_period1*100,w_peak1,w_mean1,w_conf1[0],w_conf1[1])
    
    
    
    O_period2,p_period2,m_opt2,S2,w2,w_peak2,w_mean2,w_conf2=compute_GL(T2,parallel=psx)
    
    fig, ax = plt.subplots()
    ax.plot(w2,S2)
    ax.set_ylabel('probability')
    ax.set_xlabel('f [rad/s]')
    ax.set_title('spectrum for periodic rate process')
    print ('Likelihood of periodic process =%3.2f %% most likely frequency %e mean frequency %e 95 %% confidence interval = [%e %e]\n') % (p_period2*100,w_peak2,w_mean2,w_conf2[0],w_conf2[1])
    plt.show()

