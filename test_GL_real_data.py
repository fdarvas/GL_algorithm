# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 12:05:35 2015

@author: Felix Darvas

some tests with real astronomical data
data source:
http://astrostatistics.psu.edu/datasets/Chandra_flares.html


"""

from GL_algorithm import compute_GL,compute_bin
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import inspect

def import_data(fname):
    T=np.loadtxt(fname,skiprows=1)# ignore first line, which is just column headers
    T=T[:,0] # only need first col, which contains arrival times
    return T



file_list=['COUP263.dat','COUP551.dat'] 


if __name__ == '__main__':
    # parallel execution
    px=True
    print "parallel execution\n"
    
    # find the data directory    
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    current_dir = os.path.dirname(os.path.abspath(filename))
    
    local_path="data" # this folder should containt the files in file_list
    
    
    for fname in file_list:
        fullfile=os.path.join(current_dir,local_path,fname)
        
        T=import_data(fullfile) # get astro-data from file

        t0=time.time() # benchmark performance
        O_period,p_period,m_opt,S,w,w_peak,w_mean,w_conf=compute_GL(T,parallel=px) # run GL,  parallel execution
        t1=time.time()
        print "total time used = %f with parallel execution\n" % (t1-t0)
    
        n=compute_bin(T,m=m_opt,w=w_peak,p=0) # compute resulting bin histogram
        print ('File:%s - Likelihood of periodic process =%3.2f %% most likely frequency %e mean frequency %e 95 %% confidence interval = [%e %e]\n') % (fname,p_period*100,w_peak,w_mean,w_conf[0],w_conf[1])
            
    
    
    
    # serial execution
    print "serial execution\n"
    px=False
    
    
    for fname in file_list:
        fullfile=os.path.join(current_dir,local_path,fname)
        
        T=import_data(fullfile) # get astro-data from file
        t0=time.time() # benchmark performance
        O_period,p_period,m_opt,S,w,w_peak,w_mean,w_conf=compute_GL(T,parallel=px) 
        t1=time.time()
        print "total time used = %f with serial execution\n" % (t1-t0)
        n=compute_bin(T,m=m_opt,w=w_peak,p=0) # compute resulting bin histogram
        
        #make bar graph for bin histogram
        ind=np.arange(np.size(n)) # x-axis for the bin histogram
        width = 0.35 # the width of the bars
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, n, width, color='b')
        ax.set_ylabel('bin count')
        ax.set_title('phase bin histogramm for '+fname)
        
        #plot the probability spectrum
        fig, ax = plt.subplots()
        ax.plot(w,np.log(S))
        ax.set_ylabel('log probability')
        ax.set_xlabel('f [rad/s]')
        ax.set_title('phase spectrum for periodic rate process - '+fname)
        print ('File:%s - Likelihood of periodic process =%3.2f %% most likely frequency %e mean frequency %e 95 %% confidence interval = [%e %e]\n') % (fname,p_period*100,w_peak,w_mean,w_conf[0],w_conf[1])
            
    
    plt.show()
