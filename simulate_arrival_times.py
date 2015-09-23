# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 11:50:43 2015

@author: Felix Darvas

This function simulates arrival times for either a constant rate (Possion-process)
or periodic rate.

call

T=simulate_arrival_times(N=500, rate=0.001, w=None, phase=None)

inputs:
N           - number of arrival times to simulate
rate        - event rate per sample 0<=rate<1
w           - frequency (rad/s) of periodic event.
phase       - phase (rad) of periodic event
noise_rate  - add some noise to the periodic process
call without w,phase to produce a constant rate processs

output:
numpy float array of arrival times 

"""
import numpy as np
from random import random as rand

def simulate_arrival_times(N=500, rate=0.001, w=None, phase=None,noise_rate=0.001):
    T=np.zeros(N,'float')
    k=0
    nn=0
    while nn<N:
        if w is None or phase is None:
            # constant rate process  
            if rand()<=rate:
              T[nn]=k
              nn=nn+1
        else:
          # periodic process
          # modulate rate by sine     
            u=(np.sin(w*k+phase)+1)/2 
            if rand()<=rate*u or rand()<=noise_rate:
                T[nn]=k
                nn=nn+1
        k=k+1
    return T