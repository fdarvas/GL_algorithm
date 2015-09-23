# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 17:20:28 2015

@author: Felix Darvas

computes phase-bin membership for a list of arrival times
inputs: Tlist   - numpy integer array of arrival times
        m       - number of phase bins
        w       - frequency in rad/s
        p       - phase in rad
"""

import numpy as np

def compute_bin(Tlist,m,w,p):
    n=np.zeros(m,'int')
    j=np.floor(m*np.mod(w*Tlist+p,2*np.pi)/(2*np.pi))
    j.astype(int)
    for u in range(0,m):
        n[u]=np.size(np.extract(j==u,j))
    return n