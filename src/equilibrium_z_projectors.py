"""
Calculates correlation functions for sigma z operators in the TFIM at criticalicity and equilibrium (Eg. when g=1)

This file modifies sigma_general, P_n, and P_n_correlations from equilibrium_critical_correlations.py
"""

import kblock_ising_model as kb #Useful generic functions 
from collections import Counter
import numpy as np
from numpy import cos, sin, exp, pi, sqrt, abs, pi

#Correlation Function
#Definition of this comes from sachdev
def Dn(n):
    return 2/pi/(1-2*n)
###Typical Sigma General Function
def BA(j):
    return Dn(-j+1)
###Sigma General
def sigma_general(indices):
    def remove_duplicates_in_pairs(vec):
      unique_vals, counts = np.unique(vec, return_counts=True)
      filtered_vals = unique_vals[counts % 2 != 0]
      return filtered_vals.tolist()
    #Sigma matrices on different sites commute
    indices = np.sort(indices)
    #Remove any duplicates as sigma_x^2 = 1
    indices = remove_duplicates_in_pairs(indices)

    #Get string lengths
    A_coords = np.array(indices)
    B_coords = np.array(indices)
    N = len(indices)
    C = np.zeros((N,N))
    for nx in range(N):
        for ny in range(N):
            Bx = B_coords[nx]
            Ay = A_coords[ny]
            Nd = Bx-Ay+1
            C[nx,ny] = -Dn(Nd)
    return np.linalg.det(C)

   
#Typical Pn function
def P_n(n):
    """
    P_n scales directly with the number of terms since even small odd sigma correlations need larger support to calcualte.
    Taking out odd terms works to make easier, but scaling is still 2^n.
    Minor improvements can still be had though.
    """
    indices = [i for i in range(0,n)]
    terms = kb.all_combinations(indices)
    ### Create Terms and degeneracy
    x_c = []
    for a in terms[1:]:
        a = np.array(a)
        x_c.append(tuple(a - a[0]))  # convert to tuple for hashing
    counter = Counter(x_c)
    # Extract vecs and counts as separate arrays/lists
    vecs = list(counter.keys())
    counts = np.array(list(counter.values()))
    ##
    dat = []
    for term in range(len(vecs)):
        indices = vecs[term]
        degen = counts[term]
        dat.append(sigma_general(indices)*degen)
    dat.append(1)
    #All terms have equal weight. 
    return np.sum(dat)/2**n

#P_n correlations do not calculate connected correlator at all
def P_n_correlations(n,l):
    #Term 1 
    indices = [i for i in range(0,n)]
    #Term 2 
    indices +=[i for i in range(n+l-1,2*n+l-1)]
    terms = kb.all_combinations(indices)
   # print(terms)
     ###
    x_c = []
    for a in terms[1:]:
      
        a = np.array(a)
        x_c.append(tuple(a - a[0]))  # convert to tuple for hashing
    counter = Counter(x_c)

    # Extract vecs and counts as separate arrays/lists
    vecs = list(counter.keys())
    counts = np.array(list(counter.values()))
    #
    dat = []
    for term in range(len(vecs)):
        indices = vecs[term]
        degen = counts[term]

        dat.append(sigma_general(indices)*degen)
    dat.append(1)
    #All terms have equal weight. 

    return np.sum(dat)/2**(2*n)#,terms
