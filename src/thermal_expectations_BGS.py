import sys
import os
sys.path.append("../../../../src/")
#Imports
import numpy as np
import scipy as sp
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, scatter
#Custom made imports
import single_particle_sector as sps
from time import time

###Functions
def E_p(p,L,J,h):
    #Currently only for PBC or ABC
    bc = ["ABC","PBC"]
    bc = bc[p]
    
    #Hamiltonian and spectrum
    H = sps.H_bdg(h,L,J,bc)
    E, U = la.eigh(H)
    E = E[L:]
    #Ground state correlation matrix
    G_gs = sps.G_tfim(U[:,:L])
    F = G_gs[:L,L:]
    G = G_gs[:L,:L]
    M = np.eye(L)-2*(G+F)
    #Determine parity
    n = la.det(M)*(-1)**L * (-1)**2

    return  E,U,n

def Z_p(n,E,beta):
    #Definition given in first source
    Z =  (np.prod(1+np.exp(-beta*2*E))
        +n*np.prod(1-np.exp(-2*beta*E)))
    return Z/2

def ca(n,E,beta,mu):
        #Define Positive Energies
        #Prefactor
        pf = np.exp(-2*beta*E[mu])
        #Remove site energy from Partition
        E = np.delete(E,mu)
        #Calculate
        ca_p = (np.prod(1+np.exp(-2*beta*E))+
               n*np.prod(1-np.exp(-2*beta*E)))
        #Return average of two terms
        return pf*ca_p/2
def ac(n,E,beta,mu):
        #Remove site energy from Partition

        E = np.delete(E,mu)

        ac_p = (np.prod(1+np.exp(-2*beta*E))+
               n*np.prod(1-np.exp(-2*beta*E)))
        return ac_p/2

def G_th(parities,beta,args):
    #Unpack Arguments
    L,J,h = args
    Z = 0
    G = 0
    for p in parities:
        #Spectrum
        E,U,n = E_p(p,L,J,h)
        #Partition Function
        Z += Z_p(n,E,beta)
        #Correlation matrix of Phi Phi_Dagger
        N = np.diag([ac(n,E,beta,mu) for mu in range(L)]+[ca(n,E,beta,mu) for mu in range(L)])
        G+= U@N@U.T.conj()
    return G/Z

def M(Gi):
    
    L = np.shape(Gi)[0]//2
    G = Gi[:L,:L]
    F = Gi[:L,L:]
    return np.eye(L)-2*(G+F)
