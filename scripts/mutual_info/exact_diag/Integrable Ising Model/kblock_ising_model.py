#Imports
import numpy as np
import numpy.linalg as la
from itertools import combinations
from scipy.linalg import toeplitz
from collections import Counter

#TFIM Functions
def k_vals(L):
    return (2 * np.arange(-L//2, L//2) + 1) * np.pi / (L)
def epsilon(k,g,J=1):
    """Calculates disperson from TFIM
    inputs: k (np.array), g (float), J (float)

    returns: epsilon (np.array)
    """
    return 2*J*np.sqrt( (g-np.cos(k))**2+np.sin(k)**2)

def U(k, g):
    """Returns normalized eigenvectors in k space for TFIM
    inputs: k (np.array), g (float)

    returns: [u,v] where u={uk} (np.array) and v= {vk} (np.array) are the normalized eigenvectors


    Example

    U(k,inf)[0] = u = [1 , 1 , 1....]
    """
    u = g - np.cos(k) + np.sqrt(g**2 - 2 * g * np.cos(k) + 1)
    v = np.sin(k)
    norm = np.sqrt(u**2 + v**2)
    return u / norm, v / norm

#Correlation Functions


def ca(U, k, l):
    """inputs: U = [u,v], k (np.array), l (integer)
    returns: np.array
    """
    #Unpack Variables
    u, v = U
    # surviving term against bogliouv vaccuum 
    amp = np.abs(v)**2
    #Fourier transform phase
    phase = np.exp(1j * k * l)
    return np.mean(amp * phase)
def ac(U, k, l):  
    """inputs: U = [u,v], k (np.array), l (integer)
    returns: np.array
    """
    #Unpack Variables
    u, v = U
    amp = np.abs(u)**2
    phase = np.exp(-1j * k * l)
    return np.mean(amp * phase)
def aa(U, k, l):
    """inputs: U = [u,v], k (np.array), l (integer)
    returns: np.array
    """
    #Unpack Variables
    u, v = U
    amp = np.conj(v) * u
    phase = np.exp(-1j * k * l)
    return -1j * np.mean(amp * phase)
def cc(U, k, l):
    """inputs: U = [u,v], k (np.array), l (integer)
    returns: np.array
    """
    #Unpack Variables
    u, v = U
    amp = np.conj(u) * v
    phase = np.exp(-1j * k * l)
    return 1j * np.mean(amp * phase)
def AA(args):
    """input: args = [U,k,l]"""
    return cc(*args)+ca(*args)+ac(*args)+aa(*args)
def BB(args):
    """input: args = [U,k,l]"""
    return aa(*args)-ac(*args)-ca(*args)+cc(*args)
def AB(args):
    """input: args = [U,k,l]"""
    return (-aa(*args)+ac(*args)-ca(*args)+cc(*args))
def BA(args):
    """input: args = [U,k,l]"""
    return (-aa(*args)-ac(*args)+ca(*args)+cc(*args))
def D(N,U,k):
    """Using Definition in 4.72 of Sachdev.
    
    N = i-j+1 but we will also set i = 0 and j=n so 
    N = 1-n"""
    n = 1-N
    args = [U,k,n]
    return BA(args)


def sigma_general(indices,T,k):
    def remove_duplicates_in_pairs(vec):
      unique_vals, counts = np.unique(vec, return_counts=True)
      filtered_vals = unique_vals[counts % 2 != 0]
      return filtered_vals.tolist()
    #Sigma matrices on different sites commute
    indices = np.sort(indices)
    #Remove any duplicates as sigma_x^2 = 1
    indices = remove_duplicates_in_pairs(indices)
    #Bs sit on odd sites
    odd_sites = np.array(indices[::2])
    #As site on even sites
    even_sites= np.array(indices[1::2])
    #Get string lengths
    JW_string_lengths = even_sites-odd_sites
    #Sum of string lengths is size of matrix needed
    N = sum(JW_string_lengths)
    #Fill in indices for strings
    R = []
    for i in range(0, len(indices), 2):
        start = indices[i]
        end = indices[i+1]
        R.extend(range(start, end+1))

    A_coords = [x for x in R if x not in odd_sites]
    B_coords = [x for x in R if x not in even_sites]
    N = len(A_coords)

    # Compute first column: D(B_coords[0] - A_coords + 1)
    first_col = np.array([D(B_coords[0] - A_coords[j] + 1, T, k) for j in range(N)])

    # Compute first row: D(B_coords[i] - A_coords[0] + 1)
    first_row = np.array([D(B_coords[i] - A_coords[0] + 1, T, k) for i in range(N)])

    # Construct Toeplitz matrix
    C = toeplitz(first_col, first_row)

    return la.det(C)
####PROJECTORS 


def binomial_expansion(indices):
    """
    Computes the pairs of combinations of indices.
    """
    pairs = list(combinations(indices, 2))
    return pairs

def all_combinations(indices):
    """
    Computes all combinations of all lengths of indices.
    """
    x = []
    for r in range(len(indices) + 1):
        x.extend(combinations(indices, r))
    return list(x)

def unique_elements_and_frequencies(vec):
    unique_vals, freqs = np.unique(vec, return_counts=True)
    return unique_vals, freqs


def P_n(n,U,k, even = True):
    """
    P_n scales directly with the number of terms since even small odd sigma correlations need larger support to calcualte.
    Taking out odd terms works to make easier, but scaling is still 2^n.
    Minor improvements can still be had though.
    """
    #For most cases we use periodic boundary conditions
    #That said, watch out for this definition of indices
    indices = [i for i in range(0,n)]
    terms = all_combinations(indices)

    ###
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

        if len(indices)%2 == 1:
            if even ==True:
                dat.append(0)
            else:
                constant = len(k)//2
                indices = list(indices) + [x + constant for x in indices]
                dat.append(np.sqrt(sigma_general(indices,U,k))*degen)
        else:
            dat.append(sigma_general(indices,U,k)*degen)
    dat.append(1)
    #All terms have equal weight. 
    return np.sum(dat)/2**n



def P_n_correlations(n,l,U,k, even = True):
    #Term 1 
    indices = [i for i in range(0,n)]
    #Term 2 
    indices +=[i for i in range(n+l,2*n+l)]
    terms = all_combinations(indices)
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

        if len(indices)%2 == 1:
            if even ==True:
                dat.append(0)
            else:
                constant = len(k)//2
                indices = list(indices) + [x + constant for x in indices]
                dat.append(np.sqrt(sigma_general(indices,U,k))*degen)
        else:
            dat.append(sigma_general(indices,U,k)*degen)
    dat.append(1)
    #All terms have equal weight. 

    return np.sum(dat)/2**(2*n)#,terms
