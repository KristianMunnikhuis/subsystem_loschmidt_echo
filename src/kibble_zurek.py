"""
Time-dependent equations for the Bogoliouv fermions starting at t = -inf going to t = t.

At t= 0, it is at the critical point (g=1)

At t = tau, it is in the full ordered x phase (g = 0)

This DOES NOT do free evolution. The analytical equations are derived under some specific initial conditions!

"""

################################################
#                   IMPORTS
################################################
from project_imports import * 
from mpmath import pcfd
import mpmath as mp

mp.mp.dps = 100
mp.mp.maxterms = 100000# number of hypergeometric terms before giving up
################################################
#            TIME DEPENDENT SOLUTIONS
################################################
###These Functions are given by 7.2.87-88 of "Quantum Ising Phases and Transitions"
###They represent the bogoliouv Transformation starting from t0 = - inf to a time tf with a speed related to tau. g goes from infinity to -infinity.



def u_tilde(tau, q, tf):
    """Eq. 7.2.87"""
    a = 1 + mp.cos(q)
    b = mp.sin(q)

    rot = mp.e**(-1j * 3 * mp.pi / 4)
    rot_inv = mp.e**(1j * 3 * mp.pi / 4)
    decay = mp.e**(-mp.pi * b * b * tau / 4)
    prefactor = rot * b * mp.sqrt(tau) * decay
    argument = rot_inv * 2 * (tf - a * tau) / mp.sqrt(tau)
    nu = 1j * b * b * tau - 1

    with mp.workdps(mp.mp.dps + 20):
        return prefactor * mp.pcfd(nu, argument)

def v_tilde(tau, q, t):
    """Eq. 7.2.88"""
    a = 1 + mp.cos(q)
    b = mp.sin(q)

    rot = mp.e**(-1j * 3 * mp.pi / 4)
    rot_inv = mp.e**(1j * 3 * mp.pi / 4)
    decay = mp.e**(-mp.pi * b * b * tau / 4)
    prefactor = -decay
    argument = rot_inv * 2 * (t - a * tau) / mp.sqrt(tau)
    nu = 1j * b * b * tau

    with mp.workdps(mp.mp.dps + 20):
        return prefactor * mp.pcfd(nu, argument)


def g(t,tau):
    return -t/tau + 1



# ################################################
# #            GROUND STATE SOLUTIONS (INCORRECT)
# ################################################

# def u_diag(q,g):
#     a = g+cos(q)
#     omega= sqrt(g*g+2*g*cos(q)+1)

#     denom = sqrt(2*omega*(omega-a))

#     return (a - omega)/denom

# def v_diag(q,g):
#     a = g+cos(q)
#     b = sin(q)
#     omega= sqrt(g*g+2*g*cos(q)+1)

#     denom = sqrt(2*omega*(omega-a))

#     return b/denom


################################################
#           FULL SYSTEM SOLUTIONS
################################################
def k_vals(L):
    """Generates k values >0"""
    x=(2 * np.arange(-L//2, L//2) + 1) * np.pi / (L)
    return x[x>0]


# def state(L,tau,t):
#     """Generates state as [us, vs, k]"""
#     #Generate momentum values
#     k = k_vals(L)
#     us = np.array([u_tilde(tau,q,t) for q in k])
#     vs = np.array([v_tilde(tau,q,t) for q in k])
#     return [us,vs,k]


def state(L, tau, t):
    """Generates state as [us, vs, k]"""
    k = k_vals(L)
    us = np.array([complex(u_tilde(tau, q, t)) for q in k], dtype=np.complex128)
    vs = np.array([complex(v_tilde(tau, q, t)) for q in k], dtype=np.complex128)
    k  = np.array(k, dtype=np.float64)
    return [us, vs, k]

def ground_state(L,tau,t):
    #Is supposed to generate ground state
    k = k_vals(L)
    us = np.array([v_diag(q,g(t,tau)) for q in k])
    vs = np.array([u_diag(q,g(t,tau)) for q in k])
    norm = sqrt(u**2+v**2)
    us /= norm
    vs /= norm
    return [us,vs,k]

################################################
#           CORRELATION FUNCTIONS
################################################

# note: I had an EXTREMELY tough time on this one. Pay attention to the signs and phase of the terms! G must be a real number! 
def G(r,state):
    """r = integer
    state = [us, vs, k]
    calculates {bi,aj} correlation of r=j-i
    """
    u = state[0]
    v = state[1]
    k = state[2]
    #Equation 7.2.26a
    Gr = 2*cos(k*r)*(np.abs(u)**2-np.abs(v)**2)-2j*sin(k*r)*(np.conj(u)*v-np.conj(v)*u)
    return -np.mean(Gr)/2

def D(n,state):
    """Strange Convention, n = r-1"""
    r = n+1
    return G(r,state)

def sigma_general(indices,state):
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
    C = np.zeros((N,N),dtype=complex)
    for nx in range(N):
        for ny in range(N):
            Bx = B_coords[nx]
            Ay = A_coords[ny]
            Nd = Bx-Ay+1
            C[nx,ny] = D(Nd,state)
    return np.linalg.det(C)

################################################
#           PPROJECTORS
################################################
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


def P_n(n,state, even = True):
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
                dat.append(sqrt(np.abs(sigma_general(indices,U,k)))*degen)
        else:
            dat.append(sigma_general(indices,state)*degen)
    #Manually add the 1 
    dat.append(1)
    #All terms have equal weight. 
    return np.sum(dat)/2**n

def delta(state,args,n):
    """Args = [U,k]"""
    #Caluclate the difference between KZ state and ground state given by KB
    KZ = P_n(n,state)
    GS = kb.P_n(n,*args)
    return np.abs(KZ-GS)



#P_n correlations do not calculate connected correlator at all
def P_n_correlations(n,l,state, even = True):
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

        if len(indices)%2 == 1:
            if even == True:
                dat.append(0)
            else:
                constant = len(k)//2
                indices = list(indices) + [x + constant for x in indices]
                dat.append(np.sqrt(sigma_general(indices,state))*degen)
        else:
            dat.append(sigma_general(indices,state)*degen)
    dat.append(1)
    #All terms have equal weight. 

    return np.sum(dat)/2**(2*n)#,terms
