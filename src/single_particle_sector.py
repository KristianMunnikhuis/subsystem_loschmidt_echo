import numpy as np
import scipy as sp
import numpy.linalg as la 
from itertools import combinations

#See The related TFIM_ground_state notebook for some more examples

def G_tfim(psi):
    """
    Gives Correlations matrix, in which LxL matrices G and F live.

    Input
    Psi: (L x 2L) array
    Output
    G: 2Lx2L array
    Example:
    E,U = la.eigh(H)
    G = G_tfim(U[:,:L])
    """
    return psi @ psi.T.conj() 

def h_t(t, tau, h0, hf):
    """Linearly Ramped h
    Input
    t: array or scalar
    tau: scalar
    h0: scalar
    hf: scalar
    
    Output:
    h: typeof(t)
    """
    return h0 + (hf - h0) * (t / tau)

####Correlation Functions of Sigma_x :

#Correlation Functions:
#Given by eq(239-242)
def AA(Gi):
    #239
    L = np.shape(Gi)[0]//2
    G = Gi[:L,:L]
    F = Gi[:L,L:]
    return G+(np.eye(L)-G.T)+F+F.T.conj()
def BB(Gi):
    #240
    L = np.shape(Gi)[0]//2
    G = Gi[:L,:L]
    F = Gi[:L,L:]
    return -G-(np.eye(L)-G.T)+F+F.T.conj()
def AB(Gi):
    #241
    L = np.shape(Gi)[0]//2
    G = Gi[:L,:L]
    F = Gi[:L,L:]
    return G-(np.eye(L)-G.T)-F+F.T.conj()
def BA(Gi):
    #242
    L = np.shape(Gi)[0]//2
    G = Gi[:L,:L]
    F = Gi[:L,L:]
    return -G+(np.eye(L)-G.T)-F+F.T.conj()


def construct_M_matrix(fermion_operators, correlation_function,G):
    """
    fermion_operators: list of tuples (site index, operator type)
    t: time
    correlation_function: function taking (op1, op2, t) and returning ⟨f_i f_j⟩
    """
    l = len(fermion_operators)
    M = np.zeros((l, l), dtype=complex)

    for mu in range(l):
        for nu in range(mu + 1, l):
            op1 = fermion_operators[mu]
            op2 = fermion_operators[nu]
            val = correlation_function(op1, op2,G)
            M[mu, nu] = val
            M[nu, mu] = -val

    return M

def correlation_function(op_string_1, op_string_2, G):
    """ Calculates the correlation function between two operator strings
    Input: 
    op_string_1: [scalar, string]
    op_string_2: [scalar, string]
    G: 2Lx2L matrix

    Output: 
    Scalar value of correlation function

    Example:
    op_string_1 = [1,"A"]
    op_string_2 = [5,"B"]
    returns: AB[1,5]
    """
    i, op1 = op_string_1
    j, op2 = op_string_2
    if op1 == "A" and op2 == "A":
        return AA(G)[i,j]
    elif op1 == "B" and op2 == "B":
        return BB(G)[i,j]
    elif op1 == "A" and op2 == "B":
        return AB(G)[i,j]
    elif op1 == "B" and op2 == "A":
        return BA(G)[i,j]
    
    ###GENERAL SIGMA X X 
def remove_duplicates_in_pairs(vec):
    unique_vals, counts = np.unique(vec, return_counts=True)
    filtered_vals = unique_vals[counts % 2 != 0]
    return filtered_vals.tolist()
def sigma_general(indices,Gi,L):
    """
    This way is on the order of 99% faster than using a pfaffian method. 
    With Pfaffians the matrix grows very quickly.
    Has been checked against pfaffian with over 500 strings of lengths up to 25 and error was on average 1e-16

    Calculates the expectation values of sigma_x operators put at arbitary sites"
    Inputs:
    indices = list of integers
    Gi = 2L x2L correlation matrix corresponding to quantum state of interest

    Outputs:
    complex scalar 

    Example:

    indices = [1, 2, 3, 4, 5] 
    sigma_5pt = sigma_general(indices,G)
    """
    #Correlation Matrices
    F = Gi[:L,L:]
    G = Gi[:L,:L]
    #M[x,y] = <BxAy>
    M = np.eye(L)- 2*(G+F)
    #Sigma matrices on different sites commute
    indices = np.sort(indices)
    #Remove any duplicates as sigma_x^2 = 1
    indices = remove_duplicates_in_pairs(indices)
    if len(indices)%2 == 1:
        return 0
   #     constant = L//2
    #    indices = list(indices) + [x + constant for x in indices]
     #   return np.sqrt(np.abs(sigma_general(indices,Gi,L)))
    
    #Bs sit on odd sites
    odd_sites = np.array(indices[::2])
    #As site on even sites
    even_sites= np.array(indices[1::2])
    #Get string lengths
    JW_string_lengths = even_sites-odd_sites\
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

    ###Building C
    C = np.zeros( (N,N))
    for nx in range(N):
        for ny in range(N):
            Bx = B_coords[nx]
            Ay = A_coords[ny]

            C[nx,ny] = M[Bx,Ay]

    return la.det(C)

##Depreciated Method##

# def sigma_general(indices,G):
#     """
#     Calculates the expectation values of sigma_x operators put at arbitary sites"
#     Inputs:
#     indices = list of integers
#     G = 2L x2L correlation matrix corresponding to quantum state of interest

#     Outputs:
#     complex scalar 

#     Example:

#     indices = [1, 2, 3, 4, 5] 
#     sigma_5pt = sigma_general(indices,G)
#     """
#     #Commutation lets us sort indices
#     indices = np.sort(indices)
#     #If odd, return 0 (Even projection)
#     if len(indices)%2 == 1:
#         constant = 10
#         indices = list(indices) + [x + constant for x in indices]
#         return np.sqrt(np.abs(sigma_general(indices,G)))
    
#     #Helper function to remove duplicates in the list
#     # <Sigma^2> = 1 for all indices and sigmas
#     def remove_duplicates_in_pairs(vec):
#         unique_vals, counts = np.unique(vec, return_counts=True)
#         filtered_vals = unique_vals[counts % 2 != 0]
#         return filtered_vals.tolist()
#     #Remove duplicate indices
#     indices = remove_duplicates_in_pairs(indices)
#     #<sigma^2> = 1
#     if len(indices)==0:
#         return 1
    
#     #Prepare operator string list
#     string_list = []
#     #Loop over each index
#     for i in range(len(indices)):
#         #All operators are of the form B-A
#         if i%2==0:
#             #B
#             string_list.append([indices[i],"B"])
#         else:
#             #A
#             string_list.append([indices[i],"A"])
    
#     for i in range(0,len(indices),2):
#         #Fill in JW-Strings, these commute and therefore can be added in at the end
#         for a in range(indices[i]+1,indices[i+1]):
#             #Adds "A B " terms
#             string_list.append([a,"A"])
#             string_list.append([a,"B"])
    
#     M = construct_M_matrix(string_list, correlation_function, G)


#     return np.sqrt(la.det(M))


#Projectors
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
from collections import Counter
def P_n(n,G,L):
    """
    P_n scales directly with the number of terms since even small odd sigma correlations need larger support to calcualte.
    Taking out odd terms works to make easier, but scaling is still 2^n.
    Minor improvements can still be had though.
    """
    #For most cases we use periodic boundary conditions
    #That said, watch out for this definition of indices
    indices = [i for i in range(0,n)]
    terms = all_combinations(indices)
    dat = []

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

    for term in range(len(vecs)):
        dat.append(sigma_general(vecs[term],G,L)*counts[term])
    dat.append(sigma_general([],G,L))

    #All terms have equal weight. 
    return np.sum(dat)/2**n


#### MODEL

def H_bdg(h, L, J, boundary_condition = "ABC"):
    #Generating Hamiltonian
    A = np.zeros((L, L))
    B = np.zeros((L, L))

    for j in range(L):
        A[j, j] = 2 * h
    for j in range(L - 1):
        A[j, j + 1] = A[j + 1, j] = -J
        B[j, j + 1] = -J
        B[j + 1, j] = J
    


   # Antiperiodic boundary condition (ABC) introduces a minus sign
    if boundary_condition== "ABC":
        A[0, L - 1] = A[L - 1, 0] = J
        B[L - 1, 0] = J
        B[0, L - 1] = -J

    elif boundary_condition== "PBC":
        A[0, L - 1] = A[L - 1, 0] = -J
        B[L - 1, 0] = -J
        B[0, L - 1] = J

    return 1/2*np.block([[A, B], [-B, -A]])


from scipy.integrate import solve_ivp

def TFIM_time_evolve(N_steps,tau, h0,hf, J, L, U0= None, bc = "ABC"):
    
    def rhs(t, U_flat, h_t, tau, h0, hf, J,bc):
        h = h_t(t,tau,h0,hf)
        H = H_bdg(h,L,J)
        U = U_flat.reshape(2*L,2*L)
        dUdt = -1j*2*H@U
        return dUdt.flatten()
    
    #No input vector then grab initial ground state
    if U0 == None:
        H0 = H_bdg(h0,L,J,boundary_condition= bc)
        _,U0 = la.eigh(H0)
    #ABSOLUTELY NECESSARY STEP
    print(type(U0))
    U0 = U0.astype(np.complex128) 
    #solver stuff 
    args = (h_t,tau,h0,hf,J,bc)
    t_grid = np.linspace(0,tau,N_steps)
    sol = solve_ivp(rhs, [0, tau], U0.flatten(), args=args, t_eval=t_grid, method="RK45", vectorized=False, dtype=np.complex128)
    
    
    U_t = np.array([sol.y[:, i].reshape(2*L, 2*L) for i in range(len(sol.t))])
    #print(sol.message)
    return U_t