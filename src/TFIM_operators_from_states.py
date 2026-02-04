
################################################
#                   IMPORTS
################################################
from project_imports import * 
################################################
#           HELPFUL FUNCTIONS
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


def Pfaffian(M):
    """
    Compute the Pfaffian of an antisymmetric matrix M.
    M must be (2N x 2N) and antisymmetric.
    ChatGPT wrote this....
    """
    M = np.array(M, dtype=complex)
    n = M.shape[0]

    if n % 2 != 0:
        return 0.0

    pf = 1.0 + 0j
    A = M.copy()

    for k in range(0, n-1, 2):
        # find a nonzero pivot
        if abs(A[k, k+1]) < 1e-14:
            swap_found = False
            for l in range(k+2, n):
                if abs(A[k, l]) > 1e-14:
                    # swap columns k+1 and l
                    A[:, [k+1, l]] = A[:, [l, k+1]]
                    A[[k+1, l], :] = A[[l, k+1], :]
                    pf *= -1
                    swap_found = True
                    break
            if not swap_found:
                return 0.0

        pf *= A[k, k+1]

        for i in range(k+2, n):
            for j in range(i+1, n):
                A[i, j] -= (
                    A[k, i]*A[k+1, j] - A[k, j]*A[k+1, i]
                ) / A[k, k+1]
                A[j, i] = -A[i, j]

    return pf


################################################
#           TFIM FUNCTIONS
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

def GA(r,state):
    """r = integer
    state = [us, vs, k]
    calculates {ai,aj} correlation of r=j-i
    """
    u = state[0]
    v = state[1]
    k = state[2]
    #Equation 7.2.26a
    Gr = 2*cos(k*r)*(np.abs(u)**2+np.abs(v)**2)+2j*sin(k*r)*(np.conj(u)*v+np.conj(v)*u)
    return -np.mean(Gr)/2

def GB(r,state):
    """r = integer
    state = [us, vs, k]
    calculates {bi,bj} correlation of r=j-i
    """
    u = state[0]
    v = state[1]
    k = state[2]
    #Equation 7.2.26a
    Gr = 2*cos(k*r)*(np.abs(u)**2+np.abs(v)**2)-2j*sin(k*r)*(np.conj(u)*v+np.conj(v)*u)
    return -np.mean(Gr)/2
################################################
#           Majorna Correlators
################################################
def Majorana_String_single_pair_X(pair):
    """
    Helper Function
    pair = [i, j] with i < j
    Returns the Majorana string as a list of ("A"/"B", site) tuples.
    For X projectors, any two sites are connected with a string of majorana fermions. 
   < B_i \sum (Aj Bj) A_f>
    """
    i, j = pair
    #Start the string
    string = [("B", i)]
    for s in range(i + 1, j):
        #Add z terms
        string.append(("A", s))
        string.append(("B", s))
    #End string
    string.append(("A", j))
    return string
def Create_Majorana_Strings_Z(indices):
    """Helper Function
    For Z projectors, any site simply contributes a pair of Majorna Fermions
    <(B0A0)(B1A1)...>
    """
    string = []
    for i in indices:
        #Add z terms
        string.append(("B", i))
        string.append(("A", i))
    return string

def Create_Majorana_Strings_X(pairs):
    """Helper Function
   Connects pairs of indices to create a full majorana string for sigma x operators.
    """
    lst=[]
    for pair in pairs:
        lst+= Majorana_String_single_pair_X(pair)
    return lst

def Majorana_Correlation(a,b,state):
    """Computes correlation of any two Majorana fermions
    Fermions must be in the form a= ("A/B",index)
    """
    #Unpack Variables
    S1,i = a
    S2,j = b

    if S1=="A" and S2=="A":
        #This minus sign proves to be right though it bothers me... check <sigma_0 sigma_j> in a quench to hf=0 to verify.
        return -GA(j-i,state)
    elif S1=="A" and S2=="B":
        return -G(i-j,state)
    elif S1=="B" and S2=="A":
        return G(j-i,state)
    else:
        return GB(j-i,state)

def Majorana_Covariance(Majorana_String,state):
    """
    Creates a covariance matrix out of a string of Majorna operators, the pfaffian of which is expectation value of the string using wick contractions.
    If you don't use the pfaffian, there are actually some problems with using sqrt(det(M))!
    """
    N = len(Majorana_String)
    M = np.zeros((N,N),dtype=complex)
    for a in range(N):
        for b in range(N):
            if a < b:
                M1,M2= Majorana_String[a], Majorana_String[b]
                M[a, b] = Majorana_Correlation(M1, M2,state)
                #Matrix is anti-symmetric
                M[b, a] = -M[a, b]
    return M 

################################################
#           X Operators
################################################

def sx_correlation(indices,state):
    """
    Computes the correlation of a string of sigma_x operators at sites labeled by indices. 
    This code commutes them into an ordered form (sigma_i<sigma_j<..)
    Then evaluates the expectation through the use of Majorana Fermions
    """
    #Remove Duplicates, sort (allowed via commutation algebra and sigma^2=1)
    indices = list(set(indices))
    #If there is an odd number of sites, return 0.
    if len(indices)%2==1:
        return 0
    #Seperate indices into pairs
    pairs = [indices[i:i+2] for i in range(0, len(indices), 2)]
    #Create full Majorna String
    Majorana_String = Create_Majorana_Strings_X(pairs)
    #Create full Covariance Matrix
    Majorana_cov = Majorana_Covariance(Majorana_String,state)
    #Evaluate expectation with Pfaffian
    return Pfaffian(Majorana_cov)

def Px_n(n,state):
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
      
        a = array(a)
        x_c.append(tuple(a - a[0]))  # convert to tuple for hashing
    counter = Counter(x_c)

    # Extract vecs and counts as separate arrays/lists
    vecs = list(counter.keys())
    counts = array(list(counter.values()))

    ##
    dat = []

    for term in range(len(vecs)):
        indices = vecs[term]
        degen = counts[term]

        if len(indices)%2 == 1:
            dat.append(0)
        else:
            dat.append(sx_correlation(indices,state)*degen)
    #Manually add the 1 
    dat.append(1)
    #All terms have equal weight. 
    return np.sum(dat)/2**n



def Px_n_correlations(n,l,state, even = True):
    """<P_n P_n(r)>"""
    #Term 1 
    indices = [i for i in range(0,n)]
    #Term 2 
    indices +=[i for i in range(n+l-1,2*n+l-1)]
    terms = all_combinations(indices)
   # print(terms)
     ###
    x_c = []
    for a in terms[1:]:
      
        a = array(a)
        x_c.append(tuple(a - a[0]))  # convert to tuple for hashing
    counter = Counter(x_c)

    # Extract vecs and counts as separate arrays/lists
    vecs = list(counter.keys())
    counts = array(list(counter.values()))
    #
    dat = []
    for term in range(len(vecs)):
        indices = vecs[term]
        degen = counts[term]

        if len(indices)%2 == 1:
            dat.append(0)
       
        else:
            dat.append(sx_correlation(indices,state)*degen)
    dat.append(1)
    #All terms have equal weight. 

    return np.sum(dat)/2**(2*n)#,terms


################################################
#           Z Operators
################################################
def sz_correlation(indices,state):
    """
    Computes the correlation of a string of sigma_x operators at sites labeled by indices. 
    This code commutes them into an ordered form (sigma_i<sigma_j<..)
    Then evaluates the expectation through the use of Majorana Fermions
    """
    #Remove Duplicates, sort (allowed via commutation algebra and sigma^2=1)
    indices = list(set(indices))
    #Create full Majorna String
    Majorana_String = Create_Majorana_Strings_Z(indices)
    #Create full Covariance Matrix
    Majorana_cov = Majorana_Covariance(Majorana_String,state)
    #Evaluate expectation with Pfaffian
    return Pfaffian(Majorana_cov)
    
def Pz_n(n,state):
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
        a = array(a)
        x_c.append(tuple(a - a[0]))  # convert to tuple for hashing
    

    counter = Counter(x_c)

    # Extract vecs and counts as separate arrays/lists
    vecs = list(counter.keys())
    counts = array(list(counter.values()))

    ##
    dat = []

    for term in range(len(vecs)):
        indices = vecs[term]
        degen = counts[term]

        dat.append(sz_correlation(indices,state)*degen)
    dat.append(1)
    #All terms have equal weight. 
    return np.sum(dat)/2**n

#P_n correlations do not calculate connected correlator at all
def Pz_n_correlations(n,l,state):
    #Term 1 
    indices = [i for i in range(0,n)]
    #Term 2 
    indices +=[i for i in range(n+l-1,2*n+l-1)]
    terms = all_combinations(indices)
   # print(terms)
     ###
    x_c = []
    for a in terms[1:]:
      
        a = array(a)
        x_c.append(tuple(a - a[0]))  # convert to tuple for hashing
    counter = Counter(x_c)

    # Extract vecs and counts as separate arrays/lists
    vecs = list(counter.keys())
    counts = array(list(counter.values()))
    #
    dat = []
    for term in range(len(vecs)):
        indices = vecs[term]
        degen = counts[term]

        dat.append(sz_correlation(indices,state)*degen)
    dat.append(1)
    #All terms have equal weight. 

    return np.sum(dat)/2**(2*n)#,terms
