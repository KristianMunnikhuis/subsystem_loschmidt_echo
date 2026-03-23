import sys 
sys.path.append("../../src/")
import TFIM_operators_from_states as operators
import numerical_quenches as nq
from project_imports import *

def op(n):
    """Constructs general Projector term"""
    return ["I", f"x{n}", f"z{n}"]

def distribute(L1, L2):
    """Distributes across two lists of strings
    E.g. P2 = (1+sx+sz)(1+sx+sz)=(9 terms)"""
    return [f"{u}{v}" for u in L1 for v in L2]
def distribute_list(List):
    """Distributes across an entire list."""
    L1 = List[0]
    for L2 in List[1:]:
        L1 = distribute(L1,L2)
    return L1

def reduce_string(str):
    """The Identity is removed to make our other problems expectations
    We also work in the even parity sector where <x>=0
    """
    str= str.replace("I","")
    if str.count("x")%2==1:
        #Even Parity Condition
        return ""
    return str


from re import findall
from collections import Counter

def translational_invariance(string):
    """
    Incorporates Translational invariance as a way to reduce the total number of calculations we need to make.
    ex: "x2z3"->"x0z1"
    """
    parts = findall(r'([A-Za-z]+)(\d+)', string)

    if not parts:
        return string

    x0 = int(parts[0][1])

    shifted = [(op, str(int(idx) - x0)) for op, idx in parts]

    return "".join(op + idx for op, idx in shifted)

def filter_terms(list):
    """Given a list of terms, removes any identities and any odd-partiy x terms"""
    new_list = [reduce_string(item) for item in list if reduce_string(item) != ""]
    #Translational_Invariance
    new_list = [translational_invariance(item) for item in new_list]
    
    return new_list

def degeneracy_counter(terms):
    """Counts the Degeneracy terms and seperates them, to reduce total number of computations."""
    items =list(Counter(terms).keys())
    degeneracies =list(Counter(terms).values())
    return items,degeneracies
###############################
#Majorana Fermions
###############################

def Prepare_Term_For_Majorana_Strings(term):
    """
    Prepares a string to be seperated into its indices, so that we can calculate Majorna Strings.
    """
    term = list(term)
    operators=term[::2]
    indices = term[1::2]
    z_indices=[]
    x_indices=[]
    for operator, index in zip(operators,indices):
        if operator=="z":
            z_indices.append(int(index))
        elif operator=="x":
            x_indices.append(int(index))
    return z_indices, x_indices

from TFIM_operators_from_states import Create_Majorana_Strings_Z, Create_Majorana_Strings_X
def Majorana_Strings_X(indices):
    pairs = [indices[i:i+2] for i in range(0, len(indices), 2)]
    Majorana_String = Create_Majorana_Strings_X(pairs)
    return Majorana_String


def Majorana_String_Filter(String):
    """
    We Combine the list of Majornas contributed by Z and X corerelators. where they interact is that when there is an overlap B^2=A^2propto1 as these operators are their own inverses.
    """
    counts = Counter(String)
    new_list = [x for x, c in counts.items() if c == 1]
    return new_list

def Create_Majorana_Strings_XZ(term):
    z_indices, x_indices = Prepare_Term_For_Majorana_Strings(term)
    x_string = Majorana_Strings_X(x_indices)
    z_string = Create_Majorana_Strings_Z(z_indices)
    majorana_string = Majorana_String_Filter(x_string+z_string)

    return majorana_string


def Generic_Correlation(term,state):
    """
    Calculates Majorana string for a term and then runs it through the Majoran Covariance, where we then take a pfaffian. 
    """
    majorana_string = Create_Majorana_Strings_XZ(term)
    M=operators.Majorana_Covariance(majorana_string,state)
    return operators.Pfaffian(M)


#######
#Full Projector
#######

def Unit_Vector_Coefficient(term,nhat):
    xs = term.count("x")
    zs = term.count("z")
    return nhat[0]**xs*nhat[1]**zs
def P_generic(N,nhat,state):
    sites=  np.arange(0,N)
    operators = distribute_list([op(site) for site in sites])
    operators = filter_terms(operators)
    uniques,degeneracies = degeneracy_counter(operators)
    results = Parallel(n_jobs=-1)(
    delayed(compute_term)(unique_term, degen,nhat,state)
    for unique_term, degen in zip(uniques, degeneracies))
    expectation_value = np.sum(results)+1
    return expectation_value/2**(N)



######
#Small Z Expansion
######
def filter_linear_z(terms: list,degeneracies : list,n_order : int):
    """
    Filters out any terms that have a number of z spins not appropriate with the expansion order. then also keeps degeneracy.
    """
    filtered_terms =[]
    filtered_degeneracy=[]
    for term,degen in zip(terms,degeneracies):
        if term.count("z")== n_order:
            filtered_terms.append(term)
            filtered_degeneracy.append(degen)

    return filtered_terms, filtered_degeneracy


def Small_Z_Expansion(N,state,n_order=1):
    sites=  np.arange(0,N)

    operators = distribute_list([op(site) for site in sites])

    operators = filter_terms(operators)

    uniques,degeneracies = degeneracy_counter(operators)

    uniques,degeneracies = filter_linear_z(uniques,degeneracies,n_order)

    results = [Generic_Correlation(term,state)*degen for term,degen in zip(uniques,degeneracies)]
    expectation_value = np.sum(results)
    return (expectation_value)


def compute_term(unique_term, degen,nhat,state):
    coefficient = Unit_Vector_Coefficient(unique_term, nhat)
    O = Generic_Correlation(unique_term, state)
    return (O * degen) * coefficient



    
def opx(n):
    """Constructs general Projector term"""
    return ["I", f"x{n}"]
def opxm(n):
    return ["I",f"-x{n}"]

def Px_n_down(state,n,reverse_indices):
    operators= [opx(ni) if ni not in reverse_indices else opxm(ni) for ni in range(n) ]


    return distribute(operators)

if __name__ == "__main__":
    print(Px_n_down(_,2,[0,1]))