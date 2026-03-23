import sys 

sys.path.append("../../src/")
from project_imports import *
import general_projectors as gp
import re 
def opx(n):
    """Constructs general Projector term"""
    return ["I", f"x{n}"]
def opxm(n):
    return ["I",f"-x{n}"]

def Px_n_down_terms(n,reverse_indices):
    operators= [opx(ni) if ni not in reverse_indices else opxm(ni) for ni in range(n) ]


    return gp.distribute_list(operators)




def filter_parity(terms):
    filtered_terms= []
    for term in terms:
        parity = (-1)**term.count("-")
        term= term.replace("-","")
        if parity == -1:
            term = "-1"+term
        else:
            term = term
        filtered_terms.append(term)
    return filtered_terms

def filter_terms(list):
    """Given a list of terms, removes any identities and any odd-partiy x terms"""
    new_list = [gp.reduce_string(item) for item in list if gp.reduce_string(item) != ""]
    new_list = filter_parity(new_list)
    #Translational_Invariance
    new_list = [translational_invariance(item) for item in new_list]
    
    return new_list



def translational_invariance(string):
    """
    Incorporates translational invariance.
    Example:
        "x2z3"   -> "x0z1"
        "-1x1x2" -> "-1x0x1"
    """

    # --- Preserve leading sign ---
    sign = ""
    if string.startswith("-"):
        sign = "-"
        string = string[1:]   # strip sign for processing

    parts = re.findall(r'([A-Za-z]+)(\d+)', string)

    if not parts:
        return sign + string

    x0 = int(parts[0][1])

    shifted = [(op, str(int(idx) - x0)) for op, idx in parts]

    return sign + "".join(op + idx for op, idx in shifted)


def cancel_opposites(lst):
    counts = Counter(lst)
    result = []

    for s in set(x.lstrip('-') for x in lst):
        pos = counts[s]
        neg = counts['-' + s]
        remaining = abs(pos - neg)

        if pos > neg:
            result.extend([s] * remaining)
        elif neg > pos:
            result.extend(['-' + s] * remaining)

    return result
def compute_strings(u,state):
    terms, degens = u
    summ= 1
    for term, degen in zip(terms,degens):
        if term[0] == "-":
            p = -1 
            term= term.replace("-","")
        else:
            p = 1
        
        summ+= p*degen*gp.Generic_Correlation(term,state)
    return summ

def Px_n_down(n,reverse_indices,state):
    terms = Px_n_down_terms(n,reverse_indices)
    terms = filter_terms(terms)
    terms = cancel_opposites(terms)
    degen_terms = gp.degeneracy_counter(terms)
    expectation = compute_strings(degen_terms,state)
    return expectation/2**n 




if __name__=="__main__":
    from numerical_quenches import groundstate
    gs= [100,1,0]
    for g in gs:
        psi = groundstate(g,100)
        print(f"g={g}, P_UpPdown = {Px_n_down(psi,2,[0])}")