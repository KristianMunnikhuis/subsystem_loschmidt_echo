import quspin
from quspin.operators import hamiltonian,quantum_operator
from quspin.basis import spin_basis_1d
from project_imports import *
from itertools import combinations
from functools import partial
def projector_string(ns,s,b): #ns = length of string s = operator
    #Max string length
    L = b.L
    if ns>L:
        ns=L
    else:
        pass
    
    
    def linears(start,m):
        #Find all indices in string
        ind=[i%L for i in range(start,start+ns)]
     #All combinations
        perm=list(combinations(ind,m))
      #Add each coupling to operator
        for j in range(len(perm)):
#            S_temp=[((-1)**m)*(2**m)/(L*(2**ns))] ########### Down #
            S_temp=[1/(1*(2**ns))] ###################### Up
            for i in range(len(perm[j])):
                S_temp.append(perm[j][i])
        ##Saving Couplin gs lists
            if j==0: 
                S_temp_2=[S_temp]
            else: 
                S_temp_2.append(S_temp)
        return S_temp_2
    #List to store orders of interaction
    o_type=np.empty(ns,dtype=object)
    
    for i in range(ns):
        #ex: "z"*2= "zz"
        o_type[i]=s*(i+1)

    #Iterate over each Site
    for k in range(1):
        #Iterate over different string lengths
      for i in range(ns):
          #First item defines list
          if k==0 and i==0:
              #example: i = 0 operator_list = [['x', linears(k,1)]]
              operator_list=[[o_type[i],linears(k,i+1)]]
          else:
              operator_list.append([o_type[i],linears(k,i+1)])
    S_temp_1 = [[(1/((ns*2**ns))), i] for i in range(ns)]
    #S_temp_1 = [[(1/((ns*2**(ns)))), 1] ]
    operator_list.append(["I",S_temp_1])
    operator_dict=dict(S=operator_list)
    O = quantum_operator(operator_dict,basis=b,check_symm=False,check_herm=False)
    #print(operator_list)
    return array(O.todense())


def TFIM(basis):
    L = basis.L
    NN = [[-1,i,(i+1)%L] for i in range(L)]
    single_site = [ [-1, i] for i in range(L)]

    ferro_terms = [["xx",NN]]
    para_terms = [["z",single_site]]
    op_dict = dict(H0=ferro_terms,g=para_terms)

    return quantum_operator(op_dict,basis=basis,check_herm=False,check_symm=False)