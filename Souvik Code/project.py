import os
import eigen
import multiprocess as mp
from tqdm import tqdm
from quspin.tools.lanczos import lanczos_full,expm_lanczos
from quspin.operators import hamiltonian
from quspin.basis import boson_basis_1d
import matplotlib.pyplot as plt
import numpy as np 

os.environ['OMP_NUM_THREADS'] = '22'
os.environ['OPENBLAS_NUM_THREADS'] = '22'
os.environ['MKL_NUM_THREADS'] = '22'

L = 14
N = L 
nb = 0.5
sps = 3
#int(nb*L+1)

basis = boson_basis_1d(L,nb=nb,sps=sps)
D = basis.Ns


print("Basis dimension: ",basis.Ns,"\n")

no_checks = dict(check_herm=False,check_symm=False,check_pcon=False)


def full(psi,psi_i):

    catch_out = []
    batch_size = 100
    batch_states = np.array_split(psi,batch_size,axis=1)

    for j in tqdm(range(batch_size)):

        psi_batch = np.array(batch_states[j])

        
        if j==0:
            catch_out = list(np.absolute(np.einsum('kl,l->k',np.conjugate(psi_batch.T),psi_i))**2)
        else:
            catch_out.extend(np.absolute(np.einsum('kl,l->k',np.conjugate(psi_batch.T),psi_i))**2)

    return np.asarray(catch_out)




def calc_prob(psi,ns,string):

    def site_project(psi_r,rung,p):

        norm = 1
        for alpha in range(0,sps):

            if alpha!=p:

                pr_ir = hamiltonian([["n",[[-1,rung]]],["I",[[alpha,rung]]]],[],basis=basis,**no_checks)
                psi_r = pr_ir.dot(psi_r)
                norm = norm*(alpha-p)

        return psi_r/norm


    def site_average(i,queue,psi_batch):

        init_str = np.array([char for char in string])
        
        if i+ns <= L:

            n = 0
            rung_psi_r = np.copy(psi_batch)

            while  n<=ns-1:

                rung_psi_r = site_project(rung_psi_r,(i+n),int(init_str[(i+n)]))
                n+=1

            queue.put(np.einsum('kl,lk->k',np.conjugate(psi_batch.T),rung_psi_r).real)


    sites = np.arange(0,N,1)
    avlst,catch_out = [],[]
    batch_size = 100
    batch_states = np.array_split(psi,batch_size,axis=1)
    
    processes = []
    queue = mp.Queue()

    for j in tqdm(range(batch_size)):
        
        for i in range(N):
            process = mp.Process(target=site_average, args=(i,queue,np.array(batch_states[j])))
            processes.append(process)
            process.start()

        for process in processes: process.join()    


        while not queue.empty():
            catch_out.append(queue.get())
            
        processes = []
        

        if j == 0:
            avlst = list(np.mean(np.array(catch_out),axis=0))
        else:
            avlst.extend(np.mean(np.array(catch_out),axis=0))
        
        catch_out = []

    return avlst

############################ Hamiltonian ##############################

U = 1.0
V = 0.0 
J_par = 1.0
J_perp = 0.0
delta = 0.0

int_list_0_p = [[0.5*delta*(-1)**i,i] for i in range(N)]
int_list_1_p = [[-0.5*U,i] for i in range(N)]
int_list_2_p = [[0.5*U,i,i] for i in range(N)]

v_lst = [[V,i,(i+1)] for i in range(N-1)]
hop_list_p = [[-J_par,i,(i+1)] for i in range(N-1)]
hop_list_hc_p = [[J.conjugate(),i,j] for J,i,j in hop_list_p]
hop_list_per = [[-J_perp,i,(i+2)] for i in range(N-2)]
hop_list_hc_per = [[J.conjugate(),i,j] for J,i,j in hop_list_per]


static_p = [
            ["+-",hop_list_p],
            ["-+",hop_list_hc_p],
            ["nn",int_list_2_p],
            ["nn",v_lst],
            ["n",int_list_1_p],
            ["n",int_list_0_p]
           ]

dynamic = []

H = hamiltonian(static_p,[],basis=basis,**no_checks,dtype=np.float32)

print("\nStarting diagonalization\n")
F,V = eigen.eigh(solve = False, filename="projector.hdf5")
print("\nDiagonalization complete\n")


dEtab = [F[i+1]-F[i] for i in range(D-1)]
offset = 0
tol = 1e-14
for i in range(len(dEtab)-1):
    if dEtab[i-offset]<tol:
        del(dEtab[i-offset])

    offset+=1


meanlevspace = np.mean(dEtab)
sigma = meanlevspace

########################################################################

def dos(e):
    return np.sum(np.exp(-((F-e)**2)/(2*sigma**2))*(1/(np.sqrt(2*np.pi)*sigma)))
dosv = np.vectorize(dos)


def calibrate(string,n):

    global energy 
    psi_i = np.zeros(basis.Ns,dtype=np.float32)
    psi_i[basis.index(string)] = 1.0

    energy = H.expt_value(psi_i)
    print("\n","Energy density = ",energy/L,"\n")

    omega = dosv(F)

    with open("dos_L"+str(L)+".out",'w') as f:
        np.savetxt(f, omega, delimiter=',')

#    plt.semilogy(F,omega,'o')
#    plt.axvline(x=energy,linestyle='dashed',color='black',zorder=10)
#    plt.show()

    if n<L:

        p = calc_prob(V,n,string)

    elif n == L:

        p = full(V,psi_i)

    return F/L,p



def time_average():
       
    state_str = '1100'

    msr_str1 = state_str

    psi_i = np.zeros(basis.Ns,dtype=np.float32)
    psi_i[basis.index(state_str)] = 1.0

    energy = H.expt_value(psi_i)
    print("\n","Energy= ",energy,"\n")

    print("\nStarting lanczos:\n")
    
    dt = 0.25
    kdim = 20
    nsteps = 400

    psit_la = psi_i

    psit = np.zeros((basis.Ns,nsteps),dtype=np.float32)

    for i in range(nsteps):
        
        if i==0:
            psit[:,i] = psi_i
            
        else:

            E2,V2,Q2 = lanczos_full(H,psit_la,kdim)
            psit_la = expm_lanczos(E2,V2,Q2,a=-1j*dt)
            psit[:,i] = psit_la

        print("finished evolution till: {0:f}.".format((i+1)*dt))

    print("\nStarting expectation calculation:\n")
    
    mean2s = calc_prob(psit,2,msr_str1)
    mean4s = calc_prob(psit,4,msr_str1)
    mean6s = calc_prob(psit,6,msr_str1)
    mean8s = calc_prob(psit,8,msr_str1)

    print("\nfinished delta: {0:f}.\n".format(delta))

    return mean2s,mean4s,mean6s,mean8s


"""
####################################### time average #################################
p2s,p4s,p6s,p8s = time_average()

with open("entropy2_"+str(L)+"_2.out",'w') as f:
    np.savetxt(f, [p2s], delimiter=',')
with open("entropy4_"+str(L)+"_2.out",'w') as f:
    np.savetxt(f, [p4s], delimiter=',')
with open("entropy6_"+str(L)+"_2.out",'w') as f:
    np.savetxt(f, [p6s], delimiter=',')
with open("entropy8_"+str(L)+"_2.out",'w') as f:
    np.savetxt(f, [p8s], delimiter=',')
#######################################################################################
"""

############################################################# full projector ED #####################

state_str = '01010101010101'
e,p = calibrate(state_str,L)

######################################################################################################



################################## partial projectors ################################################

e,p1 = calibrate(state_str,1)
print("\nCompleted N=",1,"\n")
e,p2 = calibrate(state_str,2)
print("\nCompleted N=",2,"\n")
e,p3 = calibrate(state_str,3)
print("\nCompleted N=",3,"\n")
e,p4 = calibrate(state_str,4)
print("\nCompleted N=",4,"\n")
e,p5 = calibrate(state_str,5)
print("\nCompleted N=",5,"\n")
e,p6 = calibrate(state_str,6)
print("\nCompleted N=",6,"\n")
e,p7 = calibrate(state_str,7)
print("\nCompleted N=",7,"\n")
e,p8 = calibrate(state_str,8)
print("\nCompleted N=",8,"\n")
e,p9 = calibrate(state_str,9)
print("\nCompleted N=",9,"\n")
e,p10 = calibrate(state_str,10)
print("\nCompleted N=",10,"\n")
e,p11 = calibrate(state_str,11)
print("\nCompleted N=",11,"\n")
e,p12 = calibrate(state_str,12)
print("\nCompleted N=",12,"\n")
e,p13 = calibrate(state_str,13)
print("\nCompleted N=",13,"\n")


#plt.plot(e,p1,'o')
#plt.plot(e,p3,'o')
#plt.plot(e,p5,'o')
#plt.plot(e,p7,'o')
#plt.plot(e,p9,'o')

#plt.axvline(x=energy/L,linestyle='dashed',color='black',zorder=10)
#plt.show()

with open("energy_mid_L14.out",'w') as f:
    np.savetxt(f, e, delimiter=',')

with open("p_1_mid_L14.out",'w') as f:
    np.savetxt(f, p1, delimiter=',')
with open("p_2_mid_L14.out",'w') as f:
    np.savetxt(f, p2, delimiter=',')
with open("p_3_mid_L14.out",'w') as f:
    np.savetxt(f, p3, delimiter=',')
with open("p_4_mid_L14.out",'w') as f:
    np.savetxt(f, p4, delimiter=',')
with open("p_5_mid_L14.out",'w') as f:
    np.savetxt(f, p5, delimiter=',')
with open("p_6_mid_L14.out",'w') as f:
    np.savetxt(f, p6, delimiter=',')
with open("p_7_mid_L14.out",'w') as f:
    np.savetxt(f, p7, delimiter=',')
with open("p_8_mid_L14.out",'w') as f:
    np.savetxt(f, p8, delimiter=',')
with open("p_9_mid_L14.out",'w') as f:
    np.savetxt(f, p9, delimiter=',')
with open("p_10_mid_L14.out",'w') as f:
    np.savetxt(f, p10, delimiter=',')
with open("p_11_mid_L14.out",'w') as f:
    np.savetxt(f, p11, delimiter=',')
with open("p_12_mid_L14.out",'w') as f:
    np.savetxt(f, p12, delimiter=',')
with open("p_13_mid_L14.out",'w') as f:
    np.savetxt(f, p13, delimiter=',')
with open("p_14_mid_L14.out",'w') as f:
    np.savetxt(f, p, delimiter=',')


#######################################################################################################


