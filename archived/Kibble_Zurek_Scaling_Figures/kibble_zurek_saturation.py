import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot
from tqdm import tqdm
import sys
sys.path.append("../../../src/")
import kibble_zurek as kz
from time import time

#PARAMETERS
L = 1000
n_max = 16
tau_lower_exponent= -2.5
tau_upper_exponent= 2.5
N_taus= 50

#16 may be wayyyyyy tooo high
n_values = np.array([ni for ni in range(2,n_max)])
taus = np.logspace(tau_lower_exponent,tau_upper_exponent,N_taus)
N_n = len(n_values)

# Preallocate arrays
KZ_array = np.zeros((N_n,N_taus))


if __name__== "__main__":
    start_time = time()
    for t_i, tau in enumerate(tqdm(taus)):
        t = tau
        state_t = kz.state(L, tau, t)
        for i, n in enumerate(n_values):
            P_kz = kz.P_n(n, state_t)
            KZ_array[i,t_i] = P_kz
    print(f"Time total = {time()-start_time}")
    np.savez(f"Data/KZ_data_nmax={n_max}_NTaus=_{N_taus}.npz",KZ_Array=KZ_array,taus=taus,n_values=n_values)

