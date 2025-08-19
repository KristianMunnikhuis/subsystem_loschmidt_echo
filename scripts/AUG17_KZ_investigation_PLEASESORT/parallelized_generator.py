import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot
from tqdm import tqdm
import sys, os, time
sys.path.append("../../src/")
import kibble_zurek as kz
from concurrent.futures import ProcessPoolExecutor

# PARAMETERS
L = 1000
n_max = 24
tau_lower_exponent = -2.5
tau_upper_exponent = 2.5
N_taus = 50

n_values = np.array([ni for ni in range(2, n_max)])
taus = np.logspace(tau_lower_exponent, tau_upper_exponent, N_taus)
N_n = len(n_values)

# Preallocate arrays
KZ_array = np.zeros((N_n, N_taus), dtype=float)

def _worker(args):
    tau, L, n_values = args
    t = tau
    state_t = kz.state(L, tau, t)
    return np.array([kz.P_n(int(n), state_t) for n in n_values], dtype=float)

if __name__ == "__main__":
    start_time = time.time()
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        rows = list(tqdm(ex.map(_worker, [(tau, L, n_values) for tau in taus]),
                         total=len(taus)))
    KZ_array[...] = np.stack(rows, axis=1)  # shape: (len(n_values), len(taus))
    print(f"Time total = {time.time()-start_time}")
    os.makedirs("Data", exist_ok=True)
    np.savez(f"Data/KZ_data_nmax={n_max}_NTaus=_{N_taus}.npz",
             KZ_Array=KZ_array, taus=taus, n_values=n_values)
