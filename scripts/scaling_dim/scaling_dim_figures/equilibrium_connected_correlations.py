import sys
import os
sys.path.append("../../../src/")
#Imports
import single_particle_sector as sps
import kblock_ising_model as kb
import numpy as np
import scipy as sp
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, scatter
from collections import Counter
from tqdm import tqdm

##Personal Imports##

from numpy import pi, cos, sin, exp, log, sqrt, array
#import equilibrium_critical_correlations as ecc
from equilibrium_critical_correlations import P_n, sigma_general, P_n_correlations, Dn

###correlator
def correlator_data(n,rmin,rmax):
    Pn = P_n(n)

    r_points = np.arange(rmin,rmax)
    ##Definition of connected correlator in a translationally invariant system
    data = array([P_n_correlations(n,r)-Pn**2 for r in tqdm(r_points)])
    return [data,r_points]

def correlator_scaled(n,rmin,rmax):
    r_points = np.arange(rmin,rmax)
    Pn = P_n(n)
    ##Definition of connected correlator in a translationally invariant system
    data = array([(P_n_correlations(n,r)-Pn**2)/Pn**2 for r in tqdm(r_points)])
    return [data,r_points]

def correlator_unconnected(n,rmin,rmax):
    r_points = np.arange(rmin,rmax)
   # Pn = P_n(n)
    ##Definition of connected correlator in a translationally invariant system
    data = array([(P_n_correlations(n,r)) for r in tqdm(r_points)])
    return [data,r_points]



if len(sys.argv) != 2:
    print("Usage: python script.py <n>")
    sys.exit(1)

try:
    n = int(sys.argv[1])

except ValueError:
    print("n must be an integer")
    sys.exit(1)

# Do whatever you want with n here
print(f"n = {n}")


folder_name = "connected_correlations"

# Create the folder if it doesn't exist
if not os.path.exists(folder_name):
    os.makedirs(folder_name)





rmin=0

rmax=200



    # get only the data part
data, _ = correlator_data(n, rmin, rmax)

    # save it
np.save(f"{folder_name}/corr_connected_n{n:02d}_r{rmin:02d}-{rmax:02d}.npy", data)
    # or text:
    # np.savetxt(f"data/correlator_n{n}_r{rmin}-{rmax}.txt", data)
print(f"n = {n} completed!")

