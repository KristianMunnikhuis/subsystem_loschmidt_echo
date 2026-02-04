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


################################################
#           FULL SYSTEM SOLUTIONS
################################################
def k_vals(L):
    """Generates k values >0"""
    x=(2 * np.arange(-L//2, L//2) + 1) * np.pi / (L)
    return x[x>0]

def state(L, tau, t):
    """Generates state as [us, vs, k]
    """
    k = k_vals(L)
    us = np.array([complex(u_tilde(tau, q, t)) for q in k], dtype=np.complex128)
    vs = np.array([complex(v_tilde(tau, q, t)) for q in k], dtype=np.complex128)
    k  = np.array(k, dtype=np.float64)
    return [us, vs, k]
