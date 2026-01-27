import sys
from common_imports import *
from scipy.integrate import solve_ivp
from numpy import conj
from numpy import linalg as la
import kibble_zurek as kz
from time import time
from joblib import Parallel, delayed
import numpy as np
import time_evolved_z_correlators as zz
#Load Pauli Matrices
s3,s1,s2 = array([[1,0],[0,-1]],dtype=complex),array([[0, 1],[1,0]],dtype=complex), array([[0,-1j],[1j,0]],dtype=complex)

#Quench Functions
#These must be in the form (t, args) to work! 
def g_linear_hold(t,args):
    """
    args = tau,tf
    A linear quench of the form 1-t/tau
    when t>tf:
    1-tf/tau
    """
    tau,tf= args
    if t<tf:
        return 1-t/tau
    else:
        return 1-tf/tau+(t-tf)/tau
def g_pulse(t,args):
    """
    args= t0,tf,g_off,g_on
    Quenches from g_off to g_on, only while t0<t<tf
    """
    t0,tf,g_off,g_on = args
    if t> t0 and t<tf:
        return g_on
    else:
        return g_off
def g_exp_slowdown(t,args):
    """
    args=A,B,C
    Returns A+B exp(-t*C)
    """
    A,B,C = args
    return A+B*exp(-t*C)

def sin_g(t,args):
    """
    args=A,B,C
    returns A+ B*sin(C*t)
    """
    A,B,C = args
    return  A+B*sin(C*t)


#TFIM functions
def a(q,g):
    return g+ cos(q)
def b(q,g):
    return sin(q)
def w(q,g):
    return sqrt(g*g+2*g*cos(q)+1)
def BdG(q,g):
    return -2*(a(q,g)*s3+b(q,g)*s2)
###

def BdG_evolve(t,U,g,q,args):
    return -1j*BdG(q,g(t,args))@U
def groundstate(g,L):
    n = np.arange(0, L//2 )
    k = (2 * n + 1) * pi / L
    us,vs=[],[]
    for ki in k:
        V0 = la.eigh(BdG(ki, g))[1][:, 1]
        us.append(V0[0])
        vs.append(V0[1])
    us,vs = array(us),array(vs)
    return [us,vs,k]


def evolve_mode(q,g,t_eval,u0,v0,args):
    V0 = array([u0,v0])
    sol = solve_ivp(
        BdG_evolve,
        [t_eval[0], t_eval[-1]],
        V0,
        args=(g,q,args),
        t_eval=t_eval,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
        max_step=0.05,
    )
    return sol.y[0], sol.y[1]

def numerical_states(L,t_eval,g,args,state_0 =None):
    """Computes a states list which is compatable with code for analytic quenches.
    
    L = System Size
    T_eval = array of time to be evolved over
    g= function for the quench parameter
    args= arguments for the quench parameter
    state_0 = optional input state object.

    Initaliazes state in g(T_eval[0]) if not otherwise given.
    """

    g0 = g(t_eval[0],args)
    print(g0)
    if state_0 ==None:
        u,v,k = groundstate(g0,L)
    else:
        u,v,k = state_0
    
    results = Parallel(n_jobs=-1, prefer="processes")(
        delayed(evolve_mode)(k[i], g, t_eval,u[i],v[i],args) for i in range(len(k))
    )

    us, vs = map(np.array, zip(*results))

    states = [[us[:, i], vs[:, i], k] for i in range(len(t_eval))]
    return states

def overlap(states,state):
    overlaps = []
    u0,v0,_ = state
    for ui,vi,_ in states:
        overlaps.append( np.prod(np.abs(np.conj(u0)*ui+np.conj(v0)*vi))**2)
    return array(overlaps)
