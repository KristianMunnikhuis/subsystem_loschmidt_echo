from itertools import combinations
from collections import Counter
import sys
sys.path.append("../../src")
from common_imports import *
from kibble_zurek import state
from kibble_zurek import sigma_general as sigmax
from kibble_zurek import P_n as Px
from time_evolved_z_correlators import sigma_general as sigmaz
from time_evolved_z_correlators import P_n as Pz
from joblib import Parallel, delayed
import numpy as np

def BA(r, taus):
    def delta(a,b):
        if a==b:
            return 1
        else: 
            return 0
    def sign(a):
        if a ==0:
            return 0
        else:
            return a/abs(a)
    
    ln = log(taus)
    xi = sqrt(taus)
    term1 = -sqrt(2)*exp(-r^2/8/Pi/taus)/pi/xi

    term2 = delta(1,r)+sign(r)*delta(1,r)

    term3 = -2 *exp(
        (-pi *r^2+1j*r*r*ln -8j*taus*taus*ln*ln )/(4*taus*ln*ln)
        )
    term4 = sqrt(1-exp(
        -pi * r* r
                       )/(
                        4*taus*ln*ln
                       )
                       
                       )*sign(r)
    term5 = term3*term4/sqrt(pi*taus*ln)
    
    return term1+term2+term5