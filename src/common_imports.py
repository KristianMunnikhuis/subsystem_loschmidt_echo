import numpy as np
from numpy import cos, sin, exp, log, pi, cosh, tanh, sinh, array, linspace, logspace,sqrt, conj, vdot
from numpy import linalg as la
from numpy.linalg import eigh 
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, scatter, xlabel,ylabel,savefig, xlim, ylim, yscale, xscale, loglog, legend,show,title, imshow
import scipy
import os 
from tqdm import tqdm
from numpy.fft import fft, fftshift,fftfreq
from joblib import Parallel, delayed
from time import time
from scipy.integrate import solve_ivp