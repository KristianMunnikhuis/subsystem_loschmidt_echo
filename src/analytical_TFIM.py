import numpy as np
import matplotlib.pyplot as plt


#TFIM Functions

#Diagonal Dispersion
def Epsilon_0(k, J, h):
    """
    Calculate the diagonal dispersion relation for the TFIM.
    
    Parameters:
    k (float): The wave vector.
    J (float): The coupling constant.
    h (float): The transverse field strength.
    
    Returns:
    float: The diagonal dispersion relation.
    """
    return 2*h-2*J*np.cos(k)

#Full Dispersion
def Epsilon_h(k, J, h):
    """
    Calculate the full dispersion relation for the TFIM.
    
    Parameters:
    k (float): The wave vector.
    J (float): The coupling constant.
    h (float): The transverse field strength.
    
    Returns:
    float: The full dispersion relation.
    """
    return np.sqrt(Epsilon_0(k,J,h)**2+4*J**2*np.sin(k)**2)


#Diagonalization Angle

def Theta(k, J, h):
    """
    Calculate the diagonalization angle for the TFIM.
    
    Parameters:
    k (float): The wave vector.
    J (float): The coupling constant.
    h (float): The transverse field strength.
    
    Returns:
    float: The diagonalization angle.
    """
 #   return (Epsilon_0(k,J,h)+2j*J*np.sin(k))/Epsilon_h(k,J,h) 
    return np.real(-1j*np.log( (Epsilon_0(k,J,h)+2j*J*np.sin(k))/Epsilon_h(k,J,h) ) )

def Delta(k,h1,h2,J):
    """
    Calculate the difference in angle  for the TFIM.
    
    Parameters:
    k (float): The wave vector.
    h1 (float): The first transverse field strength.
    h2 (float): The second transverse field strength.
    J (float): The coupling constant.
    
    Returns:
    float: The angle differnce
    """
    return Theta(k,J,h2)-Theta(k,J,h1)

#Bogoliouv Amplitudes

def A(k,t,h1,h2,J):
    th = Theta(k,J,h2)
    delta = Delta(k,h1,h2,J)
    eps = Epsilon_h(k,J,h2)
    return np.cos(th/2)*np.cos(delta/2)*np.exp(-1j*eps*t) + np.sin(th/2)*np.sin(delta/2)*np.exp(1j*eps*t)
def B(k,t,h1,h2,J):
    th = Theta(k,J,h2)
    delta = Delta(k,h1,h2,J)
    eps = Epsilon_h(k,J,h2)
    return 1j*np.sin(th/2)*np.cos(delta/2)*np.exp(1j*eps*t) - 1j*np.cos(th/2)*np.sin(delta/2)*np.exp(-1j*eps*t)

k = np.array([0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

np.exp(k)
def aa(l,t):
    phases = np.exp(-1j*l*k)
    return np.sum( A(k,t,h1,h2,J)*B(-k,t,h1,h2,J)*phases )/len(k)
def cc(l,t):
    phases = np.exp(-1j*l*k)
    return np.sum( np.conj(A(k,t,h1,h2,J))*np.conj(B(-k,t,h1,h2,J))*phases )/len(k)

def ca(l,t):
    phases = np.exp(1j*l*k)
    return np.sum( B(k,t,h1,h2,J)*np.conj(B(k,t,h1,h2,J))*phases )/len(k)    

def ac(l,t):
    phases = np.exp(1j*l*k)
    return np.sum( A(k,t,h1,h2,J)*np.conj(A(k,t,h1,h2,J))*phases )/len(k)


h1 = 0.5
h2 = 0.1
J =1 
aa(0,1)