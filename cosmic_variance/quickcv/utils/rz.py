

# Routines for computing the Alcock-Paczynski distortion

# January 2023, Christian Kragh Jespersen, rewritten from IDL


#  H(z) in units where H(0) = 1

import numpy as np 
from scipy.integrate import quad

OmegaM = 0.308
OmegaL = 0.692
OmegaQ = 0.0
wQ = 0.0
eps = 1e-10

def hubble(z):
    x = 1. + z
    OmegaK = 1. - OmegaM - OmegaL - OmegaQ
    h2 = OmegaM * x**3 + OmegaK * x**2 + OmegaL + OmegaQ * x**(3.+3.*wQ)
    return np.sqrt(h2)

# Compute d(chi)/dz where d(chi)^2 = dr^2 / (1-kr^2)

def dchi_dz(z):
    return 1./hubble(z)

# Compute coordinate distance r(z) for FRW metric

def rz(z, OmegaMatter = OmegaM,
            OmegaLambda = OmegaL,
            OmegaQ = OmegaQ,
            wQ = wQ,
            eps = eps):
    kurv = OmegaMatter + OmegaLambda + OmegaQ - 1.
    nz = len(z)
    # nz=1 
    dchi = np.zeros(nz)
    chi = np.zeros(nz)
    r = np.zeros(nz)
    z2 = z[0]
    if z2 == 0.:
        dchi[0] = 0.
    else:
        dchi[0] = quad(dchi_dz, 0., z2, epsrel=eps)[0]
    for i in range(1, nz):
        z1 = z[i-1]
        z2 = z[i]
        dchi[i] = quad(dchi_dz, z1, z2, epsrel=eps)[0]
    chi[0] = dchi[0]
    for i in range(1, nz):
        chi[i] = chi[i-1] + dchi[i]
    if abs(kurv) < 1.e-4: #flat
        r = chi
    elif kurv > 0.: #closed
        r = np.sin(chi*np.sqrt(kurv))/np.sqrt(kurv)
    else: #open
        r = np.sinh(chi*np.sqrt(-kurv))/np.sqrt(-kurv)
    return r

# print(rz([0,1]))