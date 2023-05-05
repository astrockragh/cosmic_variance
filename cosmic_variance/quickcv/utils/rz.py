

# Routines for computing the Alcock-Paczynski distortion

# January 2023, Christian Kragh Jespersen, rewritten from IDL


# H(z) in units where H(0) = 1

import numpy as np 
from scipy.integrate import quad

#get hubble parameter at redshift z relative to h at z=0
def hubble(z, OmegaM = 0.308, OmegaL = 0.692):
    x = 1. + z
    OmegaK = 1. - OmegaM - OmegaL #curvature.
    h2 = OmegaM * x**3 + OmegaK * x**2 + OmegaL
    return np.sqrt(h2)

# Compute d(chi)/dz where d(chi)^2 = dr^2 / (1-kr^2)
def dchi_dz(z):
    return 1./hubble(z)

# Compute coordinate distance r(z) for FRW metric
def rz(z, OmegaM = 0.308, OmegaL = 0.692, eps = 1e-10):
    '''Compute coordinate distance r(z) for FRW metric, given cosmological parameters OmegaM and OmegaL and redshift z
    eps parameter is the precision for scipy.integrate.quad'''
    kurv = OmegaM + OmegaL - 1.
    nz = len(z)
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