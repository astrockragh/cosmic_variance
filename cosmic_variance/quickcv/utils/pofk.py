import numpy as np

def pofk(karr, gamma=None, n=None):

    """ For a single k, return the power spectrum"""
    h=0.678
    omega0=0.308
    omegaB=0.022/h**2
    sigma8 = 0.82

    if not n:
        n = 0.96 #after Planck
    if not gamma:
        gamma = omega0*h*np.exp(-omegaB*(1.+np.sqrt(2*h)/omega0))

    #  determine sigma8 from model
    #  following Peacock, p. 528

    keff=(0.172+0.011*(np.log(gamma/0.34))**2)*1.
    q=keff/(h*gamma)
    tk=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
    sigma8squn=keff**(3.0+n)*tk**2

    q=karr/(h*gamma)   # k/h Gamma
    tk=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
    delsq=karr**n*tk**2       #took out factor of k**3
    delsq=delsq*(sigma8)**2/sigma8squn        #normalize to Borgani et al. sigma 8, omegam=0.31,omegal=0.69

    pk=2*np.pi**2*delsq

    return pk


def pofkint(k0, k1, k2, gamma=None, n=None):

    """ For 3d k, calculate a single k, return the power spectrum"""
    k = np.linalg.norm(np.array([k0, k1, k2]))
    h=0.678
    omega0=0.308
    omegaB=0.022/h**2
    sigma8 = 0.82

    if not n:
        n = 0.96 #after Planck
    if not gamma:
        gamma = omega0*h*np.exp(-omegaB*(1.+np.sqrt(2*h)/omega0))

    #  determine sigma8 from model
    #  following Peacock, p. 528

    keff=(0.172+0.011*(np.log(gamma/0.34))**2)*1.
    q=keff/(h*gamma)
    tk=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
    sigma8squn=keff**(3.0+n)*tk**2

    q=k/(h*gamma)   # k/h Gamma
    tk=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
    delsq=k**n*tk**2       #took out factor of k**3
    delsq=delsq*(sigma8)**2/sigma8squn        #normalize to Borgani et al. sigma 8, omegam=0.31,omegal=0.69

    pk=2*np.pi**2*delsq

    return pk

# print(pofk(np.array([0.1,1,10])))