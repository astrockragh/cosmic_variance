import numpy as np

def pofkint(k0, k1, k2, OmegaM = 0.308, OmegaBaryon = 0.022/(0.678)**2, sigma8 = 0.82, ns = 0.96, h = 0.678):

    """ For 3d k, calculate a single k and returns the power spectrum

    OmegaM is the matter density parameter
    OmegaBaryon is the baryon density parameter, here it is set to the Planck value, 0.022/h^2
    sigma8 is the normalization of the power spectrum
    ns is the spectral index
    h is the hubble parameter H0/100

    returns the power at k = sqrt(k0**2 + k1**2 + k2**2)
    """

    k = np.linalg.norm(np.array([k0, k1, k2]))

    gamma = OmegaM*h*np.exp(-OmegaBaryon*(1.+np.sqrt(2*h)/OmegaM))

    keff=(0.172+0.011*(np.log(gamma/0.34))**2)*1.
    q=keff/(h*gamma)
    tk=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
    sigma8squn=keff**(3.0+ns)*tk**2

    q=k/(h*gamma)   # k/h Gamma
    tk=np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**(-0.25)
    delsq=k**ns*tk**2     
    delsq=delsq*(sigma8)**2/sigma8squn #normalize to Borgani et al. 

    pk=2*np.pi**2*delsq

    return pk