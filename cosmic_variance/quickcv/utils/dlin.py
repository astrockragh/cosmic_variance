def Dlin(OmegaM, z): 
  '''Approximation for linear growth factor, assumes a flat cosmology, from Sean Carroll's paper
   OmegaM is the matter density parameter
   z is the redshift'''
  Omega = OmegaM*(1+z)**3/(1-OmegaM + (1+z)**3*OmegaM)
  gz = 2.5*Omega/(1./70. + 209.0*Omega/140. + Omega**(4./7.))
  g0 = 2.5*OmegaM/(1./70. + 209.0*OmegaM/140. + OmegaM**(4./7.))
  D = gz/(g0*(1+z))

  return D