# ;+
# ; NAME:
# ;
# ;  QUICKCV
# ;
# ; PURPOSE:
# ;
# ;   calculate cosmic variance for some geometry
# ;
# ; CATEGORY:
# ;
# ;   cosmic variance
# ;
# ;
# ; CALLING SEQUENCE:
# ;
# ;   fracerror=quickcv(side1,side2 [zarr,deltaz,sig8=sig8,npoints=npoints]) 
# ;
# ; INPUTS:
# ;
# ;   side1: length of 1 side of field on sky, in degrees (rect. geom)
# ;   side2: length of other side of field on sky, in degrees
# ;
# ; OPTIONAL INPUTS:
# ;   zarr: array of z's for which cosmic variance will be calculated; default:1
# ;   deltaz: total length of volume in z direction; default: 0.1
# ;
# ; KEYWORD PARAMETERS:
# ;   sig8: desired value of sigma_8; default: 0.77
# ;   npoints: desired number of integration points for first part of 
# ;     integration; default: 200
# ;
# ;
# ; OUTPUTS:
# ;
# ;  fracerror: array containing the fractional error in a count,
# ;  sigma/N, due to cosmic variance for objects with bias b=1 
# ;  in a volume side1 x side2 degrees x delta_z centered at redshifts
# ;  zarr.  NOTE THIS IS SIGMA, NOT VARIANCE!!!!
# ;
# ; RESTRICTIONS:
# ;
# ;  uses power spectrum in pofk.pro (gamma=0.1872, omega0=0.26), distance
# ;  relation in rz.pro 
# ;
# ;
# ; MODIFICATION HISTORY:
# ; modified by rss july 2006 to include Dlin scaling with z
# ; modified by bpm sep  2009 to wmap3 cosmology
# ;-

import numpy as np
from .utils.rz import rz
from .utils.intpk4 import intpk4
from .utils.dlin import Dlin

def quickcv(side1,side2,za,deltaz, omegam = 0.31, omegal = 0.69, sig8=0.82, acc='low'):
  """za and deltaz have to be np.arrays"""
  #sides in degrees, z is redshift bin middle, deltaz is redshift bin width, everything else is cosmology and integration parameters
  #see intpk4.py for acc options

  # assume LCDM

  rza=rz(za, OmegaMatter=omegam, OmegaLambda=omegal)
  rzmin=rz(za-deltaz/2, OmegaMatter=omegam, OmegaLambda=omegal)
  rzmax=rz(za+deltaz/2, OmegaMatter=omegam, OmegaLambda=omegal)

  radeg = 57.2958 # degrees per radians

  x1a=3000./2.*(rza * side1)/radeg
  x2a=3000./2.*(rza * side2)/radeg
  x3a=3000.*(rzmax-rzmin)/2.
  xs = np.array([x1a, x2a, x3a])
  cvi = intpk4(xs, acc =  acc)

  #cvi is the fractional variance in a count in a rectangle
  #cvar is the fractional error (sigma), NOT VARIANCE

  cvi = cvi*sig8**2 #normalize by sig8^2
  d = Dlin(omegam, za)
  cvar = np.sqrt(cvi)*d

  return cvar[0]

# print(quickcv(1,2, np.array([3]), np.array([1])))