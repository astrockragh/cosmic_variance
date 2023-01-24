from .pofk import pofkint
import numpy as np
from scipy.integrate import nquad

import time

# xss = np.array([np.array([41.126813, 82.253627, 280.86420]), np.array([45.320475, 90.640949, 208.52065]), np.array([49.203843, 98.407686, 231.35671])]) #the unit is comoving megaparsec

def intpk4(xs, acc = 'high', verbose = False):

	"""" Solver for integral in equation 2 of https://arxiv.org/pdf/astro-ph/0109130.pdf
	xs are the side lengths of the box in comoving Mpc
	acc is the accuracy of the integral, 'high' or 'low', low is about 40 times faster and has about a 2% error
	verbose is whether to print out the time taken for the integral

	This version is written by Christian Kragh Jespersen, January 2023, rewritten from the IDL code by Jonathan Baker
	"""
	if acc == 'high':
		acc = 1
	elif acc == 'low':
		acc = 5

	# kmax = 4*np.pi/xs ## this is the maximum k value for the window function that contributes significantly to the integral

	kmax = (0.4+16/xs)*1.25 #this is to match from the idl code
	kmax = kmax/acc

	def window_idl(k0, k1, k2):
		return np.sin(k0*xs[0])/(k0*xs[0])*np.sin(k1*xs[1])/(k1*xs[1])*np.sin(k2*xs[2])/(k2*xs[2]) 

	def integrand_idl(k0,k1,k2):
		return window_idl(k0, k1, k2)**2*pofkint(k0,k1,k2)

	start = time.time()

	integral = nquad(integrand_idl, [[0, kmax[0]], [0, kmax[1]], [0, kmax[2]]], opts = {"limlst": 10, "maxp1": 10})
	end = time.time()

	if verbose:
		print(f'Integration time{end - start:.2f} seconds')

	#factor of pi**3 is to normalize volume in k space, there are offsetting factors of 8 not included here
	return integral[0]/np.pi**3
