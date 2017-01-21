

import matplotlib.pyplot as plt
import numpy as np
import time

import sys as syst
from os.path import dirname
syst.path.append('./SPGL1_python_port')
from spgl1 import *
from spgl_aux import *


from fista_JeanKossaifi import Fista


#def test_fista():
'''
test fista for the lasso problem, with x non-negative
F(x) = ||Ax - y||^2 + lam*||x||_1
'''


class FISTA:
	'''
	FISTA algorithm for solving:
		minimize_x f(x) + lam*g(x)
	'''
	def __init__(self):
		# set default parameters
		self.stepsize = 5e-5 #0.00005
		self.tolerance = 1e-6
		self.maxsteps = 3000
		self.lam = 0.01
	
	def fit(self,x0,y,grad,proxy):
		'''
		x0:			initial point
		grad: 		function returning [ gradient f(x) , f(x) ], takes x, y as arguments
		proxy: 		proximity operator
		'''
		y_current = np.array(x0, dtype=np.float) 
		y_next = np.array(x0, dtype=np.float)
		xk = np.copy(y_next) # linear combination of vectors of the two last iterations
		tau_1 = 1

		for i in range(self.maxsteps):
			y_current = y_next # B_(k-1) = B_(k)
			[ gradxk , fxk ] = grad(xk,y)
			y_next = proxy( xk - self.stepsize*gradxk , self.stepsize*self.lam )

			tau_0 = tau_1 #tau_(k+1) = tau_k
			tau_1 = (1 + np.sqrt(1 + 4*tau_0**2))/2

			xk = y_next + (tau_0 - 1)/tau_1*(y_next - y_current)
	
			if abs(fxk) < self.tolerance:
				break

		return xk



def test_FISTA():
	'''
	compare the fista implementation to another FISTA implementation and to the SPGL1 solver
	'''

	### generate a problem instance
	n = 400
	m = 100
	s = 10

	A = np.random.rand(m,n)
	x = np.zeros(n)
	x[:s] = np.abs(np.random.randn(s))
	x = x/np.linalg.norm(x)
	y = np.dot(A,x)
	lam = 0.01

	stepsize = 1/np.linalg.norm(np.dot(A, A.T), 2) # Lipschitz constant; largest singular value

	### solve with FISTA version from the web

	lam = 0.01
	fis = Fista(lam,'least-square')
	fis = Fista(lam,'least-square','l11',3000)

	start_time = time.time()
	fis.fit(A,y,stepsize)
	end_time = time.time()

	xhat = fis.coefs_
	print "FISTA from the web"
	print "reconstruction error: ", np.linalg.norm(x-xhat)/np.linalg.norm(x)
	print "time elapsed: ", end_time - start_time

	### solve with my FISTA version

	def grad(xt,y):
		'''
		f(x) := ||Ax - y||^2, grad f(x) = 2 A^T(Ax-y)
		function returns:
			[ gradient f(x) ,  f(x) ]
		'''
		residual = np.dot(A,xt) - y
		return [ np.dot( A.T , residual ) , np.linalg.norm(residual)**2 ]

	def soft_thresholding(x,alpha):
		return np.sign(x) * np.maximum( np.abs(x) - alpha , 0. )

	#def proxy_nonneg(x,alpha=lam*stepsize):
	#	return np.maximum( x - alpha , 0 )


	fista = FISTA()
	fista.stepsize = stepsize
	fista.lam = 0.01 # lasso parameter

	start_time = time.time()
	xhat = fista.fit( np.zeros(n) , y, grad, soft_thresholding )
	end_time = time.time()

	print "FISTA"
	print "reconstruction error: ", np.linalg.norm(x-xhat)/np.linalg.norm(x)
	print "time elapsed: ", end_time - start_time

	### solve with SPGL1
	
	opts = spgSetParms({
				'iterations' :   10000  # Set the maximum number of iterations to a large number (default is 10* .. )
				})
	sigma = 0.02
	start_time = time.time()
	xhat,resid,grad,info = spg_bpdn(A, y, sigma, opts) # Call the SPGL1 solver 
	end_time = time.time()

	print "SPGL1"
	print "reconstruction error: ", np.linalg.norm(x-xhat)/np.linalg.norm(x)
	print "time elapsed: ", end_time - start_time



