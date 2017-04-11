


import numpy as np
import time



class AMP:
	'''
	Problem: estimate sparse x from
		y = Ax + z, A is n x N
		A normalized such that the columns have approximatelly unit norm
	AMP algorithm: (e.g. Donoho et al. ``The noise-sensitivity phase transition in CS'', 2011, p. 8)
		z^t     = y - A x^t + b_t z^{t-1},  b_t = ||x^t||_0 / n
		x^{t+1} = softthresholding( A^H z^t + x^t , theta_t )
		theta_t = tau * sigma_t
		estimate mean square error as: sigma_t = sqrt( ||z_t||_2^2 / n )
		tau is a tuning parameter
	'''
	def __init__(self):
		# set default parameters
		self.maxsteps = 100
		self.tau = 2
	
	def fit(self,x0,y, applyA, applyAH, proxy):
		'''
		x0:			initial point
		y:			measurement
		applyA: 	applyA(x) = A*x
		applyAH: 	applyAH(y) = A.T*y
		proxy: 		soft thresholding or non-negative soft thresholding
		'''
		n = y.size
		x = x0
		z = y
		for i in range(self.maxsteps):
			bt = np.sum(x>0.0)/n
			z = y - applyA(x) + bt*z
			thetat = self.tau * np.linalg.norm(z)/np.sqrt(n)
			x = proxy(  applyAH(z) + x , thetat )

		return x


def test_AMP():
	'''
	simple test
	'''
	
	def soft_thresholding(x,alpha):
		return np.sign(x) * np.maximum( np.abs(x) - alpha , 0. )

	def nonneg_proxy(x,alpha):
		return np.maximum( x - alpha , 0. )

	### generate a problem instance
	n = 400
	m = 100
	s = 10

	A = np.random.randn(m,n)
	A = 1/np.linalg.norm(A,axis=0)*A # normalize

	applyA = lambda x: np.dot(A,x) 
	applyAH = lambda y: np.dot(A.T,y) 

	x = np.zeros(n)
	x[:s] = np.abs(np.random.randn(s))
	x = x/np.linalg.norm(x)
	y = np.dot(A,x)
	lam = 0.01

	amp = AMP()

	start_time = time.time()
	xhat = amp.fit( np.zeros(n) , y, applyA, applyAH, soft_thresholding)
	end_time = time.time()
	print "AMP soft thresholding"
	print "reconstruction error: ", np.linalg.norm(x-xhat)/np.linalg.norm(x)
	print "time elapsed: ", end_time - start_time

	start_time = time.time()
	xhat = amp.fit( np.zeros(n) , y, applyA, applyAH, nonneg_proxy)
	end_time = time.time()
	print "AMP non-negative proxy"
	print "reconstruction error: ", np.linalg.norm(x-xhat)/np.linalg.norm(x)
	print "time elapsed: ", end_time - start_time

