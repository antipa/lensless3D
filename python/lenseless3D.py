


import numpy as np
import time
from os import listdir
from os.path import isfile, join
from PIL import Image # Python Imaging Library

from scipy.linalg import norm

import matplotlib.pyplot as plt

from FISTA import FISTA
from AMP import AMP


# read in the point spread functions (psfs)

def readpsfs(downsampleratio = 4):
	psf_folder = '/Users/reinhardheckel/Research/lensless3D/ConvertedImages'

	psfs_paths = [psf_folder + '/' + f for f in listdir(psf_folder) if isfile(join(psf_folder, f))]


	im = Image.open( psfs_paths[0] )
	#im.show()
	(width,height) = im.size


	#psfs = [ np.array(Image.open(p))  for p in psds_paths]
	h = [] # list of psds
	for p in psfs_paths:
		im = Image.open(p)
		im = im.resize((width/downsampleratio, height/downsampleratio), Image.ANTIALIAS)
		imarr = np.array(im)[:,:,0] # just take R color channel
		imarr = imarr * (1/np.linalg.norm(imarr,'fro')) # normalize
		h.append(imarr)

	return np.array(h) # h is an tensor consisting of NZ many images

def psfsimilarity():
	'''
	examine the similarity of the psf's:
	for varying dowsampling ratios, plot how similar the psd's are
	'''
	res = []
	for downsampleratio in [2,4,8,16]:
		print 'at ', downsampleratio
		h = readpsfs(downsampleratio)
		(NZ,NX,NY) = h.shape
		corr = []
		for psd in h:
			corr.append( np.linalg.norm(psd - h[NZ/2]) / np.linalg.norm(h[NZ/2]) )
		res.append(np.array(corr))

	for ds in res:
		plt.plot(ds)
	plt.show()



class Alenseless:
	'''
	implements matrix operations for lensless imaging
	'''
	def __init__(self,h):
		self.h = h
		# define problem size
		self.NX = h[0].shape[0]
		self.NY = h[0].shape[1]
		self.NZ = len(h)
	
	def pad(s,A):
		return np.lib.pad(A,((s.NX/2,np.ceil(s.NX/2.0)),(s.NY/2,np.ceil(s.NY/2.0))),'constant') 
	
	def crop(s,B):
		return B[s.NX/2:3*s.NX/2,s.NY/2:3*s.NY/2] # crop
	
	def applyA(self,x):
		# y = Ax = sum_z h^z \ast x
		# A consists of all shifts of h in 
		# x: variable to convolve with h, (NY, NX, NY) array
		Y = np.zeros((2*self.NX,2*self.NY),dtype=complex) 
		for psd,xplane in zip(self.h,x):
			Y += np.fft.fft2( self.pad(xplane) ) * np.fft.fft2( self.pad(psd) )
		
		return np.array(np.real( self.crop( np.fft.ifftshift( np.fft.ifft2(Y) ) ) ) )
	
	def applyAH(self,y):
		x =	[] # signal to return
		Y = np.fft.fft2(self.pad(y))
		for psf in self.h:
			H = np.conj( np.fft.fft2( self.pad(psf) ) )
			xplane = np.real( self.crop( np.fft.ifftshift( np.fft.ifft2( H*Y ) ) ) )
			x.append(xplane)
		return np.array(x)

	def estimate_Lipschitz(self, nit = 1000):
		'''
		Estimate of 1 / || A.T A || from nit power iterations
		'''
		x = np.zeros((self.NZ, self.NX, self.NY),dtype=complex) 
		x = np.random.randn(self.NZ, self.NX, self.NY)
		
		for i in range(nit):
			x =	self.applyAH( self.applyA(x) )
			x = x / np.linalg.norm(x)
		x = self.applyAH( self.applyA(x) )
		return 1.0 / np.linalg.norm(x)

	def grad(self,x,y):
		residual = self.applyA(x) - y
		return [ self.applyAH(residual) , np.linalg.norm(residual,'fro')**2 ]


# proxy's for FISTA

def proxy_nonneg(x,alpha):
	return np.maximum( x - alpha , 0 )

def soft_thresholding(x,alpha):
	return np.sign(x) * np.maximum( np.abs(x) - alpha , 0. )

### tests 

def test_simple_pointsource():
	'''
	identify the position of two point sources in a simple example with the original psf's
	'''
	downsampleratio = 16
	print 'read images downsampled by', downsampleratio, '...'
	h = readpsfs(downsampleratio)
	h = np.array( [ h[1,:,:], h[25,:,:], h[50,:,:] ] ) # a toy example
	(NZ,NX,NY) = h.shape
	A = Alenseless(h)
	# generate a test instance
	xast = np.zeros(A.h.shape)
	xast[0, NX/2, NY/2] = 1
	xast[1,NX/2-2,NY/2-2] = 1
	y = A.applyA(xast) 

	def print_cmp(xhat,xast):
		w = 3
		print "original point source: "
		print xast[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]
		print xast[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]
		print "estimated point sources: "
		print xhat[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]
		print xhat[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]

	### estimation with AMP

	amp = AMP()
	amp.maxsteps = 100
	amp.tau = 100

	print 'fit AMP...'
	x0  = np.random.rand(NZ,NX,NY) # random initializer
	xhat = amp.fit( x0 , y, A.applyA, A.applyAH, proxy_nonneg)
	print_cmp(xhat,xast)


	### estimation with FISTA

	print 'estimate Lipschitz constant...'
	stepsize = 0.8*A.estimate_Lipschitz() # factor of 0.8 to account for uncertainty in estimating the Lipschitz constant

	fista = FISTA()
	fista.maxsteps = 500
	fista.lam = 0.1
	fista.stepsize = stepsize

	print 'fit...'
	start_time = time.time()
	x0  = np.random.rand(NZ,NX,NY) # random initializer
	xhat = fista.fit( x0 , y, A.grad, proxy_nonneg )
	#xhat = fista.fit( x0 , grad, soft_thresholding)
	end_time = time.time()

	print "time: ", (end_time - start_time)
	print "rel. reconstruction error: ", np.linalg.norm( xast - xhat )/np.linalg.norm(xast)
	print "rel. residual: ", np.linalg.norm( A.applyA(xast - xhat) )/np.linalg.norm(y)

	w = 3
	print "original point source: "
	print xast[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]
	print xast[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]
	print "estimated point sources: "
	print xhat[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]
	print xhat[1,NX/2-2-w+1:NX/2-2+w,NY/2-2-w+1:NY/2-2+w]


def test_random_psf():
	'''
	a simple test case for lenseless imaging with FISTA reconstruction for Gaussian psf's
	'''

	NZ = 2
	NX = 5
	NY = 4
	# generate Gaussian psfs
	h = []
	for z in range(NZ):
		hz = np.random.randn(NX,NY)
		hz = hz / np.linalg.norm(hz)
		h.append(hz)
	h = np.array( h )

	A = Alenseless(h)

	# generate a test instance
	xast = np.zeros(A.h.shape)
	xast[0, NX/2, NY/2] = 1
	xast[1,NX/2-2,NY/2-2] = 1

	y = A.applyA(xast) 

	stepsize = 0.8*A.estimate_Lipschitz() # factor of 0.8 to account for uncertainty in estimating the Lipschitz constant

	fista = FISTA()
	fista.maxsteps = 500
	fista.lam = 0.001
	fista.stepsize = stepsize

	start_time = time.time()
	x0  = np.random.rand(NZ,NX,NY) # random initializer
	xhat = fista.fit( x0 , y, A.grad, proxy_nonneg )
	#xhat = fista.fit( x0 , grad, soft_thresholding)
	end_time = time.time()

	print "time: ", (end_time - start_time)
	print "reconstruction error: ", np.linalg.norm( xast - xhat )/np.linalg.norm(xast)


def test_estimate_Lipschitz(A,nit = 1000):
	'''
	estimate largest singular value of the matrix A.T*A with power iterations
	'''
	(m,n) = A.shape
	x = np.random.randn(n)
	
	for i in range(nit):
		x =	np.dot(A.T, np.dot(A,x) )
		x = x / np.linalg.norm(x)

	x = np.dot(A.T, np.dot(A,x) )
	return np.linalg.norm(x), np.linalg.norm( np.dot(A,A.T), 2)


