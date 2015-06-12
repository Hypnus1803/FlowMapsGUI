#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np
import math
import time
from scipy.ndimage import correlate1d

def convol1d(array,kernel,scale_factor=None):
	"""
	The convol1d function convolves an array with a kernel 1D, 
	and returns the result. Convolution is a general process 
	that can be used for various types of smoothing, signal 
	processing, shifting, differentiation, edge detection, etc.
	"""
	
	
	row = array.shape[0]
	column = array.shape[1]
	R = np.zeros([row,column])	
	m = len(kernel)
	if scale_factor == None:
		r=correlate1d(array,kernel)
		R[:,m/2:column-math.ceil(m/2.)+1]=r[:,m/2:column-math.ceil(m/2.)+1]
	kernel=kernel/float(scale_factor)
	r=correlate1d(array,kernel)
	R[:,m/2:column-math.ceil(m/2.)+1]=r[:,m/2:column-math.ceil(m/2.)+1]
	
	return R

def sconvol1d(arreglo,kernel=None,scale_factor=None,fwhm=None,std=None):
	"""
	This program will smooth a 2D array, including the edges,
	with one-dimensional kernels. Problems of this kind arise when,
	e.g. an array is to be convolved with a 2D symmetric
	gaussian, which is separable into two one-dimensional
	convolutions.
	"""
	#~ s=len(arreglo.shape)
	dims = np.ndim(arreglo)
	rows = arreglo.shape[0]
	collumns = arreglo.shape[1]

	if dims != 2:
		raise ValueError('Array must be 2-dimensional')
	if scale_factor == None:
		scale_factor = 1.
	if kernel == None:
		if (fwhm==None) and (std==None):
			raise ValueError('Convolve with what?')
		elif fwhm != None:
			std=fwhm/(2.*math.sqrt(2.*math.log(2.)))
		#~ elif  std != None:
			#~ std=std
		elif std != None:
			width=int(std*9.)
		if width%2 == 0:
			width+=1
		kernel=np.arange(float(width))-width/2
		kernel=np.exp(-(kernel*kernel)/(2.*std*std))
		kernel=kernel/(std*math.sqrt(2.*math.pi))

	else:
		width=len(kernel)
		if width%2 == 0:
			raise ValueError('Dimension of kernel must be odd')
		
	big=np.empty([arreglo.shape[0]+width-1,arreglo.shape[1]+width-1])
	
	edge=width/2
	big[edge:big.shape[0]-edge,edge:big.shape[1]-edge]=arreglo
	for i in range(0,edge):
		big[edge:big.shape[0]-edge,i]=arreglo[:,edge-1-i]
		big[edge:big.shape[0]-edge,arreglo.shape[1]+edge+i]=arreglo[:,arreglo.shape[1]-1-i]
	big=convol1d(big,kernel,scale_factor)
	big=np.rot90(big,-1)
	for i in range(0,edge):
		big[:,i] = big[:,2*edge-1-i]
		big[:,arreglo.shape[0]+edge+i] = big[:,arreglo.shape[0]+edge-1-i]
	big=convol1d(big,kernel,scale_factor)
	big=np.rot90(big,-3)
	big=big[edge:arreglo.shape[0]+edge,edge:arreglo.shape[1]+edge]

	return big



#~ start = time.time()
#~ a=np.linspace(0,math.pi,30).reshape([5,6])
#~ kernel = np.array([1,2,3,2,1])
#~ print sconvol1d(a,std=0.10616525)
#~ print (time.time() - start), " seconds"	
