#!/usr/bin/python
# -*- coding: utf8 -*-


#Class flowmaker
import numpy as np
import sunpy
from sunpy.image.rescale import resample
from astropy.io import fits
from convol import *
from smooth import *
from scipy.io import readsav as restore
from correlation import fivepoint

__all__ = ['flowmaker_2']
__authors__ = ["Jose Ivan Campos Rozo, Santiago Vargas Dominguez"]
__email__ = "hypnus1803@gmail.com"



def flowmaker_2(mov,lag,fwhm,reb,keyword=None):
	"""
	Compute flow maps and returns X and Y components for the proper 
	motion map.
	Inputs:
	-------
			mov: 3-D array with series of image
			lag: time-lag between 'references' and 'life' subseries.
			fwhm: fwhm for smoothing window in pixels.
			reb: rebinning	factor to change scale/range of November's
			method.
	
	Keywords:
	---------
			boxcar: if set, a boxcar window of width "fwhm" is used. Hence,
			        FWHM must be an odd number.
			adif: uses an absolute differences algorithm.
			corr: uses a multiplicative algorithm. Default is the sum of
				  square of the local differences.
			qfit2: uses 9 points fitting procedure.
			crossd: uses cross derivative interpolation formulae
	Example:
	--------
			>>> vx,vy=flowmaker(cube,1,8,1)
			
	"""
	shf=1
	std1=fwhm*0.424661
	std2=std1/reb
	#s1=size(mov)
	dims = np.ndim(mov)
	n_im = mov.shape[0]
	rows = mov.shape[1]
	columns = mov.shape[2]
	
	if dims != 3:
		raise ValueError('Array must be 3-dimensional! Breaked')
		
	s=[dims,int(columns/reb),int(rows/reb)]
	
	n=n_im-lag
	n_p=s[1]*s[2]
	
	cc=np.zeros((3,3,s[2],s[1]))
	for k in range(n):
		a=resample(mov[k,:,:],(s[2],s[1]),method='neighbor',minusone=False)
		b=resample(mov[k+lag,:,:],(s[2],s[1]),method='neighbor',minusone=False)
		a=a-a.sum()/n_p
		b=b-b.sum()/n_p
		
		for i in range(-1,2):
			for j in range(-1,2):
				if keyword=='adif':
					#~ cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:]+abs(shift(a,i*shf,j*shf)-shift(b,-i*shf,-j*shf))
					cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:]+abs(np.roll(np.roll(a,i*shf,axis=1),j*shf,axis=0)-np.roll(np.roll(b,-i*shf,axis=1),-j*shf,axis=0))
				if keyword=='corr':
					#~ cc[j+1,i+1,:,:]=cc[j+i,i+1,:,:]+shift(a,i*shf,j*shf)-shift(b,-i*shf,-j*shf)
					cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:]+np.roll(np.roll(a,i*shf,axis=1),j*shf,axis=0)-np.roll(np.roll(b,-i*shf,axis=1),-j*shf,axis=0)
				else:
					#~ dumb=shift(a,i*shf,j*shf)-shift(b,-i*shf,-j*shf)
					dumb = np.roll(np.roll(a,i*shf,axis=1),j*shf,axis=0)-np.roll(np.roll(b,-i*shf,axis=1),-j*shf,axis=0)
					cc[j+1,i+1,:,:]=cc[j+1,i+1,:,:]+dumb*dumb
					dumb = 0
		a = 0 
		b = 0
	
	cc[:,:,:,0]=cc[:,:,:,1]
	cc[:,:,0,:]=cc[:,:,1,:]
	cc[:,:,:,s[1]-1]=cc[:,:,:,s[1]-2]
	cc[:,:,s[2]-1,:]=cc[:,:,s[2]-2,:]
	
	for i in range(3):
		for j in range(3):
			if keyword=='boxcar':
				cc[j,i,:,:]=smoothe(cc[j,i,:,:],int(fwhm/reb))
			cc[j,i,:,:]=sconvol1d(cc[j,i,:,:],std=std2)
	if keyword=='qfit2':
		vx,vy=qfit2(cc)
	elif keyword=='crossD':
		vx,vy=crossD(cc)
	else:
		vx,vy=fivepoint(cc)

	vx=2.*shf*vx
	vy=2.*shf*vy
	vx=resample(vx,(rows,columns),center=True,method='neighbor',minusone=False)*reb
	vy=resample(vy,(rows,columns),center=True,method='neighbor',minusone=False)	*reb
	
	return vx,vy




		
