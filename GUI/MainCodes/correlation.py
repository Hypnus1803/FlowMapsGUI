#!/usr/bin/python
# -*- coding: utf8 -*-
import numpy as np

def fivepoint(cc):
	""" Purpose: 
		Measure the position of minimum or maximum in a 3x3 matrix.
			
	"""
	if cc==None:
		raise TypeError('Wrong number of parameters')
	dim = np.ndim(cc)
	n = [dim]
	for i in range(dim):
		i+=1
		n.append(cc.shape[dim-i])
	if (n[0] < 2) and (n[0] > 4):
		raise ValueError('Wrong input array')
	if n[n[0]-1] != 3 :#or n[n[0]] != 3:
		raise ValueError('Array must be CC[3,3,*,*]')
		
	if n[0] == 4:
		y=2.*cc[1,1,:,:]
		x=(cc[1,0,:,:]-cc[1,2,:,:])/(cc[1,2,:,:]+cc[1,0,:,:]-y)*0.5
		y=(cc[0,1,:,:]-cc[2,1,:,:])/(cc[2,1,:,:]+cc[0,1,:,:]-y)*0.5
	elif n[0] == 3:
		y=2.*cc[1,1,:]
		x=(cc[1,0,:]-cc[1,2,:])/(cc[1,2,:]+cc[1,0,:]-y)*0.5
		y=(cc[0,1,:]-cc[2,1,:])/(cc[2,1,:]+cc[0,1,:]-y)*0.5
	elif n[0] == 2:
		y=2.*cc[1,1]
		x=(cc[1,0]-cc[1,2])/(cc[1,2]+cc[1,0]-y)*0.5
		y=(cc[0,1]-cc[2,1])/(cc[2,1]+cc[0,1]-y)*0.5
	return np.round(x,6),np.round(y,6)

def qfit2(cc=None):
	""" Purpose:
		Measure the position of extrem value in a 3x3 matrix
	"""
	if cc==None:
		raise TypeError('Wrong number of parameters')
	dim = np.ndim(cc)
	n = [dim]
	for i in range(dim):
		i+=1
		n.append(cc.shape[dim-i])
	if (n[0] < 2) and (n[0] > 4):
		raise ValueError('Wrong input array')
	if n[n[0]-1] != 3 :#or n[n[0]] != 3:
		raise ValueError('Array must be CC[3,3,*,*]')
		
	#~ n=n[0]
	if n[0] == 4 :
		a1=cc[0,0,:,:]+cc[0,2,:,:]+cc[2,0,:,:]+cc[2,2,:,:]
		a2 = a1+cc[0,1,:,:]+cc[2,1,:,:]
		a1 = a1+cc[1,0,:,:]+cc[1,2,:,:]
		a3 = cc[0,0,:,:]-cc[0,2,:,:]-cc[2,0,:,:]+cc[2,2,:,:]
		a4 = -cc[0,0,:,:]+cc[2,2,:,:]
		a5 = a4-cc[0,1,:,:]-cc[0,2,:,:]+cc[2,0,:,:]+cc[2,1,:,:]
		a4 = a4+cc[0,2,:,:]-cc[1,0,:,:]+cc[1,2,:,:]-cc[2,0,:,:]
		a1 = .5*a1-cc[0,1,:,:]-cc[1,1:,:]-cc[2,1,:,:]
		a2 = .5*a2-cc[1,0,:,:]-cc[1,1,:,:]-cc[1,2,:,:]
	elif n[0] == 3:
		a1 = cc[0,0,:]+cc[0,2,:]+cc[2,0,:]+cc[2,2,:]
		a2 = a1+cc[0,1,:]+cc[2,1,:]
		a1 = a1+cc[1,0,:]+cc[1,2,:]
		a3 = cc[0,0,:]-cc[0,2,:]-cc[2,0,:]+cc[2,2,:]
		a4 = -cc[0,0,:]+cc[2,2,:]
		a5 = a4-cc[0,1,:]-cc[0,2,:]+cc[2,0,:]+cc[2,1,:]
		a4 = a4+cc[0,2,:]-cc[1,0,:]+cc[1,2,:]-cc[2,0,:]
		a1 = .5*a1-cc[0,1,:]-cc[1,1,:]-cc[2,1,:]
		a2 = .5*a2-cc[1,0,:]-cc[1,1,:]-cc[1,2,:]
	elif n[0] == 2:
		a1=cc[0,0]+cc[0,2]+cc[2,0]+cc[2,2]
		a2=a1+cc[0,1]+cc[2,1]
		a1=a1+cc[1,0]+cc[1,2]
		a3=cc[0,0]-cc[0,2]-cc[2,0]+cc[2,2]
		a4=-cc[0,0]+cc[2,2]
		a5=a4-cc[0,1]-cc[0,2]+cc[2,0]+cc[2,1]
		a4=a4+cc[0,2]-cc[1,0]+cc[1,2]-cc[2,0]
		a1=.5*a1-cc[0,1]-cc[1,1]-cc[2,1]
		a2=.5*a2-cc[1,0]-cc[1,1]-cc[1,2]
	print a1,a2,a3,a4,a5
	dim=((64./9)*a1*a2-a3**2)*1.5
	print dim
	cx=(a3*a5-((8./3)*a2*a4))/dim
	cy=(a3*a4-8./3*a1*a5)/dim
	return numpy.round(cx,6),numpy.round(cy,6)

def crossD(cc=None):
	""" Purppose:
				Measure the position of a minimum or maximum in a 3x3 array
	Example:
	--------
	>>> c=dist(3,3)
	
	>>> crossD(c)
	(0.369398, 0.369398)
	
	"""
	if cc==None:
		raise TypeError('Wrong number of parameters')
	dim = np.ndim(cc)
	n = [dim]
	for i in range(dim):
		i+=1
		n.append(cc.shape[dim-i])
	if (n[0] < 2) and (n[0] > 4):
		raise ValueError('Wrong input array')
	if n[n[0]-1] != 3 :#or n[n[0]] != 3:
		raise ValueError('Array must be CC[3,3,*,*]')
	
	if n[0] == 2:
		c4=cc[1,2]+cc[1,0]-cc[1,1]*2.
		c2=cc[1,2]-cc[1,0]
		c5=cc[2,1]+cc[0,1]-cc[1,1]*2.
		c3=cc[2,1]-cc[0,1]
		c6=(cc[2,2]-cc[2,0]-cc[0,2]+cc[0,0])/4.
	elif n[0]==3:
		c4=cc[1,2,:]+cc[1,0,:]-cc[1,1,:]*2.
		c2=cc[1,2,:]-cc[1,0,:]
		c5=cc[2,1,:]+cc[0,1,:]-cc[1,1,:]*2
		c3=cc[2,1:]-cc[0,1,:]
		c6=(cc[2,2,:]-cc[2,0,:]-cc[0,2,:]+cc[0,0,:])/4.
	elif n[0]==4:
		c4=cc[1,2,:,:]+cc[1,0,:,:]-cc[1,1,:,:]*2.
		c2=cc[1,2,:,:]-cc[1,0,:,:]
		c5=cc[2,1,:,:]+cc[0,1,:,:]-cc[1,1,:,:]*2
		c3=cc[2,1:,:]-cc[0,1,:,:]
		c6=(cc[2,2,:,:]-cc[2,0,:,:]-cc[0,2,:,:]+cc[0,0,:,:])/4.
	determ=0.5/(c4*c5 - c6*c6)
	x=determ*(c6*c3 - c5*c2)
	y=determ*(c6*c2 - c4*c3)
	return numpy.round(x,6),numpy.round(y,6)

def divergen(vx,vy):
	"""
		PURPOSE:
			Make divergence of a 2-D velocity map.
	>>> x=array([[  0.,   1.,   2.,   3.,   4.],
				 [  5.,   6.,   7.,   8.,   9.],
				 [ 10.,  11.,  12.,  13.,  14.],
				 [ 15.,  16.,  17.,  18.,  19.],
				 [ 20.,  21.,  22.,  23.,  24.]])
	
	>>> y=array([[  1. ,   1.5,   2. ,   2.5,   3. ],
				 [  3.5,   4. ,   4.5,   5. ,   5.5],
				 [  6. ,   6.5,   7. ,   7.5,   8. ],
				 [  8.5,   9. ,   9.5,  10. ,  10.5],
				 [ 11. ,  11.5,  12. ,  12.5,  13. ]])
	>>> divergen(x,y)
	array([[  5. ,   5. ,   5. ,   5. ,   5. ],
		   [  5.5,   5.5,   5.5,   5.5,  -8. ],
           [  5.5,   5.5,   5.5,   5.5, -15.5],
           [  5.5,   5.5,   5.5,   5.5, -23. ],
           [  5. ,   5. ,   5. ,   5. ,   5. ]])
	
	"""
	#~ ss=sizeidl(vx)
	dim = np.ndim(vx)
	rows = vx.shape[0]
	columns = vx.shape[1]
	if columns != len(vy[0,:]) or rows != len(vy[:,0]):
		raise ValueError('Dimensions of input arrays must be equal')
	div=np.roll(np.roll(vx,-1,axis=1),0,axis=0) - np.roll(np.roll(vx,1,axis=1),0,axis=0)#shift4d(vx,-1,0) - shift4d(vx,1,0)
	div[:,0]= -3.*vx[:,0] + 4.*vx[:,1] - vx[:,2]
	n=columns
	div[:,n-1]= -3.*vx[:,0] - 4.*vx[:,n-2] + vx[:,n-3]
	
	div=div + np.roll(np.roll(vy,0,axis=1),-1,axis=0) - np.roll(np.roll(vy,1,axis=1),1,axis=0)#shift4d(vy,0,-1) - shift4d(vy,0,1)
	div[0,:]= -3.*vy[0,:] + 4.*vy[1,:] - vy[2,:]
	n=rows
	div[n-1,:] = 3.*vy[n-1,:]-4.*vy[n-2,:] + vy[n-3,:]
	
	return div/2.
