#!/usr/bin/python
# -*- coding: utf8 -*-

import numpy as np
import math
from scipy.ndimage import uniform_filter

def smooth(array,width):
	"""
	The SMOOTH function returns a copy of Array smoothed with a boxcar 
	average of the specified width. The result has the same type and 
	dimensions as Array. 
	"""
	if type(width) != np.ndarray:
		if (type(width) != list) and (type(width) != tuple): 
			width = np.array([width])
		else:
			width = np.array(width)
	row = array.shape[0]
	column = array.shape[1]
	dims = np.ndim(array)
	if len(width) == 1:
		r = uniform_filter(array,size=width[0],mode="nearest")
		data_copy=array.copy()
		data_copy[math.ceil((width[0]-1)/2.):row-(width[0]+1)/2+1,math.ceil((width[0]-1)/2.):column-(width[0]+1)/2+1]=r[math.ceil((width[0]-1)/2.):row-(width[0]+1)/2+1,math.ceil((width[0]-1)/2.):column-(width[0]+1)/2+1]
	if len(width) == 2:
		r = uniform_filter(array,size=width,mode="nearest")
		data_copy=array.copy()
		data_copy[math.ceil((width[0]-1)/2.):row-(width[0]+1)/2+1,:]=r[math.ceil((width[0]-1)/2.):row-(width[0]+1)/2+1,:]
		data_copy[:,math.ceil((width[1]-1)/2.):column-(width[1]+1)/2+1]=r[:,math.ceil((width[1]-1)/2.):column-(width[1]+1)/2+1]
	return data_copy

def smoothe(arreglo,smoonum,keyword=None):
	"""
	This program will smooth an array, including the edges,
	using the smooth function (see smooth) by surrounding the array
	with duplicates of itself and then smoothing the large
	array.
	"""
	
	row = arreglo.shape[0]
	column = arreglo.shape[1]
	dim = np.ndim(arreglo)
	border = smoonum*2
	def oned():
		eg1=column + (smoonum-1)
		bigarr=np.zeros(border+column)
		bigarr[smoonum:eg1+1] = arreglo
		rott=arreglo[::-1]
		bigarr[0:smoonum]=rott[column-smoonum:column]
		bigarr[column+smoonum:column+(2*smonum)]=rott[0:smoonum]
		sbigar=smooth(bigarr,smoonum)
		done1=sbigarr[smoonum:eg1+1]
		return done1

		
	def twod():
		eg1=column+(smoonum-1)
		eg2=row+(smoonum-1)
		max1=column+border-1
		max2=row+border-1
		
		bigarr=np.zeros([(border+row),(border+column)])
		bigarr[smoonum:eg2+1,smoonum:eg1+1]=arreglo
		rott=np.rot90(arreglo[::-1],2)
		bigarr[smoonum:eg2+1,0:smoonum]=rott[:,column-smoonum:column]
		bigarr[smoonum:eg2+1,column+smoonum:max1+1]=rott[:,0:smoonum]
		rott=np.rot90(arreglo[::-1],4)
		bigarr[0:smoonum,smoonum:eg1+1]=rott[row-smoonum:row,:]
		bigarr[row+smoonum:max2+1,smoonum:eg1+1]=rott[0:smoonum,:]
		rott=np.rot90(arreglo,2)
		bigarr[0:smoonum,0:smoonum]=rott[row-smoonum:row,column-smoonum:column]
		bigarr[row+smoonum:max2+1,column+smoonum:max1+1]=rott[0:smoonum,0:smoonum]
		bigarr[row+smoonum:max2+1,0:smoonum]=rott[0:smoonum,column-smoonum:column]
		bigarr[0:smoonum,column+smoonum:max1+1]=rott[row-smoonum:row,0:smoonum]
		sbigarr=smooth(bigarr,smoonum)
		done2=sbigarr[smoonum:eg2+1,smoonum:eg1+1]
		return done2
		
	if dim == 0:
		done=-1
		print "That ain't an array! Its a scalar."

	if dim == 1:
		done=oned()
		
	if dim == 2:
		done=twod()
		
	if dim == 3:
		done=-1
		print 'Sorry, I am too tired to process such a large array.'
	
	return done	












#~ start = time.time()
#~ data3=np.sin(dist(8,7))/3
#~ width = 5
#~ print smoothe(data3,width)
#~ print data3 - smoothe(data3,width)
#~ print (time.time() - start), " seconds"	

		
	
	
		
	
