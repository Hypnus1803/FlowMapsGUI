#!/usr/bin/python
# -*- coding: utf8 -*-


import sys
import numpy as np
import time
import numpy.ma as ma
from flowmaker_2 import *
from scipy.io import readsav as restore
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from correlation import divergen
from astropy.io import fits


def flow(cube):
	"""
	Programa para generar las imagenes de los flujos vectoriales
	Las entradas se van poniendo via  consola.
	
	"""
	print 'The next variables are necessaries: lag, fwhm_arcsec, reb, pix, t_step, fwhm, kmperasec, h_m, char'
	lag=1#input('Introduce the lag between the images to be compared (number of images):')
	fwhm_arcsec=10.#0.5#input('Introduce the fwhm of the window for tracking(arcsec) (Could be 3., 0.8,1.2,0.5,1.):')
	reb=1.#input('Introduce the rebinning factor:')
	pix=1.9856#0.0544#input('Size of pixel (pix):')
	t_step=96.*60.#input('temporal sampling interval in the time series (seconds):')
	fwhm=fwhm_arcsec/pix
	print 'fwhm is calculated automatically:',fwhm
	kmperasec=725#input('Value of kilometers per arcsec:')
	h_m=150#input('mass-flux scale-heigth (November 1989, ApJ,344,494):')
	char=2.#input('Size of the Char for the histograms:')
	#Calculos derivados de los parametros anteriores
	v_limit=2*reb+reb #cota maxima velocidad en pixeles.
	delta_t=t_step*lag # time-lag in seconds
	factor=pix*kmperasec/delta_t
	print '-------------------------------------------------------'
	print 'Maximum velocity that can be detected-->',v_limit*factor,'km/s'
	print 'FWHM of the window for tracking (pix/arc)-->',fwhm,'/',fwhm_arcsec
	print '-------------------------------------------------------'
	#************************************************************
	
	num=cube.shape[0]
	aver=cube[0,:,:]*0
	print 'Number of images', num
	print 'Generando imagen promedio del cubo ...'
	for i in range(num-1):
		aver=aver+cube[i,:,:]
	aver=aver/num
	print 'Calculando flujos horizontales ...'
	#~ vx,vy=flowmaker(cube[:-1,:,:],lag,fwhm,reb)
	vx,vy=flowmaker_2(cube,lag,fwhm,reb)
	vx=vx.clip(min=-v_limit,max=v_limit) ; vy=vy.clip(min=-v_limit,max=v_limit)

	print 'Computes the divergence and vertical velocities ...'
	
	vx_kps=vx*factor 	#vx in km/s
	vy_kps=vy*factor	#vy in km/s
	div= divergen(vx_kps,vy_kps)
	vz_kps=h_m*div
	mag = np.sqrt(vx_kps*vx_kps + vy_kps*vy_kps)
	
	print 'Plotting velocities map ...'
	fig = plt.figure()
	ax = fig.add_subplot(111) 
	X,Y=np.meshgrid(np.arange(0,aver.shape[1],1),np.arange(0,aver.shape[0],1))
	#~ plt.xlim(0,(aver.shape[0]-1)*pix)
	#~ plt.ylim(0,(aver.shape[1]-1)*pix)	
	ax.imshow(aver.clip(min=-700,max=700),origin='low',cmap='gray')#, extent=(0,(aver.shape[1]-1)*pix,0,(aver.shape[0]-1)*pix))
	#~ ax.set_title('Velocity field applying LCT')
	
	M2=ma.masked_greater(aver,-150)
	M1=ma.masked_less(aver,150)
	vx1_kps = ma.masked_array(vx_kps,mask=M1.mask)
	vy1_kps = ma.masked_array(vy_kps,mask=M1.mask)
	vx2_kps = ma.masked_array(vx_kps,mask=M2.mask)
	vy2_kps = ma.masked_array(vy_kps,mask=M2.mask)
	
	plt.quiver(X[::15,::15],Y[::15,::15],vx1_kps[::15,::15],vy1_kps[::15,::15],pivot='tail',color='blue',units='dots',scale=0.003,width=0.8)
	plt.quiver(X[::15,::15],Y[::15,::15],vx2_kps[::15,::15],vy2_kps[::15,::15],pivot='tail',color='red',units='dots',scale=0.003,width=0.8)
	#~ plt.axis([0,(aver.shape[1]-1)*pix,0,(aver.shape[0]-1)*pix])
	#~ plt.quiver(X,Y,vx_kps,vy_kps,pivot='tail',color='blue',units='dots',scale=0.01,width=1.2)		
	#~ plt.quiver(X[::3,::3],Y[::3,::3],vx_kps[::3,::3],vy_kps[::3,::3],pivot='tail',color='red',units='dots',scale=0.01,width=1.2)
	ticksx = ticker.FuncFormatter(lambda vx1_kps, pos: '%.0f' % int(vx1_kps*pix))                                                                                                                                         
	ax.xaxis.set_major_formatter(ticksx)
	ticksy = ticker.FuncFormatter(lambda vy1_kps, pos: '{0:f}'.format(vy1_kps*pix))                                                                                                                                         
	ax.yaxis.set_major_formatter(ticksy)
	plt.draw()
	print (time.clock() - start)
	plt.show()
	

start = time.clock()
filename='/home/hypnus1803/Desktop/TESIS_MSc./Codes/PYTHON/cube_aligne2.fit'
#~ cube=(restore(filename)).cube
#~ cube = cube[:,:,15:265]
#~ cube = cube[61:67,67:182,32:197]		
hdu= fits.open(filename)
cube = hdu[0].data
flow(cube)	

