#!/usr/bin/python
# -*- coding: utf8 -*-



import numpy as np
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
	Script to generate the flow maps for the vector fields.
	
	"""
	print('The next variables are necessaries: lag, fwhm_arcsec, reb, pix, t_step, fwhm, kmperasec, h_m, char')
	lag=1#input('Introduce the lag between the images to be compared (number of images):')
	fwhm_arcsec=40 #0.5#input('Introduce the fwhm of the window for tracking(arcsec) (Could be 3., 0.8,1.2,0.5,1.):')
	print("Este es fwhm:",fwhm_arcsec)
	reb=1.#input('Introduce the rebinning factor:')
	pix=0.504*4#0.504#0.0544#input('Size of pixel (pix):')
	t_step=45.#input('temporal sampling interval in the time series (seconds):')
	fwhm=fwhm_arcsec/pix
	print('fwhm is calculated automatically:',fwhm)
	kmperasec=725#input('Value of kilometers per arcsec:')
	h_m=150#input('mass-flux scale-heigth (November 1989, ApJ,344,494):')
	#~ char=2.#input('Size of the Char for the histograms:')
	#Calculos derivados de los parametros anteriores
	v_limit=2*reb+reb #cota maxima velocidad en pixeles.
	delta_t=t_step*lag # time-lag in seconds
	factor=pix*kmperasec/delta_t
	print('-------------------------------------------------------')
	print('Maximum velocity that can be detected-->',v_limit*factor,'km/s')
	print('FWHM of the window for tracking (pix/arc)-->',fwhm,'/',fwhm_arcsec)
	print('-------------------------------------------------------')
	#************************************************************
	
	num=cube.shape[0]
	aver=cube[0,:,:]*0
	print('Number of images', num)
	print('Generando imagen promedio del cubo ...')

	print('Calculando flujos horizontales ...')

	vx,vy=flowmaker_2(cube,lag,fwhm,reb)
	vx=vx.clip(min=-v_limit,max=v_limit) ; vy=vy.clip(min=-v_limit,max=v_limit)

	print('Computes the divergence and vertical velocities ...')
	
	vx_kps=vx*factor 	#vx in km/s
	vy_kps=vy*factor	#vy in km/s
	div= divergen(vx_kps,vy_kps)
	vz_kps=h_m*div
	
	vx_kps= vx_kps-vx_kps.mean()
	vy_kps= vy_kps-vy_kps.mean()
	vz_kps= vz_kps-vz_kps.mean()
	
	mag = np.sqrt(vx_kps*vx_kps + vy_kps*vy_kps)
	
	return vx_kps, vy_kps, vz_kps, mag
	
	


