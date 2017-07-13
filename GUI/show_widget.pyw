# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:11:48 2015

@author: hypnus1803

"""
import _pickle as pk
#import cPickle as pk
from numpy import ma
import sys,os
sys.path.append(os.getcwd()+'/'+'MainCodes')


from scipy.io.idl import readsav
from astropy.io import fits


import matplotlib as mpl
mpl.use('TkAgg')
mpl.rcParams.update({'font.size': 10})
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.ticker as ticker  
from matplotlib.widgets import  RectangleSelector
from matplotlib.patches import Rectangle
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


import numpy as np
from flowmaker_2 import *
from correlation import divergen

from skimage import measure
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from Vector_Velocities_with_Sunpy import *
from plotOptions import *
from fitOptions import *

import sunpy
__authors__ = [u"Jose Iván Campos Rozo", u"Santiago Vargas Domínguez"]
__email__ = "jicamposr@unal.edu.co"


try:
	_fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
	def _fromUtf8(s):
		return s


# Matplotlib Figure object 1 (Main Image)
from matplotlib.figure import Figure

#~ from double_slider import QHRangeSlider


class ImageCanvas1(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# setup Matplotlib Figure and Axis
		self.fig = Figure(figsize=(10,10))
		self.ax = self.fig.add_subplot(111)
		
		self.fig.set_tight_layout(True)

		# initialization of the canvas
		FigureCanvas.__init__(self, self.fig)
		# we define the widget as expandable
		FigureCanvas.setSizePolicy(self,
								   QtGui.QSizePolicy.Preferred,
								   QtGui.QSizePolicy.Preferred)
		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)
		self.fig.set_facecolor('white')

# Matplotlib Figure object 2 (Context Image)
class ImageCanvas2(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# setup Matplotlib Figure and Axis
		self.fig = Figure(figsize=(3,2.5))
		self.ax = self.fig.add_subplot(111)
		self.ax.set_axis_off()
		# initialization of the canvas
		FigureCanvas.__init__(self, self.fig)
		# we define the widget as expandable
		FigureCanvas.setSizePolicy(self,
								   QtGui.QSizePolicy.Preferred,
								   QtGui.QSizePolicy.Preferred)
		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)
		self.fig.set_facecolor('white')

# Matplotlib Figure object 3 (Histograms)
class ImageCanvas3(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# setup Matplotlib Figure and Axis
		self.fig = Figure(figsize=(3,4))
		self.ax = self.fig.add_subplot(111)
		self.fig.set_tight_layout(True)
		
		# initialization of the canvas
		FigureCanvas.__init__(self, self.fig)
		# we define the widget as expandable
		FigureCanvas.setSizePolicy(self,
							 QtGui.QSizePolicy.Expanding,
								   QtGui.QSizePolicy.Expanding)
#								   QtGui.QSizePolicy.Preferred,
#								   QtGui.QSizePolicy.Preferred)
									 
		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)
		self.fig.set_facecolor('white')
  




class HistoOptions(QtGui.QWidget,Ui_PlotOptions):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		flags1 = QtCore.Qt.Drawer | QtCore.Qt.WindowStaysOnTopHint
		self.setWindowFlags(flags1)
		self.setupUi(self)
		self.pushButtonOkPlot.clicked.connect(self.ReturnPlotOptions)
	
	def ReturnPlotOptions(self):
		self.close()




class FitOptions(QtGui.QWidget,Ui_FormFunctions):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		flags2 = QtCore.Qt.Drawer | QtCore.Qt.WindowStaysOnTopHint
		self.setWindowFlags(flags2)
		self.setupUi(self)
		self.fitpushButton.clicked.connect(self.fitOptions)
	
	def fitOptions(self):
		self.close()




# Main class
class AdvancedWidget(Ui_FlowsMainWindow):
	def __init__(self,appref):
		self.ParentWindow = QtGui.QMainWindow()
		self.setupUi(self.ParentWindow)
		self.ParentWindow.show()

		
		self.xaxes = "x-arcsecs"
		self.yaxes = "y-arcsecs"
		self.menu1 = QtGui.QMenu()
		self.menu2 = QtGui.QMenu()
		##### Connecting maks velocities button with function Flow tab
		self.menu1.addAction('Mask with Velocities', self.Mask_with_Velocities)
		self.menu1.addAction('Mask with Intensities', self.Mask_with_Intensities)
		self.maskPushButton.setMenu(self.menu1)
		##### Menu to choose colors in flow menus
		self.menu3 = QtGui.QMenu()
		self.menu3.addAction('Color 1', self.select_colorHI)
		self.menu3.addAction('Color 2', self.select_colorLo)
		self.ColorSelector.setMenu(self.menu3)
		##### ButtonMenu for histograms
		self.rectbutton = self.menu2.addAction('Rectangle Region')#, self.Region)
		QtCore.QObject.connect(self.rectbutton,QtCore.SIGNAL('triggered()'),self.Region)
		self.RegionpushButton_2.setMenu(self.menu2)
		
		
		
		self.propertiesButton.clicked.connect(self.ShowPlotOptions)
		self.popWindowPlot = HistoOptions()
		
		
		self.fitButton.clicked.connect(self.ShowFitOptions)
		self.popWindowFit = FitOptions()
		
		
######### Menu cmap 1    ########################        
		self.menuCMAP1 = QtGui.QMenu()
		self.maps1 = sorted(m for m in plt.cm.datad)
		for item1 in self.maps1:
			entry1 = self.menuCMAP1.addAction(item1)
			QtCore.QObject.connect(entry1,QtCore.SIGNAL('triggered()'), lambda item1 = item1: self.get_cmap(item1))
			self.cmaptoolButton.setMenu(self.menuCMAP1)        

######### Menu cmap 2    ########################        
		self.menuCMAP2 = QtGui.QMenu()
		self.maps2 = sorted(p for p in plt.cm.datad)
		for item2 in self.maps2:
			entry2 = self.menuCMAP2.addAction(item2)
			QtCore.QObject.connect(entry2,QtCore.SIGNAL('triggered()'), lambda item2=item2: self.get_cmap(item2))
			self.CMAP_contourtoolButton.setMenu(self.menuCMAP2)    
		
		
		
		
		
		
		
		######## Conection signals to open tabs in workspace 3 (Operations) ###########
		openFlows = self.actionFlows
		openFlows.triggered.connect(self.Open_Flows_operation)
		
		openHistogram = self.actionHistogram
		openHistogram.triggered.connect(self.Open_Histogram_operation)
		
		openDivergences = self.actionDivergences
		openDivergences.triggered.connect(self.Open_Divergences_operation)
		
		openOverlays = self.actionOverlays
		openOverlays.triggered.connect(self.Open_Overlays_operation)
		
		openLasso = self.actionLasso
		openLasso.triggered.connect(self.Open_Lasso_operation)
		
		openMultiDim = self.actionMultiDim
		openMultiDim.triggered.connect(self.Open_MultiDim_operation)
		
		
		################################################################################
		
		
		

		######## Conection signals to close tabs in workspace 3 (Operations) ###########
		
		QtCore.QObject.connect(self.closeButton_1, QtCore.SIGNAL("clicked()"), self.close_flows_tab)
		QtCore.QObject.connect(self.closeButton_2, QtCore.SIGNAL("clicked()"), self.close_histograms_tab)
		QtCore.QObject.connect(self.closeButton_3, QtCore.SIGNAL("clicked()"), self.close_divergences_tab)
		QtCore.QObject.connect(self.closeButton_4, QtCore.SIGNAL("clicked()"), self.close_overlays_tab)
		QtCore.QObject.connect(self.closeButton_5, QtCore.SIGNAL("clicked()"), self.close_lasso_tab)
		QtCore.QObject.connect(self.closeButton_6, QtCore.SIGNAL("clicked()"), self.close_multidim_tab)
		
		################################################################################
		QtCore.QObject.connect(self.horizontalSliderLow, QtCore.SIGNAL("sliderReleased()"), self.ContrastTool)
		QtCore.QObject.connect(self.horizontalSliderHigh, QtCore.SIGNAL("sliderReleased()"), self.ContrastTool)
		

		
		
		
				 
#        
		################################################################################
		QtCore.QObject.connect(self.arcradioButton, QtCore.SIGNAL("clicked()"), self.Coords)
		QtCore.QObject.connect(self.meter_radioButton_2, QtCore.SIGNAL("clicked()"), self.Coords)
		
		openCube = self.actionOpenCube
		openCube.triggered.connect(self.cube_do_open)
		
		
		#~ openContext = self.actionOpen_Context_Image
		#~ openContext.triggered.connect(self.context_do_open)
		
		
		self.canvas1=ImageCanvas1()
		self.verticalLayout_7.addWidget(self.canvas1)
		self.widget_navigatorTool= NavigationToolbar(self.canvas1, self.ParentWindow, coordinates=True)
		self.verticalLayout_7.addWidget(self.widget_navigatorTool)
		
		self.canvas2=ImageCanvas2()
		self.verticalLayout_8.insertWidget(0,self.canvas2)
		
		self.canvas3 = ImageCanvas3()
		self.histoverticalLayout_12.insertWidget(1,self.canvas3)
		#~ self.widget_navigatorTool1= NavigationToolbar(self.canvas3, self.ParentWindow, coordinates=True)
		#~ self.histoverticalLayout_12.insertWidget(2,self.widget_navigatorTool1)

		
		#####  Map parameters  #######
		self.FWHMSelector.valueChanged[str].connect(self.select_fwhm)
		self.CadenceSelector.valueChanged[str].connect(self.select_cadence)
		self.PixelsizeSelector.valueChanged[str].connect(self.select_pixel_size)
	   
	   ##### Vector conections  ######
		self.PivotSelector.activated[str].connect(self.select_pivot)
		self.UnitsSelector.activated[str].connect(self.select_units)
		self.WidthSelector.valueChanged[str].connect(self.select_width)
		self.ScaleSelector.valueChanged[str].connect(self.select_scale)
		self.PixelSelector.valueChanged[str].connect(self.select_pix_per_vec)
		
		self.PivotSelector_2.activated[str].connect(self.select_pivot2)
		self.UnitsSelector_2.activated[str].connect(self.select_units2)
		self.WidthSelector_2.valueChanged[str].connect(self.select_width2)
		self.ScaleSelector_2.valueChanged[str].connect(self.select_scale2)
		self.PixelSelector_2.valueChanged[str].connect(self.select_pix_per_vec2)
		self.rgb_mtpl_ini=(1,0,1)
		QtCore.QObject.connect(self.pushButton, QtCore.SIGNAL("clicked()"), self.flow)
		
		#~ QtCore.QObject.connect(self.PlotRegionpushButton, QtCore.SIGNAL("clicked()"), self.plot_region)
		
#        QtCore.QObject.connect(self.ColorSelector, QtCore.SIGNAL("clicked()"), self.select_color)
		QtCore.QObject.connect(self.fullImageButton, QtCore.SIGNAL("clicked()"), self.plot_full)
		
		########## Connect functions of histograms   ################
		QtCore.QObject.connect(self.fullhistpushButton, QtCore.SIGNAL("clicked()"), self.HistogramFull)
		#~ QtCore.QObject.connect(self.MaskVelpushButton, QtCore.SIGNAL("clicked()"), self.HistoMaskVel)
		#~ QtCore.QObject.connect(self.MaskSSpushButton, QtCore.SIGNAL("clicked()"), self.HistoMaskSS)
		
		
		
		
		
		
		#~ QtCore.QObject.connect(self.ContrastPlotpushButton_2, QtCore.SIGNCutAL("clicked()"), self.CutLevels)
		QtCore.QObject.connect(self.RestorepushButton_2, QtCore.SIGNAL("clicked()"), self.Restore)
		
		######## Connect Functions of OverLays    ###################
		QtCore.QObject.connect(self.chooselow_pushButton_2, QtCore.SIGNAL("clicked()"), self.select_colorLo)
		QtCore.QObject.connect(self.hicolpushButton_3, QtCore.SIGNAL("clicked()"), self.select_colorHI)
		QtCore.QObject.connect(self.OverlaypushButton_4, QtCore.SIGNAL("clicked()"), self.Overlays)
		
		######## Connect Functions of Divergen    ###################
		QtCore.QObject.connect(self.plot_diverpushButton, QtCore.SIGNAL("clicked()"), self.DivergencePlot)
		QtCore.QObject.connect(self.ColorSelector_2, QtCore.SIGNAL("clicked()"), self.select_color)
		
		
		QtCore.QObject.connect(self.start_lassoButton, QtCore.SIGNAL("clicked()"), self.SelectionLasso)
		QtCore.QObject.connect(self.stop_lassoButton, QtCore.SIGNAL("clicked()"), self.disconnect)
		QtCore.QObject.connect(self.LassoMaskpushButton, QtCore.SIGNAL("clicked()"), self.PlotMaskLasso)
		
		QtCore.QObject.connect(self.SaveHist_pushButton, QtCore.SIGNAL("clicked()"), self.SaveHisto)
		
		
		#~ self.header = ["Key","Value"]
		#~ self.HeadertableWidget.setHorizontalHeaderLabels(self.header)

	   
		
		
		############ Connet functions of velocities #########################
		self.vminlineEdit.textChanged[str].connect(self.vel_min)
		self.vmaxlineEdit.textChanged[str].connect(self.vel_max)
		############ Connet functions of velocities #########################
		self.iminlineEdit.textChanged[str].connect(self.int_min)
		self.imaxlineEdit.textChanged[str].connect(self.int_max)
		#~ print dir(self.tableWidget)
		#~ print cualquiera
		
		#####################################################################
		#~ print dir(self.canvas3.fig)
		
		
	######### Functions to close Tabs in workspace 3 (Operation tabs) ####### 
		
	def close_flows_tab(self):
		self.MDToolstabWidget_2.removeTab(self.MDToolstabWidget_2.currentIndex())
		if self.MDToolstabWidget_2.count() == 0:
			self.Message.setText("<html><head/><body><p align=\"center\"><span style=\" font-size:24pt;\">In this workspace will be open </span></p><p align=\"center\"><span style=\" font-size:24pt;\">the avalaible operations</span></p><p align=\"center\"><br/></p></body></html>")
		
	
	def close_histograms_tab(self):
		self.MDToolstabWidget_2.removeTab(self.MDToolstabWidget_2.currentIndex())
		if self.MDToolstabWidget_2.count() == 0:
			self.Message.setText("<html><head/><body><p align=\"center\"><span style=\" font-size:24pt;\">In this workspace will be open </span></p><p align=\"center\"><span style=\" font-size:24pt;\">the avalaible operations</span></p><p align=\"center\"><br/></p></body></html>")
		
	def close_divergences_tab(self):
		self.MDToolstabWidget_2.removeTab(self.MDToolstabWidget_2.currentIndex())
		if self.MDToolstabWidget_2.count() == 0:
			self.Message.setText("<html><head/><body><p align=\"center\"><span style=\" font-size:24pt;\">In this workspace will be open </span></p><p align=\"center\"><span style=\" font-size:24pt;\">the avalaible operations</span></p><p align=\"center\"><br/></p></body></html>")
		
	def close_overlays_tab(self):
		self.MDToolstabWidget_2.removeTab(self.MDToolstabWidget_2.currentIndex())
		if self.MDToolstabWidget_2.count() == 0:
			self.Message.setText("<html><head/><body><p align=\"center\"><span style=\" font-size:24pt;\">In this workspace will be open </span></p><p align=\"center\"><span style=\" font-size:24pt;\">the avalaible operations</span></p><p align=\"center\"><br/></p></body></html>")
		
	def close_lasso_tab(self):
		self.MDToolstabWidget_2.removeTab(self.MDToolstabWidget_2.currentIndex())
		if self.MDToolstabWidget_2.count() == 0:
			self.Message.setText("<html><head/><body><p align=\"center\"><span style=\" font-size:24pt;\">In this workspace will be open </span></p><p align=\"center\"><span style=\" font-size:24pt;\">the avalaible operations</span></p><p align=\"center\"><br/></p></body></html>")
		
	def close_multidim_tab(self):
		self.MDToolstabWidget_2.removeTab(self.MDToolstabWidget_2.currentIndex())
		if self.MDToolstabWidget_2.count() == 0:
			self.Message.setText("<html><head/><body><p align=\"center\"><span style=\" font-size:24pt;\">In this workspace will be open </span></p><p align=\"center\"><span style=\" font-size:24pt;\">the avalaible operations</span></p><p align=\"center\"><br/></p></body></html>")
				   

	#########################################################################
	
	
	
	######### Functions to open Tabs in workspace 3 (Operation tabs) ####### 
	
	def Open_Flows_operation(self):
		self.MDToolstabWidget_2.addTab(self.tab_3, _fromUtf8(""))
		self.MDToolstabWidget_2.setTabText(self.MDToolstabWidget_2.indexOf(self.tab_3),"Flows Tools")
		self.Message.setText(" ")
		
	def Open_Histogram_operation(self):
		self.MDToolstabWidget_2.addTab(self.tab_4, _fromUtf8(""))
		self.MDToolstabWidget_2.setTabText(self.MDToolstabWidget_2.indexOf(self.tab_4),"Histogram Tools")
		self.Message.setText(" ")
		
	def Open_Divergences_operation(self):
		self.MDToolstabWidget_2.addTab(self.tab_5, _fromUtf8(""))
		self.MDToolstabWidget_2.setTabText(self.MDToolstabWidget_2.indexOf(self.tab_5),"Divergences Tools")
		self.Message.setText(" ")
		
	def Open_Overlays_operation(self):
		self.MDToolstabWidget_2.addTab(self.tab_6, _fromUtf8(""))
		self.MDToolstabWidget_2.setTabText(self.MDToolstabWidget_2.indexOf(self.tab_6),"Overlays Tools")
		self.Message.setText(" ")
		
	def Open_Lasso_operation(self):
		self.MDToolstabWidget_2.addTab(self.tab, _fromUtf8(""))
		self.MDToolstabWidget_2.setTabText(self.MDToolstabWidget_2.indexOf(self.tab),"Lasso Tools")
		self.Message.setText(" ")
	
		
	def Open_MultiDim_operation(self):
		self.MDToolstabWidget_2.addTab(self.tab_7, _fromUtf8(""))
		self.MDToolstabWidget_2.setTabText(self.MDToolstabWidget_2.indexOf(self.tab_7),"MultiDim Tools")
		self.Message.setText(" ")
		
	
	
	#########################################################################
		
		
		
		
	def header2table(self,array,qtable):
		qtable.setColumnCount(2)
		qtable.setRowCount(len(array))
		for row in range(len(array)):
			for column in range(2):
				qtable.setItem(row,column,QTableWidgetItem(QString("%1").arg(array[row][column])))
		
		
	
	
	def cube_do_open(self):
		cube_data = QtGui.QFileDialog.getOpenFileName(self.ParentWindow, 'Select a Cube', '/home', 
													  "Fits Files (*.fits *.fit *.fts *fits.gz *fts.gz);; IDL files (*sav *save);; All files (*)")
		
		if os.access(cube_data,os.R_OK):
			filename = str(cube_data)
		if filename[-4:]=='save' or filename[-3:]=='sav':
			self.s = readsav(filename)
			self.cube = self.s.values()[0]
		else:
			hdu = fits.open(filename)
			if len(hdu) == 1:
				self.cube = hdu[0].data
				self.hdr_cube = hdu[0].header
			elif len(hdu) == 2:
				self.cube = hdu[1].data
				self.hdr_cube = hdu[1].header
		
		
		
		self.numLowlabel.setNum(self.cube[self.spinBoxInit.value(),:,:].min())
		self.numHighlabel.setNum(self.cube[self.spinBoxInit.value(),:,:].max())
		
		self.horizontalSliderLow.setMinimum(self.cube[self.spinBoxInit.value(),:,:].min())
		self.horizontalSliderLow.setMaximum((self.cube[self.spinBoxInit.value(),:,:].min()+self.cube[self.spinBoxInit.value(),:,:].max())/2)
		
		self.horizontalSliderHigh.setMinimum((self.cube[self.spinBoxInit.value(),:,:].min()+self.cube[self.spinBoxInit.value(),:,:].max())/2+1)
		self.horizontalSliderHigh.setMaximum(self.cube[self.spinBoxInit.value(),:,:].max())
		
		self.horizontalSliderLow.setValue(self.cube[self.spinBoxInit.value(),:,:].min())
		self.horizontalSliderHigh.setValue(self.cube[self.spinBoxInit.value(),:,:].max())
		
		
		self.image = self.cube[self.spinBoxInit.value(),:,:]
		self.X,self.Y = np.meshgrid(np.arange(self.image.shape[1])*float(self.PixelsizeSelector.value()),np.arange(self.image.shape[0])*float(self.PixelsizeSelector.value()))
		self.canvas1.ax.imshow(self.image.clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_title("Vector Velocity Flows",fontsize=16)
		
		self.canvas2.ax.imshow(self.image.clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray")
		self.canvas2.ax.set_title("Context Image")
		self.xaxes = "x-arcsecs"
		self.yaxes = "y-arcsecs"
		self.canvas1.ax.set_xlabel(self.xaxes, fontsize=14)
		self.canvas1.ax.set_ylabel(self.yaxes, fontsize=14)
		
		self.canvas1.show()
		self.canvas1.draw()
		self.canvas2.show()
		self.canvas2.draw()
		self.spinBoxEnd.setValue(self.cube.shape[0])
		self.xsup_spinBox_4.setValue(self.cube.shape[2])
		self.ysup_spinBox_3.setValue(self.cube.shape[1])
		QtGui.QMessageBox.information(self.ParentWindow, 'First Steps', ''' Firstly, you need to calculate the flow \n field before using visualization tools ''',
			QMessageBox.Ok)
		
		self.header = ["Key","Value"]
		self.data_list = self.hdr_cube.items()
		
		self.header2table(np.array(self.data_list),self.HeadertableWidget)
		
		
		
		
		
		return self.cube,self.image
	
		
	
	
	#~ def context_do_open(self):
		#~ image_data = QtGui.QFileDialog.getOpenFileName(self.ParentWindow, 'Select a Context Image', '/home', 
													  #~ "Fits Files (*.fits *.fit *.fts *fits.gz *fts.gz);; IDL files (*sav *save);; All files (*)")
		
		#~ if os.access(image_data,os.R_OK):
			#~ filename = str(image_data)
		#~ hdu = fits.open(filename,checksum=True)
		#~ image = hdu[0].data
		#~ header_im = hdu[0].header
		#~ if len(header_im) > 9:
			#~ self.InfoplainTextEdit.setPlainText("Name : "+header_im["TELESCOP"]+"-"+header_im["INSTRUME"]+" "+header_im["DATE_OBS"]+"\n\n"+\
			#~ "Dimensions : "+str(header_im["NAXIS1"])+"x"+str(header_im["NAXIS2"])+"\n\n"\
			#~ "Min : "+str(image.min())+"\n"+ "Max : "+str(image.max())+"\n"+ "Pixel Size : "+str(header_im["cdelt1"])+"\n"+"Temporal Cadence : UNKNOWN")
		   
		#~ else:
			#~ self.InfoplainTextEdit.setPlainText("Name : UNKNOWN \n\n"+\
			#~ "Dimensions : "+str(header_im["NAXIS1"])+"x"+str(header_im["NAXIS2"])+"\n\n"\
			#~ "Min : "+str(image.min())+"\n"+ "Max : "+str(image.max())+"\n"+ "Pixel Size : UNKNOWN \n"+"Temporal Cadence : UNKNOWN")

		
		#~ self.canvas2.ax.imshow(image,origin="lower",cmap="gray")
		#~ self.canvas2.ax.set_title("Context Image")
		#~ self.canvas2.show()
		#~ self.canvas2.draw()
	
	###### Function to choose the color ######
	def select_color(self):
		inicial=QtGui.QColor(5,25,200)
		color = QtGui.QColorDialog.getColor(inicial,self.ParentWindow) 
		rgb=color.getRgb()
		self.rgb_mtpl_ini=np.array([rgb[0]/255.,rgb[1]/255.,rgb[2]/255.,rgb[3]/255.])
		pal=QtGui.QPalette(color)
		self.ColorSelector.setPalette(pal)
		return self.rgb_mtpl_ini
		
	def select_colorHI(self):
		inicial=QtGui.QColor(5,25,200)
		color = QtGui.QColorDialog.getColor(inicial,self.ParentWindow) 
		rgb=color.getRgb()
		self.rgb_mtpl_iniHi=np.array([rgb[0]/255.,rgb[1]/255.,rgb[2]/255.,rgb[3]/255.])
		pal=QtGui.QPalette(color)
		self.ColorSelector.setPalette(pal)
		return self.rgb_mtpl_iniHi
		
	def select_colorLo(self):
		inicial=QtGui.QColor(5,25,200)
		color = QtGui.QColorDialog.getColor(inicial,self.ParentWindow) 
		rgb=color.getRgb()
		self.rgb_mtpl_iniLo=np.array([rgb[0]/255.,rgb[1]/255.,rgb[2]/255.,rgb[3]/255.])
		pal=QtGui.QPalette(color)
		self.ColorSelector.setPalette(pal)
		return self.rgb_mtpl_iniLo
		
	######   functions to get parameters   ########
	def select_fwhm(self,fwhm1):
		return fwhm1
		
	def select_pixel_size(self,pixsize1):
		return pixsize1
	
	def select_cadence(self,cadence_time1):
		return cadence_time1
	
	def list_of_variables_to_Flow(self,fwhm2,pixsz2,cadence_time2):
		list1=[fwhm2,pixsz2,cadence_time2]
		return list1
		
  ## Functions to get vectors characteristics  #######
	def select_pivot(self,pivot1):
		return pivot1
		
	def select_units(self,unit1):
		return unit1
	
	def select_width(self,width1):
		return width1
	
	def select_scale(self,scale1):
		return scale1
	def select_pix_per_vec(self,ppv):
		return ppv
	
	#### Flow functions calcule the velocities vector for the cube ######
	def flow(self):
		cube = self.cube[self.spinBoxInit.value():self.spinBoxEnd.value(),:,:]
		list1=self.list_of_variables_to_Flow(float(self.FWHMSelector.value()),
											 float(self.PixelsizeSelector.value()),
											 int(self.CadenceSelector.value()))
		lag=1#input('Introduce the lag between the images to be compared (number of images):')
		fwhm_arcsec = list1[0]#0.5#input('Introduce the fwhm of the window for tracking(arcsec) (Could be 3., 0.8,1.2,0.5,1.):')
		reb=1.#input('Introduce the rebinning factor:')
		pix = list1[1]#0.0544#input('Size of pixel (pix):')
		t_step=list1[2]#input('temporal sampling interval in the time series (seconds):')
		fwhm=fwhm_arcsec/pix
		kmperasec=725#input('Value of kilometers per arcsec:')
		h_m=150#input('mass-flux scale-heigth (November 1989, ApJ,344,494):')
#        char=2.#input('Size of the Char for the histograms:')
		#Calculos derivados de los parametros anteriores
		v_limit=2*reb+reb #cota maxima velocidad en pixeles.
		delta_t=t_step*lag # time-lag in seconds
		factor=pix*kmperasec/delta_t
		
		
		vx,vy=flowmaker_2(cube,lag,fwhm,reb)
		
		vx=vx.clip(min=-v_limit,max=v_limit) ; vy=vy.clip(min=-v_limit,max=v_limit)
		self.vx_kps=vx*factor 	#vx in km/s
		self.vy_kps=vy*factor	#vy in km/s
		self.div= divergen(self.vx_kps,self.vy_kps)
		self.vz_kps=h_m*self.div
		self.mag = np.sqrt(self.vx_kps*self.vx_kps + self.vy_kps*self.vy_kps)
		VH= np.mean(self.mag)
		FWHM=fwhm/fwhm_arcsec
		self.InfoplainTextEdit.appendPlainText("Mean Horizontal Velocity (km/s) : "+str(VH)+"\n"+\
		"Maximum Horizontal Velocity (km/s) : "+str(np.max(self.mag))+"\n"+\
		"Minimum Horizontal Velocity (km/s) : "+str(np.min(self.mag)))
		
		self.vx_kps = self.vx_kps- self.vx_kps.mean()
		self.vy_kps = self.vy_kps- self.vy_kps.mean()
		self.vz_kps = self.vz_kps- self.vz_kps.mean()
		
		self.vminlineEdit.setText(str(np.round(self.mag.min(),2)))
		self.vmaxlineEdit.setText(str(np.round(self.mag.max(),2)))
		self.iminlineEdit.setText(str(int(cube[0,:,:].min())))
		self.imaxlineEdit.setText(str(int(cube[0,:,:].max())))
		self.lineEdit_4.setText(str(self.vz_kps.min()))
		self.lineEdit_5.setText(str(self.vz_kps.max()))
		
		
		return self.vx_kps,self.vy_kps,self.vz_kps,self.mag,self.div
	
		
		
	def list_of_variables_to_Vectors(self,ppv2,pivot2,unit2,width2,scale2):
		list2=[ppv2,pivot2,unit2,width2,scale2]
		return list2
		
	def plot_full(self):
		
		list2=self.list_of_variables_to_Vectors(int(self.PixelSelector.value()),
												str(self.PivotSelector.currentText()),
												str(self.UnitsSelector.currentText()),
												float(self.WidthSelector.value()),
												float(self.ScaleSelector.value()))
		self.canvas1.ax.clear()
		self.canvas1.ax.set_title('Velocity field applying LCT')
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:].clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)
		self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx_kps[::list2[0],::list2[0]],self.vy_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=1./list2[4],width=list2[3])
		self.canvas1.show()
		self.canvas1.draw()
	
	#Plotting just a region
	
	#~ def plot_region(self):
		#~ list2=self.list_of_variables_to_Vectors(int(self.PixelSelector.value()),
												#~ str(self.PivotSelector.currentText()),
												#~ str(self.UnitsSelector.currentText()),
												#~ float(self.WidthSelector.value()),
												#~ float(self.ScaleSelector.value()))
		#~ self.canvas1.ax.clear()
		#~ self.canvas1.ax.set_title('Velocity field applying LCT')
		#~ self.im = self.cube[self.InitspinBox.value(),self.yinf_spinBox_2.value():self.ysup_spinBox_3.value(),self.xinf_spinBox.value():self.xsup_spinBox_4.value()]
		#~ self.canvas1.ax.imshow(self.im,origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		#~ self.canvas1.ax.set_xlabel(self.xaxes)
		#~ self.canvas1.ax.set_ylabel(self.yaxes)
		#~ self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  #~ self.maskVelx1[::list2[0],::list2[0]],self.maskVely1[::list2[0],::list2[0]],
							  #~ pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=1./list2[4],
							  #~ scale_units="xy",width=list2[3])
		#~ if (self.arcradioButton.isChecked() == True):
		#~ self.canvas1.ax.axis([self.xinf_spinBox.value(),self.xsup_spinBox_4.value(),self.yinf_spinBox_2.value(),self.ysup_spinBox_3.value()])
			
		#~ self.canvas1.show()
		#~ self.canvas1.draw()
		
	
	### add menu to maskpushButton ####
	
	###########  Functions to create velocities mask ######################
	def vel_min(self,vmin):
		return vmin
	def vel_max(self,vmax):
		return vmax
	def list_val_vel(self,vmin1,vmax1):
		list3 = [vmin1,vmax1]
		return list3
	
	def Mask_with_Velocities(self):
		list3 = self.list_val_vel(float(self.vminlineEdit.text()),float(self.vmaxlineEdit.text()))
		list2=self.list_of_variables_to_Vectors(int(self.PixelSelector.value()),
												str(self.PivotSelector.currentText()),
												str(self.UnitsSelector.currentText()),
												float(self.WidthSelector.value()),
												float(self.ScaleSelector.value()))
#        
		
		self.canvas1.ax.clear()        
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:].clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)
		self.vx1_kps = np.copy(self.vx_kps)
		self.vy1_kps = np.copy(self.vy_kps)
		
		self.vx1_kps[np.logical_or(self.mag < list3[0],self.mag > list3[1])] = np.nan
		self.vy1_kps[np.logical_or(self.mag < list3[0],self.mag > list3[1])] = np.nan
		
		#~ self.vx1_kps = np.ma.masked_where(self.mag < list3[0],self.vx1_kps)
		#~ self.vy1_kps = np.ma.masked_where(self.mag < list3[0],self.vy1_kps)
		#~ self.vx1_kps = np.ma.masked_where(self.mag > list3[1],self.vx1_kps)
		#~ self.vy1_kps = np.ma.masked_where(self.mag > list3[1],self.vy1_kps)
		
		#~ self.mag1 = np.sqrt(self.vx1_kps*self.vx1_kps + self.vy1_kps*self.vy1_kps)
		self.mag1 = np.sqrt(np.power(self.vx1_kps,2)+np.power(self.vy1_kps,2))
		
		
		self.statisticsplainTextEdit.appendPlainText("Statistics with Mask of Velocities \n"+\
		
		"Mean Horizontal Velocity (km/s) : "+str(np.mean(self.mag1))+"\n"+\
		"Maximum Horizontal Velocity (km/s) : "+str(np.max(self.mag1))+"\n"+\
		"Minimum Horizontal Velocity (km/s) : "+str(np.min(self.mag1)))

		self.canvas1.ax.set_title('Mask Velocity field applying LCT')
		self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx1_kps[::list2[0],::list2[0]],self.vy1_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=1./list2[4],
							  scale_units="xy",width=list2[3])
		
		
		
		self.canvas1.show()
		self.canvas1.draw()
		
		self.maskVelx1 = self.vx1_kps
		self.maskVely1 = self.vy1_kps
		
		
		
		return self.maskVelx1,self.maskVely1,self.mag1
		
		
	def PlotMaskLasso(self):
		
		list3 = self.list_val_vel(float(self.vminlineEdit.text()),float(self.vmaxlineEdit.text()))
		
		list2=self.list_of_variables_to_Vectors(int(self.PixelSelector.value()),
												str(self.PivotSelector.currentText()),
												str(self.UnitsSelector.currentText()),
												float(self.WidthSelector.value()),
												float(self.ScaleSelector.value()))
		
		self.canvas1.ax.clear()        
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:].clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)
		
		
		#~ self.vx1 = self.maskVelx1.copy()
		#~ self.vy1 = self.maskVely1.copy()
		#~ self.M = np.zeros(self.vx1.shape, dtype='bool')
		#~ self.M[np.where(self.new_mask==0.)] = True
		#~ self.vx1 = np.ma.masked_array(self.vx1, mask = self.M)
		#~ self.vy1 = np.ma.masked_array(self.vy1, mask = self.M)
		
		
		
		self.mask = self.new_mask
		
		self.vx1 = np.copy(self.vx_kps)
		self.vy1 = np.copy(self.vy_kps)
		self.vx1 = np.ma.masked_where(self.mag < list3[0],self.vx1)
		self.vy1 = np.ma.masked_where(self.mag < list3[0],self.vy1)
		self.vx1 = np.ma.masked_where(self.mag > list3[1],self.vx1)
		self.vy1 = np.ma.masked_where(self.mag > list3[1],self.vy1)
		
		self.vx1 = np.ma.masked_where(self.mask == 0.,self.vx1)
		self.vy1 = np.ma.masked_where(self.mask == 0.,self.vy1)
		
		
		self.mag2 = np.sqrt(np.power(self.vx1,2) + np.power(self.vy1,2))
		
		self.statisticsplainTextEdit.appendPlainText("Statistics with Mask removing Sunspots regions \n"+\
		
		"Mean Horizontal Velocity (km/s) : "+str(np.mean(self.mag2))+"\n"+\
		"Maximum Horizontal Velocity (km/s) : "+str(np.max(self.mag2))+"\n"+\
		"Minimum Horizontal Velocity (km/s) : "+str(np.min(self.mag2)))
		
		self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx1[::list2[0],::list2[0]],self.vy1[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=1./list2[4],
							  scale_units="xy",width=list2[3])
							  
		
							  
		self.canvas1.show()
		self.canvas1.draw()
		
		return self.vx1,self.vy1,self.mag2


###########  Functions to create intensities mask ######################
	def int_min(self,imin):
		return imin
	def int_max(self,imax):
		return imax
	def list_val_int(self,imin1,imax1):
		list4 = [imin1,imax1]
		return list4
	
	def Mask_with_Intensities(self):
		list4 = self.list_val_int(int(self.iminlineEdit.text()),int(self.imaxlineEdit.text()))
		list2=self.list_of_variables_to_Vectors(int(self.PixelSelector.value()),
												str(self.PivotSelector.currentText()),
												str(self.UnitsSelector.currentText()),
												float(self.WidthSelector.value()),
												float(self.ScaleSelector.value()))
#        
		self.canvas1.ax.clear()
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:].clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)
		self.M2 = np.ma.masked_greater(self.cube[0,:,:],list4[0])
		self.M1 = np.ma.masked_less(self.cube[0,:,:],list4[1])
		self.vx1_kps = np.ma.masked_array(self.vx_kps,mask=self.M1.mask)
		self.vy1_kps = np.ma.masked_array(self.vy_kps,mask=self.M1.mask)
		self.vx2_kps = np.ma.masked_array(self.vx_kps,mask=self.M2.mask)
		self.vy2_kps = np.ma.masked_array(self.vy_kps,mask=self.M2.mask)
		

		self.canvas1.ax.set_title('Mask Velocity field applying LCT')
		self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx1_kps[::list2[0],::list2[0]],self.vy1_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=1./list2[4],width=list2[3])
		self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx2_kps[::list2[0],::list2[0]],self.vy2_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniLo,units=list2[2],scale=1./list2[4],width=list2[3])
		
		self.canvas1.show()
		self.canvas1.draw()
		
		self.maskIntx1 = self.vx1_kps
		self.maskInty1 = self.vy1_kps
		self.maskIntx2 = self.vx2_kps
		self.maskInty2 = self.vy2_kps
		
		return self.maskIntx1, self.maskInty1, self.maskIntx2, self.maskInty2
		

		
 ##################################################################################

	def Region(self):
		rect = Rectangle((0,0),0, 0,fill = False,edgecolor="red",linestyle="dashed")
		self.canvas1.ax.add_patch(rect) 
		def onselect(eclick,erelease):
			
			rect.set_width(erelease.xdata - eclick.xdata)
			rect.set_height(erelease.ydata - eclick.ydata)
			rect.set_xy((eclick.xdata, eclick.ydata))
			self.canvas1.draw()
			self.xpos = [eclick.xdata,erelease.xdata]
			self.ypos = [eclick.ydata,erelease.ydata]
			
			if self.pixelradioButton.isChecked() == True:
				self.data_flat = self.cube[0,min(self.ypos):max(self.ypos),min(self.xpos):max(self.xpos)].flatten()
				self.bins = int(self.lineEdit_3.text())
				self.canvas3.ax.clear()
				self.canvas3.ax.hist(self.data_flat,self.bins,facecolor='b',edgecolor='b',alpha=0.5)
				self.canvas3.ax.grid()
				self.canvas3.show()
				self.canvas3.draw()

			if self.velocityradioButton_2.isChecked() == True:
				self.vel_flat = self.mag[min(self.ypos):max(self.ypos),min(self.xpos):max(self.xpos)].flatten()
				self.bins = int(self.lineEdit_3.text())
				self.canvas3.ax.clear()
				self.canvas3.ax.hist(self.vel_flat,self.bins,facecolor='b',edgecolor='b')
				self.canvas3.ax.grid()
				self.canvas3.show()
				self.canvas3.draw()
				
			
		def toggle_selector(event):
			print(' Key pressed.')
			if event.key in ['Q', 'q'] and toggle_selector.RS.active:
				print(' RectangleSelector deactivated.')
				toggle_selector.RS.set_active(False)
			if event.key in ['A', 'a'] and not toggle_selector.RS.active:
				print(' RectangleSelector activated.')
				toggle_selector.RS.set_active(True)
		
		
		toggle_selector.RS = RectangleSelector(self.canvas1.ax, onselect, drawtype='box', 
										useblit=True,
										rectprops=dict(edgecolor='red',
						fill=False,linestyle='dashed'))
		self.canvas1.mpl_connect('key_press_event', toggle_selector)
		self.canvas1.show()
		   
	
#======================================================================#
#																	   #
#                 Histogram Tools                                      #
#																	   #
#======================================================================#
	
	
	
	def ShowPlotOptions(self):
		self.popWindowPlot.show()
		
	def ShowFitOptions(self):
		self.popWindowFit.show()
	
	
	def HistogramFull(self):
		if self.pixelradioButton.isChecked() == True:
			self.data_flat = self.cube[self.spinBoxInit.value(),:,:].flatten()
			self.bins = int(self.lineEdit_3.text())
			self.canvas3.ax.clear()            
			self.canvas3.ax.hist(self.data_flat,self.bins,facecolor='b',edgecolor='b',alpha=0.5)
			self.canvas3.ax.set_title("Intensity Distribution",fontsize=10)
			self.canvas3.ax.set_ylabel("Frequency",fontsize=9)
			self.canvas3.ax.set_xlabel("Intensity Value (Pixel value)",fontsize=9)
			self.canvas3.ax.xaxis.set_label_coords(0.5, -0.065)
			self.canvas3.ax.tick_params(axis='x', labelsize=6.5)
			self.canvas3.ax.tick_params(axis='y', labelsize=6.5)
			#~ self.canvas3.ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
			#~ self.canvas3.ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
			self.canvas3.ax.grid()
			self.canvas3.show()
			self.canvas3.draw()
			
		if self.velocityradioButton_2.isChecked() == True:
			self.vel_flat = self.mag.flatten()
			self.bins = float(self.lineEdit_3.text())
			self.canvas3.ax.clear()
			self.canvas3.ax.hist(self.vel_flat,int(self.bins),facecolor='b',edgecolor='b', alpha = 0.5)
			self.canvas3.ax.set_title("Velocity Distribution",fontsize=10)
			self.canvas3.ax.set_ylabel("Frequency",fontsize=9)
			self.canvas3.ax.set_xlabel("Velocity Value (km/s)",fontsize=9)
			self.canvas3.ax.xaxis.set_label_coords(0.5, -0.065)
			self.canvas3.ax.tick_params(axis='x', labelsize=6.5)
			self.canvas3.ax.tick_params(axis='y', labelsize=6.5)
			self.canvas3.ax.grid()
			self.canvas3.show()
			self.canvas3.draw()
			
		if (self.LOGcheckBox.isChecked() == True) and (self.pixelradioButton.isChecked() == True):
			self.data_flat = self.cube.flatten()
			self.bins = float(self.lineEdit_3.text())
			self.canvas3.ax.clear()
			self.canvas3.ax.hist(self.data_flat,int(self.bins),facecolor='b',edgecolor='b',log = True, alpha=0.5)
			self.canvas3.ax.set_title("Log Intensity Distribution",fontsize=10)
			self.canvas3.ax.grid()
			self.canvas3.show()
			self.canvas3.draw()
			
		if (self.LOGcheckBox.isChecked() == True) and (self.velocityradioButton_2.isChecked() == True):
			self.vel_flat = self.mag.flatten()
			self.bins = float(self.lineEdit_3.text())
			self.canvas3.ax.clear()
			self.canvas3.ax.hist(self.vel_flat,int(self.bins),facecolor='b',edgecolor='b',log = True, alpha=0.5)
			self.canvas3.ax.set_title("Log Velocity Distribution",fontsize=10)
			self.canvas3.ax.grid()
			self.canvas3.show()
			self.canvas3.draw()
			

	
	def SaveHisto(self):
		filename = QtGui.QFileDialog.getSaveFileName(self.ParentWindow, 'Save Image', os.getenv('HOME'))
		filename = str(filename)
		self.canvas3.fig.savefig(filename,format="png")
			
	
	#~ def CutLevels(self):
		#~ self.canvas1.ax.clear()
		#~ self.canvas1.ax.imshow(self.cube[self.InitspinBox.value(),:,:].clip(min=float(self.LOClineEdit.text()),max = float(self.HIClineEdit.text())),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		#~ self.canvas1.ax.set_xlabel(self.xaxes)
		#~ self.canvas1.ax.set_ylabel(self.yaxes)
		#~ self.canvas1.show()
		#~ self.canvas1.draw()
		
	
	
	def Restore(self):
		self.canvas1.ax.clear()
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:],origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)
		self.canvas1.show()
		self.canvas1.draw()
#        
	
		
	
	def Overlays(self):
		self.hi_value = float(self.lineEdit_10.text())

		self.lo_value = float(self.lineEdit_7.text())
	  
		self.opacity = self.doubleSpinBox_3.value()        
		self.hi_color = self.rgb_mtpl_iniHi
  
		self.lo_color = self.rgb_mtpl_iniLo

		self.canvas1.ax.clear()
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:].clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value()),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)

		
		self.patches1 = []
		self.patches2 = []
		
		self.contours1 = measure.find_contours(self.cube[self.spinBoxInit.value(),:,:], self.hi_value)
		self.contours2 = measure.find_contours(self.cube[self.spinBoxInit.value(),:,:], self.lo_value)
		
		if (self.arcradioButton.isChecked() == True):	
			for p in range(len(self.contours1)):
				self.verts = mpath.Path(np.roll(self.contours1[p]*0.504,1,axis=1))
				self.patch = mpatches.PathPatch(self.verts)
				self.patches1.append(self.patch)

			for p in range(len(self.contours2)):
				self.verts = mpath.Path(np.roll(self.contours2[p]*0.504,1,axis=1))
				self.patch = mpatches.PathPatch(self.verts)
				self.patches2.append(self.patch)
		
		if (self.meter_radioButton_2.isChecked() == True):	
			for p in range(len(self.contours1)):
				self.verts = mpath.Path(np.roll(self.contours1[p]*0.504*725./1000.,1,axis=1))
				self.patch = mpatches.PathPatch(self.verts)
				self.patches1.append(self.patch)

			for p in range(len(self.contours2)):
				self.verts = mpath.Path(np.roll(self.contours2[p]*0.504*725./1000.,1,axis=1))
				self.patch = mpatches.PathPatch(self.verts)
				self.patches2.append(self.patch)
		
		self.collection1 = PatchCollection(self.patches1, facecolor=self.hi_color,edgecolor=self.hi_color,alpha=self.opacity) 
		self.collection2 = PatchCollection(self.patches2, facecolor=self.lo_color,edgecolor=self.lo_color,alpha=self.opacity)
		
		self.canvas1.ax.add_collection(self.collection1)
		self.canvas1.ax.add_collection(self.collection2)
		
		
		self.canvas1.show()
		self.canvas1.draw()
	
	
	def ContrastTool(self):
		self.smax = float(self.numHighlabel.text())
		self.smin = float(self.numLowlabel.text())
		
		self.canvas1.ax.clear()
		im=self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:].clip(min=self.smin,max=self.smax),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)

		
		self.canvas1.show()
		self.canvas1.draw()
		#.clip(min=self.horizontalSliderLow.value(),max=self.horizontalSliderHigh.value())
	
	
	
	
	def Coords(self):
		if (self.arcradioButton.isChecked() == True):
			self.X,self.Y = np.meshgrid(np.arange(self.xinf_spinBox.value(),self.xsup_spinBox_4.value())*float(self.PixelsizeSelector.value()),np.arange(self.yinf_spinBox_2.value(),self.ysup_spinBox_3.value())*float(self.PixelsizeSelector.value()))
			self.xaxes = "x-arcsecs"
			self.yaxes = "y-arcsecs"
		
		if (self.meter_radioButton_2.isChecked() == True):
			self.X,self.Y = np.meshgrid(np.arange(self.xinf_spinBox.value(),self.xsup_spinBox_4.value())*float(self.PixelsizeSelector.value())*725/1000,np.arange(self.yinf_spinBox_2.value(),self.ysup_spinBox_3.value())*float(self.PixelsizeSelector.value())*725/1000)
			self.xaxes = "x-Mm"
			self.yaxes = "y-Mm"
		print(self.xaxes,type(self.xaxes))
		print(self.PixelsizeSelector.value())
		print(self.image.shape[1],self.image.shape[0])
		return self.X,self.Y,self.xaxes,self.yaxes
	   
	def select_pivot2(self,pivot1):
		return pivot1
		
	def select_units2(self,unit1):
		return unit1
	
	def select_width2(self,width1):
		return width1
	
	def select_scale2(self,scale1):
		return scale1
	def select_pix_per_vec2(self,ppv):
		return ppv
	def list_of_variables_to_Vectors2(self,ppv2,pivot2,unit2,width2,scale2):
		list5=[ppv2,pivot2,unit2,width2,scale2]
	def get_cmap(self,item):
		self.cmap = item
		return self.cmap

	def DivergencePlot(self):
		
		self.canvas1.ax.clear()  
		self.canvas1.ax.imshow(self.cube[self.spinBoxInit.value(),:,:],origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		self.canvas1.ax.set_xlabel(self.xaxes)
		self.canvas1.ax.set_ylabel(self.yaxes)
		if self.vecradioButton.isChecked() == True:
			list2=self.list_of_variables_to_Vectors(int(self.PixelSelector_2.value()),
												str(self.PivotSelector_2.currentText()),
												str(self.UnitsSelector_2.currentText()),
												float(self.WidthSelector_2.value()),
												float(self.ScaleSelector_2.value()))
		
			self.alpha = float(self.alphalineEdit.text())
			self.loval = float(self.lineEdit_4.text())
			self.hival = float(self.lineEdit_5.text())
			self.canvas1.ax.imshow((self.vz_kps*float(self.Amplificatorlabel.text())).clip(min = self.loval,max = self.hival),origin="lower",
							   cmap=mpl.cm.get_cmap(self.cmap),extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()),
								alpha = self.alpha)
		
			self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx_kps[::list2[0],::list2[0]],self.vy_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_ini,units=list2[2],scale=1./list2[4],width=list2[3]) 
								  
		elif self.contourradioButton.isChecked() == True:
			print("selected contour")
			self.alpha2 = float(self.doubleSpinBox.value())
			self.vmin = float(self.vamindoubleSpinBox_2.value())
			self.vmax = float(self.doubleSpinBox_2.value())
			print("Linea de origin", self.lineEdit_6.text())
			
			self.origin=self.lineEdit_6.text()
			self.levels = None
			if self.l1_radioButton.isChecked() == True:
				self.levels = [float(self.sc_lineEdit_7.text())]
			elif self.l2_radioButton_2.isChecked() == True:
				list_levels = self.lineEdit_8.text().split(",")
				for i in range(len(list_levels)):
					list_levels[i] = float(list_levels[i]) 
				self.levels = list_levels
			elif self.l3_radioButton_3.isChecked() == True:
				range_levels = self.lineEdit_9.text().split(",")
				for i in range(0,1):
					range_levels = np.linspace(float(range_levels[i]),float(range_levels[i+1]),float(range_levels[i+2]))
				self.levels = range_levels
			CS = self.canvas1.ax.contour(self.X, self.Y, self.div, levels=self.levels,alpha=self.alpha2,
										 origin = self.origin, vmin = self.vmin,vmax = self.vmax,
										 cmap=mpl.cm.get_cmap(self.cmap),extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))

		self.canvas1.show()
		self.canvas1.draw()





	






	def SelectionLasso(self):
		"""Select indices from a matplotlib collection using `LassoSelector`.

		Parameters
		----------
		ax : :class:`~matplotlib.axes.Axes`
			Axes to interact with.

		data : Solar map

		"""
		self.canvas1.ax.clear()
		
		self.im = self.cube[self.spinBoxInit.value(),self.yinf_spinBox_2.value():self.ysup_spinBox_3.value(),self.xinf_spinBox.value():self.xsup_spinBox_4.value()]
		self.data = self.canvas1.ax.imshow(self.im,origin="lower",cmap="gray")#,extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		
		
		self.canvas1.draw_idle()
		
		

		####################### NO borrar, esta sección funciona bien con la imagen de contextoooo ################

		def onselect(verts):

			l = self.data.get_array()
			ny,nx = l.shape
			
			self.arrayTest = np.ones([ny,nx])
			
			for i in verts:

				self.arrayTest[i[::-1]] = 0.
			


			for j in range(ny):
				indx = np.where(self.arrayTest[j,:]==0)[0]

				if len(indx) > 1:
					self.arrayTest[j,indx[0]:indx[-1]]=0
			

			self.data.set_data(self.data.get_array()*self.arrayTest)
			
			
			self.canvas1.draw_idle()
		################################################################################################################
		
		
			
		self.im = self.cube[self.spinBoxInit.value(),self.yinf_spinBox_2.value():self.ysup_spinBox_3.value(),self.xinf_spinBox.value():self.xsup_spinBox_4.value()]
		#~ y,x = self.im.shape
		print(self.im.mean())
		self.data = self.canvas1.ax.imshow(self.im,origin="lower",cmap="gray")#,extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
		
		
		
		
		
		if (self.arcradioButton.isChecked() == True):
			self.ticks11 = ticker.FuncFormatter(lambda x, pos:"{0:g}".format(int(x*0.504)))
			self.ticks12 = ticker.FuncFormatter(lambda y, pos:"{0:g}".format(int(y*0.504)))
			self.canvas1.ax.xaxis.set_major_formatter(self.ticks11)
			self.canvas1.ax.yaxis.set_major_formatter(self.ticks12)
			self.canvas1.ax.set_xlabel("x-arcsecs")
			self.canvas1.ax.set_ylabel("y-arcsecs")
		
		if (self.meter_radioButton_2.isChecked() == True):
			self.ticks11 = ticker.FuncFormatter(lambda x, pos:"{0:g}".format(int(x*0.504*(725./1000.))))
			self.ticks12 = ticker.FuncFormatter(lambda y, pos:"{0:g}".format(int(y*0.504*(725./1000.))))
			self.canvas1.ax.xaxis.set_major_formatter(self.ticks11)
			self.canvas1.ax.yaxis.set_major_formatter(self.ticks12)
			self.canvas1.ax.set_xlabel("x-Mm")
			self.canvas1.ax.set_ylabel("y-Mm")
			
		#~ print (self.data.get_array()).shape    # Original funciona con el imshow
		
		self.lasso = LassoSelector(self.canvas1.ax, onselect,useblit=True,lineprops = dict(color='blue', linestyle='-',
				 linewidth = 2, alpha=0.5))
		self.canvas1.show()

		self.canvas1.draw()
		
		raw_input('Press any key to exit')
		
		
		
		return self.lasso
		

	def disconnect(self):
		self.lasso.disconnect_events()
		#test para adquirir estadisticas
		print((self.data.get_array()).mean())
		self.new_mask = self.data.get_array()
		
		#~ self.canvas1.draw_idle()
		
		
		
		return self.new_mask
		















	

		





		
	
	
	
	

if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	f = AdvancedWidget(app)
	sys.exit(app.exec_())
