# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:11:48 2015

@author: hypnus1803
"""
import pickle as pk
import sys,os
sys.path.append(os.getcwd()+'/'+'MainCodes')
import time
from scipy.io.idl import readsav
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.widgets import  RectangleSelector
from matplotlib.patches import Rectangle
import numpy as np
from flowmaker_2 import *
from correlation import divergen
from astropy.io import fits
from scipy.io import readsav as restore
from PyQt4 import QtCore, QtGui
from Vector_Velocities_with_Sunpy import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

# Matplotlib Figure object 1
from matplotlib.figure import Figure
class ImageCanvas1(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# setup Matplotlib Figure and Axis
		self.fig = Figure()
		self.ax = self.fig.add_subplot(111)

		# initialization of the canvas
		FigureCanvas.__init__(self, self.fig)
		# we define the widget as expandable
		FigureCanvas.setSizePolicy(self,
								   QtGui.QSizePolicy.Preferred,
								   QtGui.QSizePolicy.Preferred)
		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)
		self.fig.set_facecolor('white')

# Matplotlib Figure object 2
class ImageCanvas2(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# setup Matplotlib Figure and Axis
		self.fig = Figure(figsize=(3,2.5))
		self.ax = self.fig.add_subplot(111)

		# initialization of the canvas
		FigureCanvas.__init__(self, self.fig)
		# we define the widget as expandable
		FigureCanvas.setSizePolicy(self,
								   QtGui.QSizePolicy.Preferred,
								   QtGui.QSizePolicy.Preferred)
		# notify the system of updated policy
		FigureCanvas.updateGeometry(self)
		self.fig.set_facecolor('white')

# Matplotlib Figure object 3
class ImageCanvas3(FigureCanvas):
	"""Class to represent the FigureCanvas widget"""
	def __init__(self):
		# setup Matplotlib Figure and Axis
		self.fig = Figure(figsize=(2.5,4))
		self.ax = self.fig.add_subplot(111)

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
#        
        
        QtCore.QObject.connect(self.arcradioButton, QtCore.SIGNAL("clicked()"), self.Coords)
        QtCore.QObject.connect(self.meter_radioButton_2, QtCore.SIGNAL("clicked()"), self.Coords)
        openCube = self.actionOpen_Cube
        openCube.triggered.connect(self.cube_do_open)
        openContext = self.actionOpen_Context_Image
        openContext.triggered.connect(self.context_do_open)
        QtCore.QObject.connect(self.WriteButton, QtCore.SIGNAL("clicked()"), self.WriteData)
        
        self.canvas1=ImageCanvas1()
        self.verticalLayout_7.addWidget(self.canvas1)
        self.widget_navigatorTool= NavigationToolbar(self.canvas1, self.ParentWindow, coordinates=True)
        self.verticalLayout_7.addWidget(self.widget_navigatorTool)
        
        self.canvas2=ImageCanvas2()
        self.verticalLayout_8.insertWidget(0,self.canvas2)
        
        self.canvas3 = ImageCanvas3()
        self.histoverticalLayout_12.insertWidget(1,self.canvas3)
        self.widget_navigatorTool1= NavigationToolbar(self.canvas3, self.ParentWindow, coordinates=True)
        self.histoverticalLayout_12.insertWidget(2,self.widget_navigatorTool1)

        
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
#        QtCore.QObject.connect(self.ColorSelector, QtCore.SIGNAL("clicked()"), self.select_color)
        QtCore.QObject.connect(self.fullImageButton, QtCore.SIGNAL("clicked()"), self.plot_full)
        
        ########## Connect functions of histograms   ################
        QtCore.QObject.connect(self.fullhistpushButton, QtCore.SIGNAL("clicked()"), self.HistogramFull)
        QtCore.QObject.connect(self.ContrastPlotpushButton_2, QtCore.SIGNAL("clicked()"), self.CutLevels)
        QtCore.QObject.connect(self.RestorepushButton_2, QtCore.SIGNAL("clicked()"), self.Restore)
        
        ######## Connect Functions of OverLays    ###################
        QtCore.QObject.connect(self.chooselow_pushButton_2, QtCore.SIGNAL("clicked()"), self.select_colorLo)
        QtCore.QObject.connect(self.hicolpushButton_3, QtCore.SIGNAL("clicked()"), self.select_colorHI)
        QtCore.QObject.connect(self.OverlaypushButton_4, QtCore.SIGNAL("clicked()"), self.Overlays)
        
        ######## Connect Functions of Divergen    ###################
        QtCore.QObject.connect(self.plot_diverpushButton, QtCore.SIGNAL("clicked()"), self.DivergencePlot)
        QtCore.QObject.connect(self.ColorSelector_2, QtCore.SIGNAL("clicked()"), self.select_color)

       
        
        
        ############ Connet functions of velocities #########################
        self.vminlineEdit.textChanged[str].connect(self.vel_min)
        self.vmaxlineEdit.textChanged[str].connect(self.vel_max)
        ############ Connet functions of velocities #########################
        self.iminlineEdit.textChanged[str].connect(self.int_min)
        self.imaxlineEdit.textChanged[str].connect(self.int_max)
#       
        
        #####################################################################
                   
        
    def track_coord(self,event):
        self.xpos = event.xdata
        self.ypos = event.ydata
        return self.xpos,self.ypos
        
    
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
            self.cube = hdu[0].data
            self.hdr_cube = hdu[0].header
        self.image = self.cube[0,:,:]
        self.X,self.Y = np.meshgrid(np.arange(self.image.shape[1])*0.16,np.arange(self.image.shape[0])*0.16)
        self.canvas1.ax.imshow(self.image,origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
        self.canvas1.ax.set_xlabel(self.xaxes)
        self.canvas1.ax.set_ylabel(self.yaxes)
        self.canvas1.show()
        self.canvas1.draw()
        QtGui.QMessageBox.information(self.ParentWindow, 'First Steps', ''' Firstly, you need to calculate the flow \n field before using visualization tools ''',
			QMessageBox.Ok)
        return self.cube,self.image
        
    
    
    def context_do_open(self):
        image_data = QtGui.QFileDialog.getOpenFileName(self.ParentWindow, 'Select a Context Image', '/home', 
                                                      "Fits Files (*.fits *.fit *.fts *fits.gz *fts.gz);; IDL files (*sav *save);; All files (*)")
        
        if os.access(image_data,os.R_OK):
            filename = str(image_data)
        hdu = fits.open(filename)
        image = hdu[0].data
        header_im = hdu[0].header
        if len(header_im) > 9:
            self.InfoplainTextEdit.setPlainText("Name : "+header_im["TELESCOP"]+"-"+header_im["INSTRUME"]+" "+header_im["DATE_OBS"]+"\n\n"+\
            "Dimensions : "+str(header_im["NAXIS1"])+"x"+str(header_im["NAXIS2"])+"\n\n"\
            "Min : "+str(image.min())+"\n"+ "Max : "+str(image.max())+"\n"+ "Pixel Size : "+str(header_im["cdelt1"])+"\n"+"Temporal Cadence : UNKNOWN")
           
        else:
            self.InfoplainTextEdit.setPlainText("Name : UNKNOWN \n\n"+\
            "Dimensions : "+str(header_im["NAXIS1"])+"x"+str(header_im["NAXIS2"])+"\n\n"\
            "Min : "+str(image.min())+"\n"+ "Max : "+str(image.max())+"\n"+ "Pixel Size : UNKNOWN \n"+"Temporal Cadence : UNKNOWN")
    
        self.canvas2.ax.imshow(image,origin="lower",cmap="gray")
        self.canvas2.ax.set_title("Context Image")
        self.canvas2.show()
        self.canvas2.draw()
    
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
        cube = self.cube
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
        
        start = time.clock()
        vx,vy=flowmaker_2(cube,lag,fwhm,reb)
        print (time.clock() - start)
        vx=vx.clip(min=-v_limit,max=v_limit) ; vy=vy.clip(min=-v_limit,max=v_limit)
        self.vx_kps=vx*factor 	#vx in km/s
        self.vy_kps=vy*factor	#vy in km/s
        self.div= divergen(self.vx_kps,self.vy_kps)
        self.vz_kps=h_m*self.div
        self.mag = np.sqrt(self.vx_kps*self.vx_kps + self.vy_kps*self.vy_kps)
        VH= np.mean(self.mag)
        FWHM=fwhm/fwhm_arcsec
        self.InfoplainTextEdit.appendPlainText("Maximum velocity (km/s) : "+str(v_limit*factor)+"\n"+\
        "FWHM of the window for tracking (pixel/arc) : "+str(FWHM)+"\n"+\
        "Mean Horizontal Velocity (km/s) : "+str(VH))
        self.vminlineEdit.setText(str(np.round(self.mag.min(),2)))
        self.vmaxlineEdit.setText(str(np.round(self.mag.max(),2)))
        self.iminlineEdit.setText(str(int(cube[0,:,:].min())))
        self.imaxlineEdit.setText(str(int(cube[0,:,:].max())))
        self.lineEdit_4.setText(str(self.div.min()))
        self.lineEdit_5.setText(str(self.div.max()))
        print self.div.min()
        print self.div.max()
        print self.div.mean()
        return self.vx_kps,self.vy_kps,self.vz_kps,self.mag,self.div
    def WriteData(self):
        filename = QtGui.QFileDialog.getSaveFileName(self.ParentWindow, 'Save File', os.getenv('HOME'))
        struct = {"vx_kps":self.vx_kps,"vy_kps":self.vy_kps,
                       "vz_kps":self.vz_kps,"divergences":self.div}
        output = open(filename,"wb")
        pk.dump(struct,output)
        output.close()
      
        
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
        self.canvas1.ax.imshow(self.cube[0,:,:],origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
        self.canvas1.ax.set_xlabel(self.xaxes)
        self.canvas1.ax.set_ylabel(self.yaxes)
        self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx_kps[::list2[0],::list2[0]],self.vy_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=list2[4],width=list2[3])
        self.canvas1.show()
        self.canvas1.draw()
    
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
        self.canvas1.ax.imshow(self.cube[0,:,:],origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
        self.canvas1.ax.set_xlabel(self.xaxes)
        self.canvas1.ax.set_ylabel(self.yaxes)
        self.vx1_kps = np.copy(self.vx_kps)
        self.vy1_kps = np.copy(self.vy_kps)
        self.vx1_kps = np.ma.masked_where(self.mag < list3[0],self.vx1_kps)
        self.vy1_kps = np.ma.masked_where(self.mag < list3[0],self.vy1_kps)
        self.vx1_kps = np.ma.masked_where(self.mag > list3[1],self.vx1_kps)
        self.vy1_kps = np.ma.masked_where(self.mag > list3[1],self.vy1_kps)

        self.canvas1.ax.set_title('Mask Velocity field applying LCT')
        self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx1_kps[::list2[0],::list2[0]],self.vy1_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=list2[4],width=list2[3])

        
        self.canvas1.show()
        self.canvas1.draw()



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
        self.canvas1.ax.imshow(self.cube[0,:,:],origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
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
							  pivot=list2[1],color=self.rgb_mtpl_iniHi,units=list2[2],scale=list2[4],width=list2[3])
        self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx2_kps[::list2[0],::list2[0]],self.vy2_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_iniLo,units=list2[2],scale=list2[4],width=list2[3])
        
        self.canvas1.show()
        self.canvas1.draw()

       

        
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
                self.bins = float(self.lineEdit_3.text())
                self.canvas3.ax.clear()
                self.canvas3.ax.hist(self.data_flat,self.bins,facecolor='b',edgecolor='b')
                self.canvas3.ax.grid()
                self.canvas3.show()
                self.canvas3.draw()

            if self.velocityradioButton_2.isChecked() == True:
                self.vel_flat = self.mag[min(self.ypos):max(self.ypos),min(self.xpos):max(self.xpos)].flatten()
                self.bins = float(self.lineEdit_3.text())
                self.canvas3.ax.clear()
                self.canvas3.ax.hist(self.vel_flat,self.bins,facecolor='b',edgecolor='b')
                self.canvas3.ax.grid()
                self.canvas3.show()
                self.canvas3.draw()
                
            
        def toggle_selector(event):
            print ' Key pressed.'
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print ' RectangleSelector deactivated.'
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print ' RectangleSelector activated.'
                toggle_selector.RS.set_active(True)
        
        
        toggle_selector.RS = RectangleSelector(self.canvas1.ax, onselect, drawtype='box', 
										useblit=True,
										rectprops=dict(edgecolor='red',
                        fill=False,linestyle='dashed'))
        self.canvas1.mpl_connect('key_press_event', toggle_selector)
        self.canvas1.show()
           
    
    
    def HistogramFull(self):
        if self.pixelradioButton.isChecked() == True:
            self.data_flat = self.cube.flatten()
            self.bins = float(self.lineEdit_3.text())
            self.canvas3.ax.clear()            
            self.canvas3.ax.hist(self.data_flat,self.bins,facecolor='b',edgecolor='b')
            self.canvas3.ax.set_ylabel("Frequency",fontsize=10)
            self.canvas3.ax.set_xlabel("Inensity Value (Pixel value)",fontsize=10)
            self.canvas3.ax.tick_params(axis='x', labelsize=8)
            self.canvas3.ax.tick_params(axis='y', labelsize=8)
            self.canvas3.ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
            self.canvas3.ax.grid()
            self.canvas3.show()
            self.canvas3.draw()
            
        if self.velocityradioButton_2.isChecked() == True:
            self.vel_flat = self.mag.flatten()
            self.bins = float(self.lineEdit_3.text())
            self.canvas3.ax.clear()
            self.canvas3.ax.hist(self.vel_flat,self.bins,facecolor='b',edgecolor='b')
            self.canvas3.ax.grid()
            self.canvas3.show()
            self.canvas3.draw()
            
        if (self.LOGcheckBox.isChecked() == True) and (self.pixelradioButton.isChecked() == True):
            self.data_flat = self.cube.flatten()
            self.bins = float(self.lineEdit_3.text())
            self.canvas3.ax.clear()
            self.canvas3.ax.hist(self.data_flat,self.bins,facecolor='b',edgecolor='b',log = True)
            self.canvas3.ax.grid()
            self.canvas3.show()
            self.canvas3.draw()
            
        if (self.LOGcheckBox.isChecked() == True) and (self.velocityradioButton_2.isChecked() == True):
            self.vel_flat = self.mag.flatten()
            self.bins = float(self.lineEdit_3.text())
            self.canvas3.ax.clear()
            self.canvas3.ax.hist(self.vel_flat,self.bins,facecolor='b',edgecolor='b',log = True)
            self.canvas3.ax.grid()
            self.canvas3.show()
            self.canvas3.draw()
            
    
    def CutLevels(self):
        self.canvas1.ax.clear()
        self.canvas1.ax.imshow(self.cube[0,:,:].clip(min=float(self.LOClineEdit.text()),max = float(self.HIClineEdit.text())),origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
        self.canvas1.ax.set_xlabel(self.xaxes)
        self.canvas1.ax.set_ylabel(self.yaxes)
        self.canvas1.show()
        self.canvas1.draw()
        
    
    
    def Restore(self):
        self.canvas1.ax.clear()
        self.canvas1.ax.imshow(self.cube[0,:,:],origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
        self.canvas1.ax.set_xlabel(self.xaxes)
        self.canvas1.ax.set_ylabel(self.yaxes)
        self.canvas1.show()
        self.canvas1.draw()
#        
    
        
    
    def Overlays(self):
        hi_value = float(self.lineEdit_10.text())

        lo_value = float(self.lineEdit_7.text())
      
        self.opacity = self.doubleSpinBox_3.value()        
        self.hi_color = self.rgb_mtpl_iniHi
  
        self.lo_color = self.rgb_mtpl_iniLo

        self.canvas1.ax.clear()
        self.canvas1.ax.imshow(self.image,origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
        self.canvas1.ax.set_xlabel(self.xaxes)
        self.canvas1.ax.set_ylabel(self.yaxes)

        self.M2 = np.ma.masked_greater(self.image,lo_value)
        self.M1 = np.ma.masked_less(self.image,hi_value)
        self.vx1_kps = np.ma.masked_array(self.vx_kps,mask=self.M1.mask)
        self.vy1_kps = np.ma.masked_array(self.vy_kps,mask=self.M1.mask)
        self.vx2_kps = np.ma.masked_array(self.vx_kps,mask=self.M2.mask)
        self.vy2_kps = np.ma.masked_array(self.vy_kps,mask=self.M2.mask)

        self.canvas1.ax.set_title('Overlay Mask')
        self.canvas1.ax.quiver(self.X[::3,::3],self.Y[::3,::3],self.vx1_kps[::3,::3],self.vy1_kps[::3,::3],
							  pivot="tail",color=self.rgb_mtpl_iniHi,units="dots",scale=15.0,width=1.9,alpha = self.opacity)
        self.canvas1.ax.quiver(self.X[::3,::3],self.Y[::3,::3], self.vx2_kps[::3,::3],self.vy2_kps[::3,::3],
							  pivot="tail",color=self.rgb_mtpl_iniLo,units="dots",scale=15.0,width=1.9,alpha = self.opacity)
        
        self.canvas1.show()
        self.canvas1.draw()
        
    def Coords(self):
        if (self.arcradioButton.isChecked() == True):
            self.X,self.Y = np.meshgrid(np.arange(self.image.shape[1])*float(self.PixelsizeSelector.value()),np.arange(self.image.shape[0])*float(self.PixelsizeSelector.value()))
            self.xaxes = "x-arcsecs"
            self.yaxes = "y-arcsecs"
        
        if (self.meter_radioButton_2.isChecked() == True):
            self.X,self.Y = np.meshgrid(np.arange(self.image.shape[1])*float(self.PixelsizeSelector.value())*725/1000,np.arange(self.image.shape[0])*float(self.PixelsizeSelector.value())*725/1000)
            self.xaxes = "x-Mm"
            self.yaxes = "y-Mm"
        print self.xaxes,type(self.xaxes)
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
        self.canvas1.ax.imshow(self.image,origin="lower",cmap="gray",extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()))
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
            self.canvas1.ax.imshow((self.div*float(self.Amplificatorlabel.text())).clip(min = self.loval,max = self.hival),origin="lower",
                               cmap=mpl.cm.get_cmap(self.cmap),extent=(self.X.min(),self.X.max(),self.Y.min(),self.Y.max()),norm=mc.Normalize(vmin=self.div.min(), vmax=self.div.max()),
                                alpha = self.alpha)
        
            self.canvas1.ax.quiver(self.X[::list2[0],::list2[0]],self.Y[::list2[0],::list2[0]],
							  self.vx_kps[::list2[0],::list2[0]],self.vy_kps[::list2[0],::list2[0]],
							  pivot=list2[1],color=self.rgb_mtpl_ini,units=list2[2],scale=list2[4],width=list2[3])     
        elif self.contourradioButton.isChecked() == True:
            print "selected contour"
            self.alpha2 = float(self.doubleSpinBox.value())
            self.vmin = float(self.vamindoubleSpinBox_2.value())
            self.vmax = float(self.doubleSpinBox_2.value())
            print "Linea de origin", self.lineEdit_6.text()
            
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
    

        
        
        

        
    
    
    
    

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    f = AdvancedWidget(app)
    sys.exit(app.exec_())
