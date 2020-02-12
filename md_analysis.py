#/usr/bin/env python
# -*- coding: utf-8 -*-
#md_analysis.py


import os,glob,sys

import mdtraj as md
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

trajs = glob.glob("*.xtc")

trj_obj = []

#=======================================================================
class md_analysis:
	def __init__(self,name,tpl):
		
		self.trj_obj = md.load(name,top=tpl)
		self.rg      = 0
		self.rmsd    = 0
		self.name    = name
		self.time    = 0
			
	#-------------------------------------------------------------------
	
	def get_rmsd_rg(self):
		self.rg   = md.compute_rg(self.trj_obj)
		self.rmsd = md.rmsd(self.trj_obj,self.trj_obj)
		self.time = self.trj_obj.time 
	
	#-------------------------------------------------------------------
	
	def write_data(self):
		
		r_data     = "time "
		r_data    += self.name + " \n"
		r_rmsd_txt = open("rmsd",'w')
		r_rg_txt   = open("trj_rg_r",'w')
		r_script   =""
		
		for i in range(len(self.rmsd)):
			r_data += str(self.time[i]) + "  " + str(self.rmsd[i]) +"\n"
			
		r_rmsd_txt.write(r_data)
		r_rmsd_txt.close()
		
		r_data     = "time "
		r_data    += self.name + " \n"
		
		for i in range(len(self.rg)):
			r_data += str(self.time[i]) + "  " + str(self.rg[i]) +"\n"
			
		r_rg_txt.write(r_data)
		r_rg_txt.close()

	#-------------------------------------------------------------------
	
	def plot_rmsd_rg(self):
		
        # Color by the Probability Density Function. 
        # Kernel density estimation is a way to estimate 
        # the probability density function (PDF) of a random 
        # variable in a non-parametric way
		fig, (ax1) = plt.subplots(nrows=1)

        # Setting data
		x = self.rg
		y = self.rmsd

        # Calculate the point density
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		
        # Sort the points by density, so that the densest points are plotted last
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]

        # Setting plot type 
		pdf = ax1.scatter(x, y, c = z, s = 50, edgecolor = '')

        # Plot title
		ax1.set_title('RG' + ' by ' + 'RMSD')

        # Hide right and top spines
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		ax1.yaxis.set_ticks_position('left')
		ax1.xaxis.set_ticks_position('bottom')
						
        # Set x and y limits
		xmin = x.min()
		xmax = x.max()
		ymin = y.mean() -0.5
		ymax = y.mean() + 0.1        
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)

        # Set x and y labels
		plt.xlabel("RG")
		plt.ylabel("RMSD")

        # Adding the color bar 
		colbar = plt.colorbar(pdf)
		colbar.set_label('Probability Density Function')     
		plt.show()
		
		
		
