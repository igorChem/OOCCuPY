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
		y = self.rg
		x = self.rmsd

        # Calculate the point density
		xy = np.vstack([x,y])
		z = gaussian_kde(xy)(xy)
		
        # Sort the points by density, so that the densest points are plotted last
		idx = z.argsort()
		x, y, z = x[idx], y[idx], z[idx]

        # Setting plot type 
		pdf = ax1.scatter(x, y, c = z, s = 50, edgecolor = '')

        # Plot title
		#ax1.set_title('RG' + ' by ' + 'RMSD')

        # Hide right and top spines
		ax1.spines['right'].set_visible(False)
		ax1.spines['top'].set_visible(False)
		ax1.yaxis.set_ticks_position('left')
		ax1.xaxis.set_ticks_position('bottom')
						
        # Set x and y limits
		xmin = x.min() 
		xmax = x.mean() + 1
		ymin = y.min() -0.05
		ymax = y.max() +0.05       
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)

        # Set x and y labels
		plt.ylabel("RG")
		plt.xlabel("RMSD")

        # Adding the color bar 
		colbar = plt.colorbar(pdf)
		colbar.set_label('Probability Density Function')   
		
		#printing varible stats 
		print("printing mean of RMSD: " + str( x.mean() ) )
		print("printing mean of RG: " + str( y.mean() ) )		
		
		plt.show()
		
	def get_frame(self):
		
		mean_rmsd  = np.mean(self.rmsd)
		mean_rg    = np.mean(self.rg)

		diff_1     = 0 
		diff_2     = 0
		diff_total = []

		for i in range(len(self.rmsd)):
			diff_1 = (self.rmsd[i] - mean_rmsd)**2 
			diff_2 = (self.rg[i] - mean_rg)**2 
			a = diff_1 + diff_2
			diff_total.append(a)
		
		chosen_frame = 0	
		
		for i in range(1,len(diff_total)):			
			if diff_total[i] < diff_total[i-1]:
				chosen_frame = i 
		
		print("chosen frame: "+str(chosen_frame))
		print("RMSD: "+ str(self.rmsd[chosen_frame]))
		print("RG: "  + str(self.rg[chosen_frame]))
		
		ls_ = glob.glob("*.pdb")
		if self.name[:-4]+"_trj.pdb" in ls_:
			print("trj ok")     			
		else:
			os.system("/usr/bin/gmx trjconv -f "+ self.name + " -s " + self.name[:-4]+".tpr -pbc mol -dt 10 -o "+ self.name[:-4]+"_trj.pdb")
		modnum = 0
		j      = 0
		init   = 0
		fin    = 0
		pdb_file = ""
		pdb_out  = open(self.name[:-4] +"_"+str(chosen_frame) + ".pdb",'w')
		
		with open(self.name[:-4]+"_trj.pdb") as file_object:
			for line in file_object:
				if not init > 0:
					line2 = line.split()
					if line2[0] == "MODEL" and line2[1] == str(chosen_frame+1):
						init = j
				elif init > 0:
					line2 = line.split()
					pdb_file +=line
					if line2[0] == "TER":
						break
				j+=1
		
		pdb_out.write(pdb_file)
		pdb_out.close()
			
					
				
			
			
		
		
		
