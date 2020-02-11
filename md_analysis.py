#/usr/bin/env python
# -*- coding: utf-8 -*-
#md_analysis.py


import os,glob,sys

try:
	import mdtraj as md
except:
	pass
	 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

trajs = glob.glob("*.xtc")

trj_obj = []

#=======================================================================
class md_analysis:
	def __init__(self,name,tpl):
		
		self.trj_obj = md.load(traj,top=tpl)
		self.rg      = 0
		self.rmsd    = 0
		self.name    = name
		self.time    = 0
			
	#-------------------------------------------------------------------
	
	def get_rmsd_rg(self):
		self.rg   = md.calculate_rg(self.trj_obj)
		self.rmsd = md.rmsd(self.trj_obj,self.trj_obj)
		self.time = self.trj_obj.time 
	
	#-------------------------------------------------------------------
	
	def write_rmsd(self):
		
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
		pass
