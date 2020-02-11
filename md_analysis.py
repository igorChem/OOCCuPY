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
		self.rg      = []
		self.rmsd    = []
			
	#-------------------------------------------------------------------
	
	def get_rmsd_rg(self):
	

	def write_rmsd(self):
		
		r_data     = "time "
		r_rmsd_txt = open("rmsd",'w')
		r_rg_txt   = open("trj_rg_r",'w')
		r_script   =""
		
		for j in range(len(self.trj_obj)):	
			r_data +=self.names[j]+" "
		r_data +="\n"
		
		for i in range(len(self.trj_obj[0])):
			r_data +=str(self.trj_obj[0].time[i]) +" "
			for j in range(len(self.trj_obj)):			
				r_data += str(rmsd[j][i]) +" "
			r_data += "\n"
		
		print("rmsd tables ok")
		input()
		
		r_rmsd_txt.write(r_data)
		r_rmsd_txt.close()
		
		r_data = "time "
		
		for j in range(len(self.trj_obj)):	
			r_data +=self.names[j]+" "
		r_data +="\n"
		
		for i in range(len(self.trj_obj[0])):
			r_data +=str(self.trj_obj[0].time[i]) +" "
			for j in range(len(self.trj_obj)):			
				r_data += str(rg[j][i]) +" "
			r_data += "\n"
			
		
		r_rg_txt.write(r_data)
		r_rg_txt.close()
		
