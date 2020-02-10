#/usr/bin/env python
# -*- coding: utf-8 -*-
#md_analysis.py


import os,glob,sys
import mdtraj as md 
import numpy as np 

trajs = glob.glob("*.xtc")

trj_obj = []

#=======================================================================
class md_analysis:
	def __init__(self):
		
		self.trajs   = glob.glob("*.xtc")
		self.trj_obj = []
		self.names = []
		
	#-------------------------------------------------------------------
	
	def load_trajs(self):
	
		for i in range(len(trajs)):
			t = md.load(trajs[i],top=trajs[i][:-4]+".gro")
			self.trj_obj.append(t)
			self.names.append(trajs[i][:4])
			
	#-------------------------------------------------------------------
	
	def plot_rmsd_rg(self):
		
		rmsd = []
		rg   = []
		for i in range(len(self.trj_obj)):
			rm = md.rmsd(self.trj_obj[i],self.trj_obj[i])
			RG = md.compute_rg(self.trj_obj[i])
			rmsd.append(rm)
			rg.append(RG)
		
		r_data     = "time "
		r_rmsd_txt = open("trj_rmsd_r",'w')
		r_rg_txt = open("trj_rg_r",'w')
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
			
		r_rmsd_txt.write(r_data)
		
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

		
		
		
