#!/usr/bin/env python
# -*- coding: utf-8 -*-
# fmoparser.py

#=======================================================================

import os 
import glob 
import numpy as np
from cube_class import * 
from pdb_class import *

#=======================================================================
class Fragment:
	
	
	def __init__(self):
		self.name = ''
		self.index = 0
		self.charge = 0
		self.Natoms = 0		
		self.chargeN1 = []		
		self.chargeN2 = []
		self.chargeN3 = []
		self.atomsN = []
		self.xcoord = []
		self.ycoord = []
		self.zcoord = []
		self.energy = []
		
	
	def get_energy(self):
		pass
	
	def __sub__(self):
		pass
		
	def get_coord(self):
		pass

class fmo_parser:
	
	def __init__(self,logfile,totCharge=0):
		self.name = logfile
		self.datfile = logfile[:-4] + '.dat'
		self.HOMOn = 0
		self.LUMOn = 0
		self.Energy = 0 
		self.totCharge = 0
		self.chargesN1 = []
		self.chargesN2 = []
		self.chargesN3 = []			
		self.Frag = []
		self.nbody = 2 
		
	def get_front_orb(self):
		
		log = open(self.name,'r')
		
		for line in log: 
			line = line.split()
			if len(line) == 8:
				if line[1] == 'ORBITALS' and line[3] == 'OCCUPIED':
					self.HOMOn += int(line[0])
		self.LUMOn = self.HOMOn + 1
		
		log.close()
	
	def cube_from_dat(self):
				
		phrase1 = 'DENSITY: full system, created by GAMESS (FMO).'
		phrase2 = '$END'
		phrase3 = 'GAMESS CUBE FORMAT: ELECTRON DENSITY'
					
		ir_init = 0
		ir_fin  = 0
		cube_info =''
		
		with open(self.datfile,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					ir_init=i
				elif phrase3 in line:
					ir_init=i					
				elif phrase2 in line:
					ir_fin=i
					
		with open(self.datfile,'r') as text:
			for (i,line) in enumerate(text):
				if i >= ir_init and i <=ir_fin :
					cube_info += line
					 
		cube_file0 = open(self.name[:-4]+'.cube','w')
		cube_file0.write(cube_info)
		cube_file0.close()
		
	def get_frag_stat(self):
		
		phrase1 = 'Fragment statistics'
		phrase2 = 'CPU     0: STEP CPU TIME=     0.12 TOTAL CPU TIME=       63.2 (    1.1 MIN)'
		
		stat_init = 0
		stat_fin  = 0
		
		with open(self.name,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					stat_init=i					
				elif phrase2 in line:
					stat_fin=i
					

					
		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i >= stat_init and i <= stat_fin:
					line2 = line.split()
					if len(line2) == 14:						
						frg = Fragment()
						frg.name   = line2[1]
						frg.index  = line2[0]
						frg.charge = line2[2]
						frg.Natoms = line2[3]
						self.Frag.append(frg)					
					
	def get_charges(self):
		
		phrase1 = 'IAT  IFG   Z       Q(1)        Q(2)        Q(3)'
		phrase2 = 'Done with FMO properties.'
		
		chg_init = 0
		chg_fin  = 0
		
		with open(self.name,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					chg_init=i					
				elif phrase2 in line:
					chg_fin=i
					
		print(chg_init,chg_fin,i)
					
		if self.nbody == 2:
			lline = 5
		elif self.nbody ==3:
			lline = 6
			
		IAT = []
		IFG = []	
		
					
		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i >= chg_init and chg_fin <= stat_fin:
					line2 = line.split()					
					if len(line2) == lline:						
						IAT.append(line2[0])
						IFG.append(line2[1])
						self.ChargesN1 = line2[3]
						self.ChargesN2 = line2[4]
					elif len(line2) == lline:
						IAT.append(line2[0])
						IFG.append(line2[1])
						self.ChargesN1 = line2[3]
						self.ChargesN2 = line2[4]
						self.chargesN3 = line2[5]

def cond_Fukui(file1,file2,chg):
	pass

def cube_Fukui_elec(file1,file2,chg):
	
	prot_neutro = Cube(Typ='Elec_Den')
	prot_anion = Cube(Typ='Elec_Den')
	
	prot_neutro.Read_Elec_Cube(file1)
	prot_anion.Read_Elec_Cube(file2)
		
	prot_neutro.scalar3d = np.absolute((-prot_neutro.scalar3d) - (-prot_anion.scalar3d))
	prot_neutro.scalar3d = prot_neutro.scalar3d/chg
	prot_neutro.write_cubeElec(file1[:-5]+'Fukui_elec')

'''	
def cube_Fukui_Nuc(file1,file2,chg):
	
	prot_neutro = Cube(Typ='Elec_Den')
	prot_cation = Cube(Typ='Elec_Den')
	
	prot_neutro.Read_Elec_Cube(file1)
	prot_anion.Read_Elec_Cube(file2)
'''	
	
	
#cube_Fukui_elec('rta_Fukui0.cube','rta_Fukui20.cube',4)	



'''
for i in range(len(rta.Frag)):
	print(rta.Frag[i].name)


rta = 
rta.get_front_orb()
rta.cube_from_dat()
'''


