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

atomnumber = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':36,'BR':36,'CL':17}
numberatom = ['H','He','Li','Be','B','C','N','O','F','Ar','Na','Mg','Al','a','P','S']

#=======================================================================

class Fragment:
	
	
	def __init__(self):
		self.name = ''
		self.index = 0
		self.charge = 0
		self.Natoms = 0	
		self.dipole = []				
		self.chargeN2 = []
		self.chargeN3 = []
		self.energy = 0		
	
	def get_energy(self):
		pass
	
	def __sub__(self):
		pass
		
	def get_coord(self):
		pass

#========================================================================

class fmo_parser:
	
	def __init__(self      ,
				logfile    ,
				totCharge=0):

		self.name      = logfile
		self.datfile   = logfile[:-4] + '.dat'		
		self.Energy    = 0 
		self.EnergyS   = 0 
		self.nAtoms    = 0 
		self.nFrag     = 0 
		self.atoms     = []
		self.totCharge = 0		
		self.chargesN2 = []
		self.chargesN3 = []			
		self.Frag      = []
		self.nbody     = 2 
		self.deltQ     = 0
		self.deltAQ    = 0
		self.solvent   = False
		
		file_g = open(self.logfile,'r')

		for line in file_g:
			line2 = line.split()
			if not line2:
				continue
			if len(line2) > 2:
				if line2[1] == "absolute" and line2[3] == "transf.":
					self.deltQ = float(line2[6])
				elif line2[1] == "amount" and line2[4] == "transf.":
					self.deltAQ = float(line2[6])
				elif line2[0] == "Free" and line2[4] == "solvent=":
					self.EnergyS = float(line2[5])
				elif line2[0] == "Total" and line2[1] == "Energy":
					self.Energy = float(line2[3])
	
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
						self.nFrag += 1 
						self.Frag.append(frg)

		prhase3 = "One-body FMO properties."
		prhase4 = "Frontier molecular orbital (FMO!) properties based on Koopmans' theorem."

		fragp_init = 0
		fragp_fin = 0 

		with open(self.name,'r') as text:
			for (i, line) in enumerate(text):
				if phrase3 in line:
					fragp_init=i
				elif phrase4 in line:
					fragp_fin=i	

		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i >= fragp_init and i <= fragp_fin:
					line2 = line.split()					
					if len(line2) == 5:
						self.Frag[i].energy = float(line2[1])
					

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
		elif self.nbody == 3:
			lline = 6
			
		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i >= chg_init and chg_fin <= stat_fin:
					line2 = line.split()					
					if len(line2) == lline:						
						a = pdb_atom()
						a.num = int(line2[0])
						a.resNum = int(line2[1])
						a.element = numberatom[int(line2[2])]
						a.charge = float(line2[4])	
						self.atoms.append(a)				
						self.ChargesN2 = float(line2[4])
						self.nAtoms += 1
					elif len(line2) == lline:
						a = pdb_atom()
						a.num = line2[0]
						a.resNum = line2[1]
						a.element = numberatom[int(line2[2])]						
						a.charge = float(line2[5])
						self.atoms.append(a)
						self.nAtoms += 1
						self.ChargesN2 = float(line2[4])
						self.chargesN3 = float(line2[5])

#======================================================================

class global_rd:

	def __init__(self  ,
				 neutro,
				 cation,
				 anion):

		self.neutro = neutro
		self.cation = cation
		self.anion = anion 
		self.hardness = 0.00
		self.softness = 0.00
		self.electrophilicity = 0.00 
		self.IP = 0.00
		self.EA = 0.00
		self.chem_pot = 0.00 
		self.frag_hardness = []
		self.frag_softness = []
		self.frag_elect = [] 
		self.frag_chempot = []


	def calc_finite(self):


		energy_neutro = 0
		energy_cation = 0
		energy_anion  = 0

		if neutro.solvent == True:		
			energy_neutro = self.neutro.EnergyS
			energy_cation = self.cation.EnergyS
			energy_anion  = self.anion,EnergyS
		elif neutro.solvent == False:
			energy_neutro = self.neutro.Energy
			energy_cation = self.cation.Energy
			energy_anion  = self.anion.Energy


		self.IP = energy_cation - energy_neutro
		self.EA = energy_neutro - energy_anion
		self.chem_pot = energy_anion - energy_cation
		self.hardness = (-self.EA) - (-self.IP)
		self.softness = 1/self.hardness
		self.electrophilicity = (self.softness**2)*self.chem_pot

	def rd_frag(self):

		for i in range(self.neutro.nFrag):
			frag_hardness.append((self.anion.Frag[i].energy+self.cation.Frag[i].energy - 2*self.neutro.Frag[i].energy)/2)
			frag_softness.append(1/((self.anion.Frag[i].energy+self.cation.Frag[i].energy - 2*self.neutro.Frag[i].energy)/2))
			frag_chempot.append((self.anion.Frag[i].energy-self.cation.Frag[i].energy)/2)
			frag_elect.append((frag.softness[i]**2)*frag_chempot[i])

	def write_glob(self,
				   inpname):

		glob_file = open(inpname,'w')

		glob_text = ''
		glob_file.write(glob_text)


#====================================================================

class local_rd:

	def __init__(self
				neutro,
				cation, 
				anion):

		self.globRD = None
		self.fukuiES = []
		self.fukuiNS = []
		self.fukuiRS = []
		self.deltFukui = []
		self.softnesMu = []
		self.electMU = []
		self.relativeES = []
		self.relativeNS = [] 

		self.globRD = global_rd(self.neutro,self.cation,self.anion)
		self.globRD.calc_finite()

	def cond_Fukui(self):

		for i in range(self.neutro.nAtoms):
			self.fukuiES.append( self.anion.atoms.charge[i] - self.neutro.atoms.charge[i] )
			self.fukuiNS.append( self.neutro.atoms.charge[i] - self.cation.atoms.charge[i] )
			self.fukuiRS.append( (self.anion.atoms.charge[i] - self.cation.atoms.charge[i])/2 )
			self.deltFukui.append( self.fukuiES[i] - self.fukuiNS[i] )
			self.softnesMu.append( self.globRD.softness*self.deltFukui[i])
			self.electMU.append( self.globRD.electrophilicity*self.deltFukui[i])
			self.relativeNS.append ( self.fukuiNS[i]/self.fukuiES[i] )
			self.relativeES.append ( self.fukuiES[i]/self.fukuiNS[i] )

	def writeLRD(self   , 
				 inpman):
		pass





