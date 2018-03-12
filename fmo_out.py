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
		
		file_g = open(self.name,'r')

		for line in file_g:
			line2 = line.split()
			if not line2:
				continue
			if len(line2) > 2:
				if line2[1] == "absolute" and line2[3] == "transf.":
					self.deltQ = float(line2[6])
				elif line2[1] == "amount" and line2[4] == "transf.":
					self.deltAQ = float(line2[7])
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
		phrase2 = ' Close fragment pairs, distance relative to vdW radii'
		
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

		phrase3 = "One-body FMO properties."
		phrase4 = "Frontier molecular orbital (FMO!) properties based on Koopmans' theorem."

		fragp_init = 0
		fragp_fin = 0 

		with open(self.name,'r') as text:
			for (i, line) in enumerate(text):
				if phrase3 in line:
					fragp_init=i+4 
				elif phrase4 in line:
					fragp_fin=i	
		
		energies = []
		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i >= fragp_init and i <= fragp_fin:
					line2 = line.split()
					try:
						if len(line2) == 5:
							energies.append( float(line2[1]) )
						elif len(line2) == 4:
							energies.append( float(line2[1]) )
					except:
						if len(line2) == 5: 
							energies.append( float(line2[2]) )
						elif len(line2) == 4:
							energies.append( float(line2[1]) )	

		print(self.nFrag,len(energies))
		for i in range(self.nFrag):			
			self.Frag[i].energy=energies[i]
						
					

	def get_charges(self):
		

		phrase1 =''
		if self.solvent == True:
			phrase1 = 'IAT  IFG   Z  surface cover,%   q(ASC)       Q(1)        Q(2)        Q(3)'
		else:
			phrase1 = 'IAT  IFG   Z       Q(1)        Q(2)        Q(3)'

		phrase2 = 'Done with FMO properties.'
		phrase3 = 'n-body  Mulliken atomic spin populations S(n)'
		
		chg_init = 0
		chg_fin  = 0
		chg_fin2 = 0 
		
		with open(self.name,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					chg_init=i										
				elif phrase2 in line:
					chg_fin=i
				elif phrase3 in line:
					chg_fin2=i

		if chg_fin2 > 0 and chg_fin2 < chg_fin:
			chg_fin = chg_fin2

		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i >= chg_init+1 and i <= chg_fin:
					line2 = line.split()					
					if len(line2) == 5 and self.solvent == False:																	
						a = pdb_atom()
						a.num = int(line2[0])
						a.resNum = int(line2[1])
						a.element = numberatom[int(float(line2[2]))-1]
						a.charge = float(line2[4])							
						self.atoms.append(a)				
						self.ChargesN2 = float(line2[4])
						self.nAtoms += 1
						#print(len(line2),i)
					elif len(line2) == 6 and self.solvent == False:
						a = pdb_atom()
						a.num = line2[0]
						a.resNum = line2[1]
						a.element = numberatom[int(line2[2])-1]						
						a.charge = float(line2[5])
						self.atoms.append(a)
						self.nAtoms += 1
						self.ChargesN2 = float(line2[4])
						self.chargesN3 = float(line2[5])
						#print(len(line2),i)
					elif len(line2) == 8:													
						a = pdb_atom()
						a.num = int(line2[0])
						a.resNum = int(line2[1])
						a.element = numberatom[int(float(line2[2]))-1]
						a.charge = float(line2[7])							
						self.atoms.append(a)				
						self.ChargesN2 = float(line2[7])
						self.nAtoms += 1
						#print(len(line2),i)
					elif len(line2) == 9:																	
						a = pdb_atom()
						a.num = int(line2[0])
						a.resNum = int(line2[1])
						a.element = numberatom[int(float(line2[2]))-1]
						a.charge = float(line2[8])							
						self.atoms.append(a)				
						self.ChargesN2 = float(line2[7])
						self.chargesN3 = float(line2[8])
						self.nAtoms += 1
						#print(len(line2),i)
						

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

		if self.neutro.solvent == True:		
			energy_neutro = self.neutro.EnergyS
			energy_cation = self.cation.EnergyS
			energy_anion  = self.anion.EnergyS
		elif self.neutro.solvent == False:
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
			self.frag_hardness.append(((self.anion.Frag[i].energy+self.cation.Frag[i].energy) - 2*self.neutro.Frag[i].energy)/2)
			self.frag_softness.append(1/((self.anion.Frag[i].energy+self.cation.Frag[i].energy - 2*self.neutro.Frag[i].energy)/2))
			self.frag_chempot.append((self.anion.Frag[i].energy-self.cation.Frag[i].energy)/2)
			self.frag_elect.append((self.softness**2)*self.frag_chempot[i])

	def write_glob(self,
				   inpname):

		glob_file = open(inpname,'w')

		glob_text = 'name of the logfile for neutro state {0} \n'.format(self.neutro.name)
		glob_text += 'Info about the different charge states calculations\n'
		glob_text += 'netro state info: \n'
		glob_text += '------\n'
		glob_text += 'Total energy: {0} | Total energy+solvent: {1} \n detla Q: {2} | delta abs Q: {3} \n'.format(self.neutro.Energy,self.neutro.EnergyS,self.neutro.deltQ,self.neutro.deltAQ)
		glob_text += 'cation state info: \n'
		glob_text += '------\n'
		glob_text += 'Total energy: {0} | Total energy+solvent: {1} \n detla Q: {2} | delta abs Q: {3} \n'.format(self.cation.Energy,self.cation.EnergyS,self.cation.deltQ,self.cation.deltAQ)
		glob_text += 'anion state info: \n'
		glob_text += '------\n'
		glob_text += 'Total energy: {0} | Total energy+solvent: {1} \n detla Q: {2} | delta abs Q: {3} \n'.format(self.anion.Energy,self.anion.EnergyS,self.anion.deltQ,self.anion.deltAQ)
		glob_text += 'chem_pot hardness softness electrophilicity \n'
		glob_text += '{0} {1} {2} {3} \n '.format(self.chem_pot,self.hardness,self.softness,self.electrophilicity)
		glob_text += '---------------------------------------\n'
		glob_text += 'gloabal descriptors for the fragments \n'
		glob_text += 'frag_name chem_pot hardness softness electrophilicity \n'

		for i in range(len(self.frag_softness)):
			glob_text += '{0} {1} {2} {3} {4} \n'.format(self.neutro.Frag[i].name,self.frag_chempot[i],self.frag_hardness[i],self.frag_softness[i],self.frag_elect[i])

		glob_file.write(glob_text)
		glob_file.close()



#====================================================================

class local_rd:

	def __init__(self,
				neutro,
				cation, 
				anion):

		self.neutro = neutro
		self.cation = cation
		self.anion = anion
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
			self.fukuiES.append( self.anion.atoms[i].charge - self.neutro.atoms[i].charge )
			self.fukuiNS.append( self.neutro.atoms[i].charge - self.cation.atoms[i].charge )
			self.fukuiRS.append( (self.anion.atoms[i].charge - self.cation.atoms[i].charge)/2 )
			self.deltFukui.append( self.fukuiES[i] - self.fukuiNS[i] )
			self.softnesMu.append( self.globRD.softness*self.deltFukui[i])
			self.electMU.append( self.globRD.electrophilicity*self.deltFukui[i])
			'''
			if self.fukuiES[i] != 0.0:
				self.relativeNS.append ( self.fukuiNS[i]/self.fukuiES[i] )
			else:
				self.relativeES[i] = 10
			if self.fukuiNS[i] != 0.0: 
				self.relativeES.append ( self.fukuiES[i]/self.fukuiNS[i] )
			else:
				 self.relativeNS[i] = 10
			'''
	def writeLRD(self   , 
				 inpman):

		lrd_file = open(inpman,'w')
		lrd_text =''
		lrd_text += 'local reactivity descriptors by atom \n '
		lrd_text += 'atom fragment_num  fukuiES, fukuiNS,  fukuiRS, deltFukui, softnesMu, electMU \n'
		for i in range(len(self.neutro.atoms)):
			lrd_text += '{0} {1} {2:<05.5f} {3:<05.5f} {4:<05.5f} {5:<05.5f} {6:<05.5f} {7:<05.5f}  \n'.format(self.neutro.atoms[i].element,self.neutro.atoms[i].resNum,self.fukuiES[i],self.fukuiNS[i],self.fukuiRS[i],self.deltFukui[i],self.softnesMu[i],self.electMU[i]) 

		lrd_file.write(lrd_text)
		lrd_file.close()
		





