#!/usr/bin/env python
# -*- coding: utf-8 -*-
# primordia_module.py

import os, glob

class primordia:

	def __init__(self,ctrs,grid=0,charge=1):
		self.neutro_list   = []
		self.cation_list   = []
		self.anion_list    = []
		self.program       = []
		self.homo_cube     = []
		self.lumo_cube     = []
		self.elecdens_cube = []
		self.cation_cube   = []
		self.anion_cube    = []
		self.charge        = charge
		self.grid          = grid
		self.ctrs          = ctrs
		
		if ctrs == 1 or ctrs == 3:
			self.neutro_list = glob.glob("*neutro.out")
		elif ctrs == 2:
			self.neutro_list = glob.glob("*neutro.out")
			self.cation_list = glob.glob("*cation.out")
			self.anion_list  = glob.glob("*anion*.out")		
		elif ctrs == 4:
			self.neutro_list = glob.glob("neutro.aux")			
			self.homo_cube = glob.glob("*HOMO.cube")
			self.lumo_cube = glob.glob("*LUMO.cube")
		elif ctrs == 5 or ctrs == 6:
			self.neutro_list = glob.glob("*neutro.aux")
			self.cation_list = glob.glob("*cation.aux")
			self.anion_list  = glob.glob("*anion.aux")	
		elif ctrs == 7:
			self.neutro_list = glob.glob("*neutro,aux")
			self.cation_list = glob.glob("*cation.aux")
			self.anion_list  = glob.glob("*anion.aux")
			self.elecdens    = glob.glob("*neutro.cube")
			self.cation_cube = glob.glob("*cation.cube")
			self.anion_cube  = glob.glob("*anion.cube")			
				
	def executable_primorida(self):
		pass
				
	def executable_compPrimordia(self):
		pass		
		
	def input_autoprimorida(self):		
		
		filled = open("auto_inp",'w')
		fil_text = "{0} {1} {2} {3} {4} \n".format(len(self.neutro_list),self.grid,self.charge,self.ctrs)
		
		if ctrs == 1: 
			for i in range(len(self.neutro_list)):
				fil_text += "{0} mopac \n".format(self.neutro_list[i])

		
		filled.write(fil_text)
		filled.close()
		
		
		
		
		
		
