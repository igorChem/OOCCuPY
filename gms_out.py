#!/usr/bin/env python
# -*- coding: utf-8 -*-
# gms_out.py


from pdb_class import*
from xyzclass import*
from cube_class import*


class gms_out:

	def __init__(self,logfile):
		self.name = logfile;
		self.dat_file = logfile[:-4]
		self.energy = 0
		self.enegyPCM = 0
		self.total_charge = 0
		self.atoms = []
		self.numOfatoms = 0
		self.homo_en = 0
		self.lumo_en = 0 


	def parse_log(self):

		logFile = open(logfile,'r')

		for line in logfile:
			line2 = line.split()
			if len(line2) == 5:
				if line2[0] == "CHARGE" and line2[2] == "MOLECULE":
					self.total_charge = line2[4]
			elif len(line2) == 11:
				if line2[0] == "INTERNAL" and line2[3] == "SOLVENT":
					self.enegyPCM=line2[9]
			elif len(line2) == 8:
				if line2[0] == "FINAL" and line2[2] == "ENERGY":
					self.energy =  line2[4] 


		phrase1 = 'ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE'
		phrase2 = 'BOND ORDER AND VALENCE ANALYSIS     BOND ORDER THRESHOLD=0.050'
		
		
		ch_init = 0
		ch_fin  = 0
				
		with open(self.datfile,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					ch_init=i								
				elif phrase2 in line:
					ch_fin=i
					
		with open(self.datfile,'r') as text:
			for (i,line) in enumerate(text):
				if i >= ch_init and i <=ch_fin :
					line2 = line.split()
					if len(line2) > 3:
						atom = pdb_atom()
						atom.num = line2[0]
						atom.element = line2[1]
						atom.charge = line2[3]
						self.atoms.append(atom)
						self.numOfatoms +=1

	def parse_dat(self):

		#separate the cube information and writes to a file
		phrase1 = 'GAMESS CUBE FORMAT: ELECTRON DENSITY'
		phrase2 = '$END'		
					
		ir_init = 0
		ir_fin  = 0
		cube_info =''
		
		with open(self.datfile,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					ir_init=i								
				elif phrase2 in line:
					ir_fin=i
					
		with open(self.datfile,'r') as text:
			for (i,line) in enumerate(text):
				if i >= ir_init and i <=ir_fin :
					line2 = line.split()
					if len(line2) > 3:
						cube_info += line
					 
		cube_file0 = open(self.name[:-4]+'.cube','w')
		cube_file0.write(cube_info)
		cube_file0.close()

		

	def write_report(self,file_name):
		

		report_file = open(name+".rep",'r')

		report_text = ""
		report_text += "qunatum report of {0} \n".format(self.name)
		report_text += "number of atoms in molecule {0} and total charge {1}\n".format(self.numOfatoms,self.total_charge)
		report_text += "Total energy = {0} | Total energy in solvent = {1}\n".format(self.energy,self.enegyPCM)
		report_text += "atoms and its muliken charges \n "

		for atom in self.atoms:
			report_text += "{0} {1} {2} \n".format(atom.num,atom.element,atom.charge)

		report_file.write(report_text)
		report_file.close()




