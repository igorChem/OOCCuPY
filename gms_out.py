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
		self.energy = 0;
		self.enegyPCM = 0; 
		self.total_charge = 0; 
		self.atoms = []
		self.numOfatoms = 0;

	def parse_log(self):

		logFile = open(logfile,'r')

		for line in logfile:
			line2 = line.split()
			if len(line2) == 5:
				self.total_charge = 

		

	def parse_dat(self):

		#separate the cube information and writes to a file
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

		

	def write_report(self,file_name):
		pass



