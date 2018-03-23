#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mopac_out.py

#class to parse mopac output info 

from pdb_class import*
from xyz_class import*
from cube_class import*


class mopac_out:

	def __init__(self,outfile, method="RM1"):
		self.name = outfile
		self.aux_name = outfile[:-4] + "aux"
		self.energy = 0
		self.heat = 0
		self.energy = 0
		self.atoms = []
		self.homo_en = 0
		self.lumo_en = 0
		self.method = method
		self.numOfatoms = 0

	def parse_out(self):

		out_file = open(self.name,'r')

		for line in out_file:
			line2 = line.split()
			if len(line2) == 9:
				if line2[1] == "HEAT" and line2[3] == "FORMATION":
					self.heat = line2[7]
			elif len(line2) == 5:
				if line2[0] == "TOTAL" and line2[1] == "ENERGY":
					self.energy = line2[3]
			elif len(line2) == 7:
				if line2[0] == "HOMO" and line2[1] == "LUMO":
					self.homo_en = line2[5]
					self.lumo_en = line2[6]


		phrase1 = 'ATOM NO.   TYPE          CHARGE      No. of ELECS.   s-Pop       p-Pop'
		phrase2 = 'DIPOLE           X         Y         Z       TOTAL'
		
		
		ch_init = 0
		ch_fin  = 0
				
		with open(self.name,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					ch_init=i								
				elif phrase2 in line:
					ch_fin=i
					
		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				if i > ch_init and i <ch_fin :
					line2 = line.split()
					if len(line2) > 3:
						atom = pdb_atom()
						atom.num = line2[0]
						atom.element = line2[1]
						atom.charge = line2[2]
						print(line2[3],i)
						self.atoms.append(atom)
						self.numOfatoms +=1


	def write_report(self):

		report_file = open(self.name+".rep",'w')

		report_text = ""
		report_text += "qunatum report of {0} \n".format(self.name)
		report_text += "number of atoms in molecule {0}\n".format(self.numOfatoms)
		report_text += "Total energy = {0} | Heat of formation = {1}\n".format(self.energy,self.heat)
		report_text += "atoms and its muliken charges \n "

		for atom in self.atoms:
			report_text += "{0} {1} {2} \n".format(atom.num,atom.element,atom.charge)

		report_file.write(report_text)
		report_file.close()

a = mopac_out("1l2yRM1neutro.out")
a.parse_out()
a.write_report()


