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
		self.gap = 0
		self.method = method
		self.numOfatoms = 0
		self.overlap = []
		self.m_dens = []
		self.MO = []
		self.dipole = []
		self.Tdipole = [] 
		self.symmetries = []
		self.AOzetas = []
		self.AOatomIndices = []
		self.AOpqn = []


	def parse_out(self):

		out_file = open(self.name,'r')

		for line in out_file:
			line2 = line.split()
			if len(line2) == 10:
				if line2[1] == "HEAT" and line2[3] == "FORMATION":
					#print(line)
					self.heat = float(line2[5])
			elif len(line2) == 5:
				if line2[0] == "TOTAL" and line2[1] == "ENERGY":
					self.energy = float(line2[3])
			elif len(line2) == 7:
				if line2[0] == "HOMO" and line2[1] == "LUMO":
					self.homo_en = float(line2[5])
					self.lumo_en = float(line2[6])
					self.gap = self.lumo_en - self.homo_en


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
						#print(line2[3],i)
						self.atoms.append(atom)
						self.numOfatoms +=1

		out_file.close();

	def parse_aux(self):


		atom_typ_in   = 0
		atom_typ_fin  = 0
		atom_indx_in  = 0
		atom_indx_fin = 0
		coords_in     = 0
		coords_fin    = 0
		atom_symtype_in = 0
		atom_symtype_fin = 0
		ao_zeta_in  = 0
		ao_zeta_fin = 0

		out_file = open(self.name,'r')
		indx = 0
		for line in out_file: 
			line2 = line.split()
			if line2 > 0:
				if line2[0][:8] == 'ATOM_EL[': 
					if line2[0][12] == ']':
						self.numOfatoms = int(line2[0][9:12])
						atom_typ_in = indx
					else:
						self.numOfatoms = int(line2[0][9:13])
						atom_typ_in = indx
				elif line2[0][:10] == 'ATOM_CORE[':
					atom_typ_fin = i
				elif line2[0][:17] == 'ATOM_X:ANGSTROMS[':
					coords_in = i
				elif line2[0][:12] == 'AO_ATOMINDEX[':
					coords_fin = i
					atom_indx_in = i
				elif line2[0][:13] == 'ATOM_SYMTYPE[':
					atom_indx_fin = i
					atom_symtype_in = i
				elif line2[0][:8] == 'AO_ZETA[':
					atom_symtype_fin = i
					ao_zeta_in = i
				elif line2[0][:9] == 'ATOM_PQN[':
					pass

			indx +=1


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


def all_out():

	list = glob.glob('*.out')
	text = 'name heat_of_formation gap \n'

	for out in list:
		obj = mopac_out(out)
		obj.parse_out()
		text+= "{0} {1} {2} \n".format(obj.name,obj.heat,obj.gap)  

	filerep = open("reportmopac",'w')
	filerep.write(text)
	filerep.close()



'''
a = mopac_out("1l2yRM1neutro.out")
a.parse_out()
a.write_report()
'''

