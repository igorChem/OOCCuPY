#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mopac_out.py

#class to parse mopac output info

from pdb_class import*
#from xyz_class import*
#from cube_class import*

class mopac_out:

	def __init__(self,outfile, method="RM1"):
		self.name = outfile
		self.aux_name = outfile[:-4] + ".aux"
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
		self.AOpqn = []
		self.eigenv = []


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
		atom_pqn_in = 0
		atom_pqn_fin = 0
		charges_in   = 0
		charges_fin  = 0
		overlap_in   = 0
		overlap_fin  = 0
		denity_in   = 0
		density_fin  = 0
		eigenvalue_in = 0
		eigenvalue_fin = 0


		aux_file = open(self.aux_name,'r')
		indx = 0
		for line in aux_file:
			line2 = line.split()
			if line2 > 0:
				if line2[0][:8] == 'ATOM_EL[':
					if line2[0][len(line2[0])-2] == ']':
						self.numOfatoms = int(line2[0][9:len(line2[0])-2])
						atom_typ_in = indx
 				elif line2[0][:10] == 'ATOM_CORE[':
					atom_typ_fin = indx
				elif line2[0][:17] == 'ATOM_X:ANGSTROMS[':
					coords_in = indx
				elif line2[0][:13] == 'AO_ATOMINDEX[':
					coords_fin = indx
					atom_indx_in = indx
				elif line2[0][:13] == 'ATOM_SYMTYPE[':
					atom_indx_fin = indx
					atom_symtype_in = indx
				elif line2[0][:8] == 'AO_ZETA[':
					atom_symtype_fin = indx
					ao_zeta_in = indx
				elif line2[0][:9] == 'ATOM_PQN[':
					ao_zeta_fin = indx
					atom_pqn_in = indx
				elif line2[0][:14] == 'NUM_ELECTRONS=':
					atom_pqn_fin = indx
				elif line2[0][:13] == 'ATOM_CHARGES[':
					charges_in = indx
				elif line2[0][:15] == 'OVERLAP_MATRIX[':
					charges_fin = indx
					overlap_in = indx
				elif line2[0][:21] == 'TOTAL_DENSITY_MATRIX[' or line2[0][:21] == 'ALPHA_DENSITY_MATRIX[' or line2[0][:20] == 'BETA_DENSITY_MATRIX[':
					overlap_fin = indx
					density_in = indx
				elif line2[0][:25] == 'BETA_M.O.SYMMETRY_LABELS[' or line2[0][:20] == 'M.O.SYMMETRY_LABELS[':
					density_fin = indx
				elif line2[0][:36] == 'ALPHA_MOLECULAR_ORBITAL_OCCUPANCIES[' or line2[0][:12] == 'EIGENVALUES[':
					eigenvalue_in = indx
				elif line2[0][:30] ==  'MOLECULAR_ORBITAL_OCCUPANCIES[':
					eigenvalue_fin = indx
			indx +=1

		aux_file.close()

		l = 0
		m = 0
		with open(self.name,'r') as text:
			for (i,line) in enumerate(text):
				line2 = line.split()
				if i > atom_typ_in and i < atom_typ_fin:
					if len(line2) > 1:
						for k in range(len(line2)):
							atom = pdb_atom()
							atom.element = line2[k]
							self.atoms.append(atom)
				elif i > coords_in and i < coords_fin:
					if len(line2) == 3:
						self.atoms[l].xcoord = float(line2[0])
						self.atoms[l].ycoord = float(line2[1])
						self.atoms[l].zcoord = float(line2[2])
						l += 1
				elif i > atom_indx_in and i < atom_indx_fin:
					if len(line2) > 1:
						for k in range(len(line2)):
							orb = AOorbital()
							orb.aoindx = int(line2[k])
							self.atoms[int(line2[k])].orbs.append(orb)
				elif i > atom_symtype_in and i < atom_symtype_fin:
					if len(line2) > 1:
						for k in range(len(line2)):
							self.symmetries.append(line2[k])
				elif i > ao_zeta_in and i < ao_zeta_fin:
					if len(line2) > 1:
						for k in range(len(line2)):
							self.AOzetas.append(line2[k])
				elif i > atom_pqn_in and i < atom_pqn_fin:
					if len(line2) > 1:
						for k in range(len(line2)):
							self.AOpqn.append(line2[k])
				elif i > charges_in and i < charges_fin:
					if len(line2) > 1:
						for k in range(len(line2)):
							self.atoms[m].charge = line2[k]
							m += 1
				elif i > overlap_in and i < overlap_fin:
					for k in range(len(line2)):
						self.overlap.append(line2[k])
				elif i > density_in and i < density_fin:
					for k in range(len(line2)):
						self.m_dens.append(float(line2[k]))
				elif i > eigenvalue_in and i < eigenvalue_fin:
					for k in range(len(line2)):
						self.eigenv.append(float(line2[k]))

		zz = 0
		for x in range(len(self.atoms)):
			for y in  range(len(self.atoms[x].orbs)):
				self.atoms[x].orbs[y].pqn = self.AOpqn[zz]
				self.atoms[x].orbs[y].symmetry = self.symmetries[zz]
				self.atoms[x].orbs[y].zeta = self.AOzetas[zz]
				zz +=1

	def write_report(self):

		report_file = open(self.name+".rep",'w')

		report_text = ""
		report_text += "{0} \n".format(self.name)
		report_text += "{0} \n".format(self.numOfatoms)
		report_text += "{0} \n".format(self.energy)

		for i in range(len(self.atoms)):
			if len(self.atoms[i].orbs) == 1:
				report_text += "{0} {1} {2} {3} {4} \n".format(self.atoms[i].element,self.atoms[i].xcoord,self.atoms[i].ycoord,self.atoms[i].zcoord,self.atoms[i].charge)
			elif len(self.atoms[i].orbs) == 4:
				report_text += "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}\n".format(atom.element,
															  atom.xcoord,
															  atom.ycoord,
															  atom.zcoord,
															  atom.charge,
															  atom.orbs[0].symmetry,
															  atom.orbs[0].pqn,
															  atom.orbs[0].zeta,
															  atom.orbs[1].symmetry,
															  atom.orbs[1].pqn,
															  atom.orbs[1].zeta,
															  atom.orbs[2].symmetry,
															  atom.orbs[2].pqn,
															  atom.orbs[2].zeta,
															  atom.orbs[3].symmetry,
															  atom.orbs[3].pqn,
															  atom.orbs[3].zeta)

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
