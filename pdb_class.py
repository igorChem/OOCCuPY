#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pdb_Reader.py

#=======================================================================

#load modules

import os

#=======================================================================

# residues data dictionaries

RES_N0C ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}
RES_N0H ={'GLY':3,'ALA':5,'VAL':9,'PHE':9,'ILE':11,'LEU':11,'PRO':7,'MET':9,'ASP':4,'GLU':6,'LYS':13,'ARG':13,'SER':5,'THR':7,'TYR':9,'CYS':5,'ASN':6,'GLN':8,'HIS':8,'TRP':10}
RES_N0N ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}
RES_N0O ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}
RES_N0S ={'GLY':2,'ALA':3,'VAL':5,'PHE':9,'ILE':6,'LEU':6,'PRO':5,'MET':5,'ASP':4,'GLU':5,'LYS':6,'ARG':6,'SER':3,'THR':4,'TYR':9,'CYS':3,'ASN':4,'GLN':5,'HIS':6,'TRP':11}

#=======================================================================

# main classes coding

class AOorbital:

    def __init__(self):
        self.symmetry = ''
        self.pqn      = 0
        self.zeta     = 0.0
        self.aoindx   = 0


class pdb_atom:

	def __init__(self):
		self.name      = ''
		self.Type      = ''
		self.element   = ''
		self.ptype     = ''
		self.xcoord    = []
		self.ycoord    = []
		self.zcoord    = []
		self.num       = []
		self.resNum    = []
		self.resTyp    = ''
		self.name      = ''
		self.charge    = 0
        self.orbs      = []

class residue:

	def __init__(self):

		obj = pdb_atom

		self.name        = ''
		self.typ         = ''
		self.num         = ''
		self.alfaC       = obj
		self.Nitrogen    = obj
		self.carb		 = obj
		self.oxygen      = obj
		self.carbB       = obj
		self.r_atoms     = []
		self.atomsNum    = []
		self.hydrogen    = 0
		self.charge      = 0
		self.side_chain  = []

	def reorg(self):

		self.r_atoms.append(self.Nitrogen)
		self.r_atoms.append(self.alfaC)
		self.r_atoms.append(self.carb)
		self.r_atoms.append(self.oxygen)
		self.r_atoms.append(self.carbB)

		for i in range(len(self.side_chain)):
			self.r_atoms.append(self.side_chain[i])
			#print(self.r_atoms[i].ptype)



class protein:

	def __init__(self         ,
				 name         ,
				 amber = False):
		self.name  = name
		self.chain = []
		self.resN  = 0
		self.atoms = []
		self.amber = amber

	def pdb_parse(self      ,
				 filename   ):

		pdb_file = open(filename,'r')

		for line in pdb_file:
			line = line.split()
			if self.amber == False:
				if line[0]=='ATOM':
					if not line[3] == "WAT":
						a = pdb_atom()
						a.num = line[1]
						a.ptype = line[2]
						a.resTyp = line[3]
						a.resNum = line[5]
						a.xcoord = float(line[6])
						a.ycoord = float(line[7])
						a.zcoord = float(line[8])
						a.Type = line[11]
						a.element = a.Type[0][0]
						a.name = a.Type + str(a.num)
						self.atoms.append(a)
					else:
						print("water atom")
			elif self.amber == True:
				if line[0]=='ATOM':
					if not line[3] == "WAT":
						a = pdb_atom()
						a.num = line[1]
						a.ptype = line[2]
						a.resTyp = line[3]
						a.resNum = line[4]
						a.xcoord = round(float(line[5]),3)
						a.ycoord = round(float(line[6]),3)
						a.zcoord = round(float(line[7]),3)
						a.Type = line[10]
						a.element = a.Type[0][0]
						a.name = a.Type + str(a.num)
						self.atoms.append(a)


		return(self.atoms)

		pdb_file.close()


	def residue_def(self,reorg = False):

		for i in range(len(self.atoms)):
			if self.atoms[i].ptype == 'N':
				a=residue()
				a.name = self.atoms[i].resTyp + self.atoms[i].resNum
				a.typ = self.atoms[i].resTyp
				a.num = self.atoms[i].resNum
				a.atomsNum.append(self.atoms[i].num)
				self.chain.append(a)

		for i in range(len(self.chain)):
			for j in range(len(self.atoms)):
				if self.chain[i].num == self.atoms[j].resNum:
					self.chain[i].atomsNum.append(self.atoms[j].num)
					if self.atoms[j].ptype == 'N':
						self.chain[i].Nitrogen = self.atoms[j]
					elif self.atoms[j].ptype == 'CA':
						self.chain[i].alfaC = self.atoms[j]
					elif self.atoms[j].ptype =='C':
						self.chain[i].carb = self.atoms[j]
					elif not self.chain[i].typ =='GLY' and self.atoms[j].ptype == 'CB':
						self.chain[i].carbB = self.atoms[j]
					elif self.chain[i].typ =='GLY' and self.atoms[j].ptype == 'H':
						self.chain[i].carbB = self.atoms[j]
					elif self.atoms[j].ptype == 'O' or self.atoms[j].ptype == 'OC1':
						self.chain[i].oxygen = self.atoms[j]
					elif self.atoms[j].element =='H':
						self.chain[i].hydrogen +=1
						self.chain[i].side_chain.append(self.atoms[j])
					else:
						self.chain[i].side_chain.append(self.atoms[j])


		if reorg == True:
			for i in range(len(self.chain)):
				self.chain[i].reorg()

		for i in range(len(self.chain)):
			for j in range(len(self.chain[i].r_atoms)):
				print(i,j,self.chain[i].typ,i+j,self.chain[i].r_atoms[j].ptype)

			cnt = 0
			for i in range(len(self.chain)):
				for j in range(len(self.chain[i].r_atoms)):
					self.atoms[cnt] = self.chain[i].r_atoms[j]
					cnt +=1


		self.resN = len(self.chain)

	def charge_res(self):


		for i in range(len(self.chain)):

			# for asp and glu
			if self.chain[i].typ == 'ASP' or self.chain[i].typ == 'GLU':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ].typ +1:
						self.chain[i].charge = -1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ].typ +2:
						self.charge = 0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ].typ +3:
						self.chain[i].charge = 1
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = -2
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +1:
						self.chain[i].charge = -1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +2:
						self.chain[i].charge = 0
					else:
						continue
				else:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = -1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge = 0
					else:
						continue
			#for lys and arg
			elif self.chain[i].typ == 'LYS' or self.chain[i].typ == 'ARG':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge = 1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+2:
						self.chain[i].charge = 2
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = 0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +1:
						self.chain[i].charge = 1
					else:
						continue
				else:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = 1
					else:
						continue
			#for his
			elif self.chain[i].typ == 'HIS':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge =1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] +2:
						self.chain[i].charge = 2
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]-1:
						self.chain[i].charge=-1
					else:
						continue
				else:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge = 1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]-1:
						self.chain[i].charge=0
					else:
						continue
			#for cys
			elif self.chain[i].typ == 'CYS':
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen   == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge = 0
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+2:
						self.chain[i].charge = 1
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=-1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=0
					else:
						continue
				else:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ] or self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]-1:
						self.chain[i].charge=0
					else:
						continue
			#for any residue
			else:
				if self.chain[i].num == self.chain[0].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+2:
						self.chain[i].charge=1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=0
					else:
						continue
				elif self.chain[i].num == self.chain[-1].num:
					if self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]:
						self.chain[i].charge=-1
					elif self.chain[i].hydrogen == RES_N0H[self.chain[i].typ]+1:
						self.chain[i].charge=0
					else:
						continue
				else:
					self.chain[i].charge = 0


	def write_xyz(self):

		input_text = '{0} \n \n'.format(len(self.atoms))

		for atom in self.atoms:
			input_text += '{0} {1} {2} {3} \n'.format(atom.element,atom.xcoord,atom.ycoord,atom.zcoord)

		xyz =open(self.name +'.xyz','w')
		xyz.write(input_text)
		xyz.close()

		return(input_text)

	def write_pdb(self,filename):

		input_text =''

		if self.amber == True:
			chain = ''
		else:
			chain = 'A'

		i=1
		for atom in self.atoms:
			input_text += "ATOM {0:>4} {1:<5} {2:<4} {3} {4:<4} {5:<05.3f} {6:<05.3f} {7:<05.3f}  1.00  0.00  {8:>10} \n".format(i,atom.ptype,atom.resTyp,chain,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.element)
			i+=1

		pdb =open(filename,'w')
		pdb.write(input_text)
		pdb.close()

	def mopac_mop(self           ,
				  mode = "SP"    ,
				  max_time = "1D",
				  solvent =True  ,
				  cutoff = 9.0   ,
				  method ="PM7"  ):


		if solvent == True:
			sol = "EPS=78.4 RSOLV1.3"
		else:
			sol = ""

		input_text = "{0} 1SCF XYZ PDB PL T={1} TIMES MOZYME {2} CUTOFF={3} \n\n".format(method,max_time,sol,cutoff)

		chain = ""

		i=1
		for atom in self.atoms:
			input_text += "ATOM {0:>6} {1:<5} {2:<4} {3} {4:<4} {5:<05.3f} {6:<05.3f} {7:<05.3f}  1.00  0.00  {8:>10} \n".format(i,atom.ptype,atom.resTyp,chain,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.element)
			i+=1


		input_file =open(self.name[:-4] + ".mop",'w')
		input_file.write(input_text)
		input_file.close()
