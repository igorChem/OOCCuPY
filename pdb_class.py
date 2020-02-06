#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pdb_class.py

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


ions = ["K+","Cl-","CL-","Na+","Mg+"]
#=======================================================================

# main classes coding

class pdb_atom:
	"""

	"""
	def __init__(self):
		self.name      = ''
		self.Type      = ''
		self.element   = ''
		self.ptype     = ''
		self.xcoord    = 0
		self.ycoord    = 0
		self.zcoord    = 0
		self.num       = 0
		self.resNum    = 0
		self.chain_t   = ''
		self.resTyp    = ''
		self.name      = ''
		self.charge    = 0
		self.occ       = 0
		self.bfactor   = 0

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

	def __init__(self,name,reorg=False):
		self.name           = name
		self.chain          = []
		self.resN           = 0
		self.atoms          = []
		self.total_charge   = 0
		self.waters         = []
		self.up_vertice     = [0,0,0]
		self.down_vertice   = [0,0,0]
		self.protein_center = [0,0,0]
		

		#---------------------------------------------------------------#
		# pdb parser 
		#---------------------------------------------------------------#
		
		i = 1
		pdb_file = open(self.name,'r')
		for line in pdb_file:
			if line[0:4]=="ATOM" or line[0:6]=="HETATM":
				a = pdb_atom()
				a.num     = i
				a.ptype   = line[12:16]
				a.resTyp  = line[17:20]
				a.chain_t = line[21:22]
				a.resNum  = int(line[22:26])
				a.xcoord  = float(line[30:38])
				a.ycoord  = float(line[38:46])
				a.zcoord  = float(line[46:54])
				a.occ     = float(line[56:60])
				a.bfactor = float(line[61:66])
				a.name    = a.Type + str(a.num)
				a.element = a.ptype[0:2]
				if a.element[0] =="1" or a.element[0]=="2" or a.element[0]=="3":
					a.element = "H"
				elif a.element == "He":
					a.element = "H"	
				elif a.ptype == "OXT":
					a.ptype = "O"					
				self.atoms.append(a)
				
			i+=1
		pdb_file.close()
		
		#---------------------------------------------------------------
		# residue definition 
		# --------------------------------------------------------------
		
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
		
		
		#--------------------------
		#define geometric properties of the protein
		#--------------------------
		xc =[]
		yc =[]
		zc =[]
		for i in range(len(self.atoms)):
			xc.append(self.atoms[i].xcoord) 
			yc.append(self.atoms[i].ycoord) 
			zc.append(self.atoms[i].zcoord)
		
		self.up_vertice[0]   = max(xc)
		self.up_vertice[1]   = max(yc)
		self.up_vertice[2]   = max(zc)
		self.down_vertice[0] = min(xc)
		self.down_vertice[1] = min(yc)
		self.down_vertice[2] = min(zc)
		
		self.protein_center[0] = (self.up_vertice[0] + self.down_vertice[0])/2
		self.protein_center[1] = (self.up_vertice[1] + self.down_vertice[1])/2
		self.protein_center[2] = (self.up_vertice[2] + self.down_vertice[2])/2
		
		
		#print properties  
		print("PDB file: "+self.name)
		print("Number of atoms: "+str(len(self.atoms)))
		print("Number of residues: "+str(len(self.chain)))
		print(self.protein_center)
		print(self.up_vertice)
		print(self.down_vertice)
		
	#-------------------------------------------------------------------	
	
	
	def remove_atom(self,i):
		del self.atoms[i]
	
	def remove_residue(self,i):
		for j in range(self.chain[i].atomsNum[0],self.chain[i].atomsNum[-1]):
			for k in range(len(self.atoms)):
				if self.atoms[k].num == j:
					del self.atoms[k]
		del self.chain[i]		
	

	def prune_pdb(self):
		a = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp[0]=="B":
				a.append(i)
				
		for i in sorted(a,reverse=True):
			del self.atoms[i]
	
	def prune_water(self,radius):
		for i in range(len(self.chain)):
			if self.chain[i].resTyp == "HOH" or  self.chain[i].resTyp == "WAT" or self.chain[i].resTyp == "SOL":
				pass
			

	def split_complex(self,lign):
		lig = []
		atoms_swap = []
		a   = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp==lign:
				lig.append(self.atoms[i])
				a.append(i)
				
		for i in sorted(a,reverse=True):
			del self.atoms[i]
				
		input_text ="HEADER {0} pdb file\n".format(self.name)

		i=1
		for atom in lig:
			input_text += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
			i+=1

		pdb = open(self.name[:-4]+"_lig.pdb",'w')
		pdb.write(input_text)
		pdb.close()

	def remove_waters(self):
		a = []
		for i in range(len(self.atoms)):
			if self.atoms[i].resTyp=="WAT" or self.atoms[i].resTyp=="HOH" or self.atoms[i].resTyp=="SOL":
				a.append(i)
				
		for i in sorted(a,reverse=True):
			del self.atoms[i]
				


		
		
		
	

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
		xyz=open(self.name[:-4] +'.xyz','w')
		
		for atom in self.atoms:
			if atom.element[0] == "H":				
				input_text += '{0} {1} {2} {3} \n'.format(atom.element[0],atom.xcoord,atom.ycoord,atom.zcoord)
			else:
				input_text += '{0} {1} {2} {3} \n'.format(atom.element,atom.xcoord,atom.ycoord,atom.zcoord)


		xyz.write(input_text)
		xyz.close()

	def write_pdb(self,filename):

		input_text ="HEADER {0} pdb file\n".format(self.name)

		i=1
		for atom in self.atoms:
			input_text += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
			i+=1

		pdb = open(filename,'w')
		pdb.write(input_text)
		pdb.close()

	def mopac_mop(self           ,
				  mode = "SP"    ,
				  mozyme = True  ,
				  solvent =True  ,
				  method ="PM7"  ):

		sol  = ""
		mozy = ""
		if solvent:
			sol = "EPS=78.4 RSOLV1.3"
		
		if mozyme:
			mozy = "mozyme"


		input_text = "{0} 1SCF large aux allvecs {1} {2}\n\n".format(method,mozy,sol,cutoff)
		chain = ""

		i=1
		for atom in self.atoms:
			input_text += "{0} {1}  1  {2} 1 {3} \n".format(atom.element,atom.xcoord,atom.ycoord,atom.zcoord)
			i+=1


		input_file =open(self.name[:-4] + ".mop",'w')
		input_file.write(input_text)
		input_file.close()
