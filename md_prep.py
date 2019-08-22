#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

from pdb_class import *
from gmx_module import *
import parmed as pmd
import glob
import sys,os 


#Cofactor list
cofac_list = ["ATP","ADN","atp"]
#ATP atoms list
atp_list = ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]

#names of bin of amber and gromacs
path_amber = "/home/igorchem/programs/amber18/bin"
Reduce     = path_amber + "/reduce "
pdb4       = path_amber + "/pdb4amber "
antech     = path_amber + "/antechamber "
parmchk    = path_amber + "/parmchk2 "
tleap      = path_amber + "/tleap "

#=======================================================================
#***********************************************************************
def my_replace(fl,old,new):
				   	
	'''Function Doc
	Move to another lib file 
	'''
	with open(fl, 'r+') as f:
		s = f.read()
		s = s.replace(old, new)
		f.write(s)
			
#***********************************************************************
def fix_cofac_atoms(lig):	
	'''Function Doc	
	'''
	new_atoms   = []
	result_text = ""
	LIG         = ""
	
	cofac = protein(lig)
	cofac.pdb_parse()
	
	if lig[:3]  == "ATP":
		LIG  = "atp"
	pdb_res  = open(LIG+".pdb",'w')
	
	for atom in cofac.atoms:
		atom.resType = "atp"
		if atom.ptype[-1:] == "'":
			atom.ptype   = atom.ptype[:-1]+"*"		
		new_atoms.append(atom)
	
	i=0
	for atom in new_atoms:
		result_text+= "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,LIG,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
		i+=1
	pdb_res.write(result_text)
	pdb_res.close()
	return(LIG)
	
#***********************************************************************
def pdb_cat(pdb1,pdb2):		
	'''	Function Doc
	Move to pdb_class file
	'''
	
	pdba = protein(pdb1) 
	pdbb = protein(pdb2) 
	pdba.pdb_parse()
	pdbb.pdb_parse()	

	result_text = "HEADER complex of {} {}".format(pdb1,pdb2)
	
	pdb_res = open(pdb1[:-4]+"_comp.pdb","w")
	
	i = 0
	for atom in pdba.atoms:
		result_text+= "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
		i+=1
	for atom in pdbb.atoms:
		result_text+= "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
		i+=1	
		
	pdb_res.write(result_text)
	pdb_res.close
	

#=======================================================================
			
class md_prep:
	'''	Class Doc
	'''
	#-------------------------------------------
	
	def __init__(self,pdb):
		'''Method Doc
		'''
		self.pdb         = pdb
		self.current_pdb = pdb 
		self.lig         = "none"
		self.net_charge  = 0
		self.lig_charge  = 0 
		
	#-------------------------------------------
	
	def prepare_lig(self,lign,chg=0,rwat=True):	
		'''Method Doc
		'''
		
		lig              = []			
		self.lig         = lign
		pdb              = protein(self.pdb)
		self.lig_charge  = int(chg)
		
		pdb.pdb_parse()
		
		if rwat:
			pdb.remove_waters()
			
		for atom in pdb.atoms:
			if atom.resTyp == lign:
				lig.append(atom)
			
		text_lig = "HEADER LIG\n"
		lig_pdb = open(lign+".pdb",'w')
		i=0
		for atom in lig:
			text_lig += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
			i+=1
		lig_pdb.write(text_lig)
		lig_pdb.close()
		
		print("=======================================================")
		print("Removing the ligand atoms from the provided PDB")
		print("grep -v "+self.lig+" "+self.pdb+" > "+self.pdb[:-4]+"_w.pdb ")
		os.system("grep -v "+self.lig+" "+ self.pdb+" > "+self.pdb[:-4]+"_w.pdb ")
		os.system("sed 's/OXT/O  /' "+ self.pdb[:-4]+"_w.pdb > "+self.pdb[:-4]+"_wl.pdb ")
		
		print("=======================================================")
		print("Removing and adding hydrogens in the ligand pdb") 
		print(Reduce +"-Trim "+self.lig+".pdb > " + self.lig+"_h.pdb")
		os.system(Reduce + "-Trim "+self.lig+".pdb > " + self.lig+"_h.pdb")
		os.system(Reduce +self.lig+"_h.pdb > " + self.lig+".pdb")
		
		if self.lig in cofac_list:
			print("=======================================================")
			print("Ligand parameters will be loaded instead of created with ANTECHAMBER")
			self.lig = fix_cofac_atoms(self.lig+".pdb")
			print(self.lig)
			fl = os.listdir('.')
			if self.lig +".frcmod" in fl:
				print("FRCMOD OK...")
			else:				
				print("FRCMOD not found")
				sys.exit()
			if self.lig +".lib" in fl:		
				print("LIB OK...")
			else:
				print("LIB file not found")
				sys.exit()	
		else:
			print("=======================================================")
			print("Ligand parameters will be created with ANTECHAMBER.")		
			
			par = False
			fl = os.listdir('.')
			if self.lig+".frcmod" in fl:
				print("Found parameters for " + self.lig)
				print("FRCMOD file found for this ligand, antechamber parametrization will be skipped!") 
				par = True					
			if not par:
				print("===================================================")
				print("Run ANTECHAMBER:")
				print(antech+" -i "+self.lig+".pdb -fi pdb -o "+self.lig+".mol2 -fo mol2 -c bcc -nc "+chg)
				os.system(antech + " -i " + self.lig+".pdb -fi pdb -o " + self.lig+".mol2  -fo mol2 -c bcc -nc "+chg )
				os.system("rm ANTECHAMBER*")
		
				print("===================================================")
				print("Run Pamchek and generate frcmod")				
				print(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
				print(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
				os.system(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
				input()
				print("=======================================================")	
				print("Creating tleap input to save ligand library")
				tleap_in = "source leaprc.gaff2 \n"
				tleap_in += self.lig+" = loadmol2 "+self.lig+".mol2\n"
				tleap_in += "check "+self.lig+"\n"
				tleap_in += "loadamberparams " +self.lig+ ".frcmod\n"
				tleap_in += "saveoff " +self.lig+" "+self.lig+".lib \n"		
				tleap_in += "saveamberparm prot prmtop inpcrd\n"
				tleap_in += "quit"
		
				tleap_file = open('tleap_in','w')
				tleap_file.write(tleap_in)
				tleap_file.close()
				print("=======================================================")
				print("Run tleap and save the library with parameter ligands.")
				print(tleap + " -f tleap_in")
				os.system(tleap + " -f tleap_in" )
		
		self.current_pdb = self.pdb[:-4] +"_wl.pdb"

	def build_complex(self):
		
		print("=======================================================")
		print("Preparing Receptor/enzyme!")
		print(Reduce +" -Trim "+self.current_pdb+  " > " +self.current_pdb[:-4]+"_p.pdb") 
		os.system(Reduce + self.current_pdb+  " > " +self.current_pdb[:-4]+"_p.pdb")
		print("=======================================================")
		print(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		os.system(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		
		print("=======================================================")
		print("Concatenating Receptor/Enzyme with ligand/substrate")
		pdb_cat(self.current_pdb[:-4]+"_c.pdb",self.lig+".pdb")
		os.rename(self.current_pdb[:-4]+"_c_comp.pdb",self.pdb[:-4]+"_comp.pdb")
		self.current_pdb = self.pdb[:-4]+"_comp.pdb"
		
		print("=======================================================")	
		#Creating tleap input to save ligand library
		tleap_in =  "source oldff/leaprc.ff99SB \n"		
		tleap_in += "source leaprc.gaff2 \n"
		tleap_in += "loadamberparams " + self.lig + ".frcmod\n"
		tleap_in += "source leaprc.water.tip3p \n"
		tleap_in += "loadoff " + self.lig + ".lib\n"
		tleap_in += "complex = loadPdb " + self.current_pdb+ " \n"
		tleap_in += "solvatebox complex TIP3PBOX 12.0 \n"
		tleap_in += "addions2 complex Na+ 0\n"
		tleap_in += "addions2 complex Cl- 0\n"
		tleap_in += "savePdb complex "+self.current_pdb+"\n"
		tleap_in += "saveamberparm complex " +self.pdb[:-4]+".prmtop "+ self.pdb[:-4] +".inpcrd\n"
		tleap_in += "quit"

		tleap_file = open('tleap_in','w')
		tleap_file.write(tleap_in)
		tleap_file.close()
		
		print(tleap + " -f tleap_in")
		os.system(tleap + " -f tleap_in")
		
		#---------------------------------------------------------------
		
		fl = os.listdir('.')
		gromp = False
		print(self.pdb[:-4]+".top")
		if self.pdb[:-4]+".top" in fl:
			print("Found gromacs parameters for " + self.pdb)
			gromp = True
		elif not self.pdb[:-4]+".prmtop" in fl:
			print("prmtop not in folder for gromacs conversion")
			sys.exit()
		else:
			print("===================================================")
			print("Saving topologies for gromacs")
			ap = pmd.load_file(self.pdb[:-4]+".prmtop",self.pdb[:-4] +".inpcrd")
			ap.save(self.pdb[:-4]+".top")
			ap.save(self.pdb[:-4]+".gro")

	def min_gromacs(self):
		
		gromacs_inp()
		os.system("sed 's/WAT/SOL/' "+self.current_pdb+" > "+self.current_pdb[:-4]+"_t.pdb")
		os.rename(self.current_pdb[:-4]+"_t.pdb",self.current_pdb)
		os.system("sed 's/WAT/SOL/' "+self.pdb[:-4]+".top > "+self.pdb[:-4]+"_t.top")
		os.rename(self.pdb[:-4]+"_t.top",self.pdb[:-4]+".top")	
		
		print("=======================================================")
		print("Preparing the gromacs input files for structure minimization.")
		text_to_run = "/usr/bin/gmx" +" grompp -f em.mdp -c "+self.pdb[:-4]+ " -p "+ self.pdb[:-4] +".top -o em.tpr -maxwarn 50"
		os.system(text_to_run)
		
		print("=======================================================")
		print("Running minimization in gromacs.")
		text_to_run = "/usr/bin/gmx" + " mdrun -v -deffnm em"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Writting minimized structure pdb")
		print("gmx editconf -f em.gro -o "+self.pdb[:-4]+"_min.pdb")
		os.system("gmx editconf -f em.gro -o "+self.pdb[:-4]+"_min.pdb")	
		
	def equilibration(self):
		pass
	
	def production(self):
		pass 
		
	def organize_files():
		pass
		
	def process_traj():
		pass
		
		

