#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

from pdb_class import *
from gmx_module import *
import parmed as pmd
import glob
import sys,os 


#Cofactor list
cofac_list = ["ATP","atp","NADP","NADH","NAD","ADP"]
#ATP atoms list
atp_list = ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]

#names of bin of amber and gromacs
path_amber = "/home/igorchem/programs/amber18/bin"
Reduce     = path_amber + "/reduce "
pdb4       = path_amber + "/pdb4amber "
antech     = path_amber + "/antechamber "
parmchk    = path_amber + "/parmchk2 "
tleap      = path_amber + "/tleap "
pymol      = "/home/igorchem/programs/bin/pymol "
path_cofac = "/home/igorchem/OOCCuPY/cofac/"
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
	LIG         = lig
	
	cofac = protein(lig)
	cofac.pdb_parse()
	
	
	if lig[:3]  == "ATP":
		LIG  = "atp"
		os.system("cp " + path_cofac +"atp* .")
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
def pdb_cat(pdb1,pdb2,flname):		
	'''	Function Doc
	Move to pdb_class file
	'''
	
	pdba = protein(pdb1) 
	pdbb = protein(pdb2) 
	pdba.pdb_parse()
	pdbb.pdb_parse()	

	result_text = "HEADER complex of {} {}".format(pdb1,pdb2)
	
	pdb_res = open(flname,"w")
	
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
		self.lig         = []
		self.net_charge  = 0
		self.lig_charge  = []
		self.num_lig     = 0 
		
	#-------------------------------------------
	
	def prepare_lig(self,nlig,lign,chg,mult,rwat=True,lig_hy="False"):
		'''Method Doc
		'''
		
		lig_h = False
		if lig_hy == "T":
			lig_h = True
		pdb              = protein(self.pdb)		
		self.num_lig     = int(nlig)
		self.current_pdb = self.pdb
		
		print("Paramters parsed:\n")
		print("pdb file: " + self.pdb)
		print("Num of ligands: "+"self.num_lig")
		
		for i in range(len(lign)):
			print("Lig #" + str(i)+": " + lign[i]) 
			print("Lig charge #" + str(i)+": " + str(chg[i])) 
			print("Lig multiplicity #" + str(i)+": " + str(mult[i])) 
			
		pdb.pdb_parse()
		if rwat:
			pdb.remove_waters()
			
		for j in range(self.num_lig):
			lig = []
			self.lig_charge.append(int(chg[j]))
			self.lig.append(lign[j])
					
			for atom in pdb.atoms:
				if atom.resTyp == lign[j]:
					lig.append(atom)
			
			text_lig = "HEADER LIG\n"
			lig_pdb = open(lign[j]+".pdb",'w')
			i=0
			
			for atom in lig:
				text_lig += "ATOM {0:6} {1:4} {2:2} {3:<1} {4:<7} {5:7.3f} {6:7.3f} {7:7.3f} {8:>5.2f} {9:>4.2f} \n".format(i,atom.ptype,atom.resTyp,atom.chain_t,atom.resNum,atom.xcoord,atom.ycoord,atom.zcoord,atom.occ,atom.bfactor)
				i+=1			
			lig_pdb.write(text_lig)
			lig_pdb.close()
			
			print("=======================================================")
			print("Removing the ligand atoms from the provided PDB")
			if self.num_lig == 1:
				print("grep -v "+self.lig[j]+" "+self.pdb+" > "+self.pdb[:-4]+"_w.pdb ")
				os.system("grep -v "+self.lig[j]+" "+ self.pdb+" > "+self.pdb[:-4]+"_w.pdb ")
				self.current_pdb  = self.pdb[:-4]+"_w.pdb"
			elif self.num_lig > 1:
				print("grep -v "+self.lig[j]+" "+self.current_pdb+" > "+self.pdb[:-4]+"_w.pdb ")
				os.system("grep -v "+self.lig[j]+" "+ self.current_pdb+" > "+self.pdb[:-4]+"_"+str(j)+"_.pdb")
				self.current_pdb = self.pdb[:-4]+"_"+str(j)+"_.pdb"
			
			print("=======================================================")
			if lig_h:
				print("Removing and adding hydrogens in the ligand pdb")
				print(Reduce +"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")
				os.system(Reduce+"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")
				if not self.lig[j] in cofac_list:
					print(Reduce +self.lig[j]+"_h.pdb > " + self.lig[j]+".pdb")
					os.system(Reduce +self.lig[j]+"_h.pdb > " + self.lig[j]+".pdb")
				else:
					if self.lig[j] in cofac_list:
						os.system("cp " +self.lig[j]+"_h.pdb "+ self.lig[j]+".pdb")
			else:
				if self.lig[j] in cofac_list:
					print("Removing and adding hydrogens in the ligand pdb")
					print(Reduce +"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")
					os.system(Reduce+"-Trim "+self.lig[j]+".pdb > " + self.lig[j]+"_h.pdb")	
					os.system("cp " +self.lig[j]+"_h.pdb "+ self.lig[j]+".pdb")
					
			os.system("cat < "+ self.lig[j]+".pdb")	
			
			if self.lig[j] in cofac_list:
				os.system( "cp " + path_cofac + "*lib "+	os.getcwd() )		
				os.system( "cp " + path_cofac + "*frcmod "+	os.getcwd() )
				os.system( "cp " + path_cofac + "*prep   "+	os.getcwd() )
				
				print("=======================================================")
				print("Ligand parameters will be loaded instead of created with ANTECHAMBER")
				self.lig[j] = fix_cofac_atoms(lign[j]+".pdb")				
				print(self.lig[j])
				fl = os.listdir('.')
				print(fl)
				if self.lig[j][:-4] +".frcmod" in fl:
					print("FRCMOD OK...")
				else:				
					print("FRCMOD not found")				
					sys.exit()	
				if self.lig[j] +".lib" in fl:		
					print("LIB OK...")
				else:
					print("=======================================================")	
					print("Creating tleap input to save ligand library")
					tleap_in =  "source leaprc.gaff2\n"
					tleap_in += "loadamberparams " +self.lig[j][:-4]+ ".frcmod\n"
					tleap_in +=  self.lig[j][:-4]+" = loadPdb "+self.lig[j]+"\n"
					tleap_in += "check "+self.lig[j][:-4]+"\n"					
					tleap_in += "saveoff " +self.lig[j][:-4]+" "+self.lig[j][:-4]+".lib \n"		
					tleap_in += "quit"
					tleap_file = open("tleap_in"+"_"+self.lig[j][:-4],'w')
					tleap_file.write(tleap_in)
					tleap_file.close()
					print("=======================================================")
					print("Run tleap and save the library with parameter ligands.")
					print(tleap + " -f tleap_in")
					os.system(tleap + " -f tleap_in"+"_"+self.lig[j][:-4])
				
			else:
				print("=======================================================")
				print("Ligand parameters will be created with ANTECHAMBER.")		
			
				par = False
				fl = os.listdir('.')
				if self.lig[j]+".frcmod" in fl:
					print("Found parameters for " + self.lig[j])
					print("FRCMOD file found for this ligand, antechamber parametrization will be skipped!") 
					self.lig[j] = self.lig[j]+".pdb"
					par = True
										
				if not par:
					print("===================================================")
					print("Run ANTECHAMBER:")
					print(antech+" -i "+self.lig[j]+".pdb -fi pdb -o "+self.lig[j]+".mol2 -fo mol2 -c bcc -nc "+chg[j]+" -m "+mult[j])
					os.system(antech + " -i " + self.lig[j]+".pdb -fi pdb -o " + self.lig[j]+".mol2  -fo mol2 -c bcc -nc "+chg[j]+" -m "+mult[j] )
					os.system("rm ANTECHAMBER*")
					print("===================================================")
					print("Run Pamchek and generate frcmod")				
					print(parmchk+" -i "+self.lig[j]+".mol2 -f mol2 -o " +self.lig[j]+".frcmod")
					print(parmchk+" -i "+self.lig[j]+".mol2 -f mol2 -o " +self.lig[j]+".frcmod")
					os.system(parmchk + " -i "+ self.lig[j]+".mol2 -f mol2 -o " + self.lig[j]+".frcmod")
										
					print("=======================================================")	
					print("Creating tleap input to save ligand library")
					tleap_in = "source leaprc.gaff2 \n"
					tleap_in += "loadamberparams " +self.lig[j]+ ".frcmod\n"
					tleap_in +=  self.lig[j]+" = loadmol2 "+self.lig[j]+".mol2\n"
					tleap_in += "check "+self.lig[j]+"\n"					
					tleap_in += "saveoff " +self.lig[j]+" "+self.lig[j]+".lib \n"		
					tleap_in += "quit"
		
					tleap_file = open("tleap_in"+"_"+self.lig[j],'w')
					tleap_file.write(tleap_in)
					tleap_file.close()
					print("=======================================================")
					print("Run tleap and save the library with parameter ligands.")
					print(tleap + " -f tleap_in")
					os.system(tleap + " -f tleap_in"+"_"+self.lig[j])
					self.lig[j] = self.lig[j]+".pdb"
				
		os.system("sed 's/OXT/O  /' "+self.current_pdb+" > "+self.pdb[:-4]+"_wl.pdb ")
		self.current_pdb = self.pdb[:-4] +"_wl.pdb"

	def build_complex(self,addH=True):
		
		print("=======================================================")
		print("Preparing Receptor/enzyme!")
		print(Reduce +" -Trim "+self.current_pdb+  " > " +self.current_pdb[:-4]+"_p.pdb") 
		os.system(Reduce+" -Trim "+self.current_pdb+" > "+self.current_pdb[:-4]+"_p.pdb")
		print("=======================================================")
		print(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		os.system(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		#os.system(pdb4 +  self.current_pdb+" > " +self.current_pdb[:-4]+"_c.pdb")
		
		self.current_pdb = self.current_pdb[:-4]+"_c.pdb"
		
		print("=======================================================")
		print("Concatenating Receptor/Enzyme with ligand/substrate")
		for i in range(self.num_lig):
			print(self.lig[i])
			pdb_cat(self.current_pdb,self.lig[i],self.current_pdb)	
						
		os.rename(self.current_pdb,self.pdb[:-4]+"_comp.pdb")
		self.current_pdb = self.pdb[:-4]+"_comp.pdb"

		print("=======================================================")	
		#Creating tleap input to save ligand library
		tleap_in =  "source oldff/leaprc.ff99SB \n"		
		tleap_in += "source leaprc.gaff2 \n"		
		tleap_in += "source leaprc.water.tip3p \n"
		for i in range(self.num_lig):
			tleap_in += "loadamberparams " + self.lig[i][:-4] + ".frcmod\n"
			tleap_in += "loadoff " + self.lig[i][:-4] + ".lib\n"
		tleap_in += "complex = loadPdb " + self.current_pdb+ " \n"
		tleap_in += "solvatebox complex TIP3PBOX 12.0 \n"
		tleap_in += "check complex \n"
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
		os.system("sed 's/WAT/SOL/' "+self.pdb[:-4]+".gro > "+self.pdb[:-4]+"_t.gro")
		os.rename(self.pdb[:-4]+"_t.gro",self.pdb[:-4]+".gro")
		
		print("=======================================================")
		print("Preparing the gromacs input files for structure minimization.")
		text_to_run = "/usr/bin/gmx_d" +" grompp -f em.mdp -c "+self.pdb[:-4]+ " -p "+ self.pdb[:-4] +".top -o em.tpr -maxwarn 50"
		os.system(text_to_run)
		
		print("=======================================================")
		print("Running minimization in gromacs.")
		text_to_run = "/usr/bin/gmx_d" + " mdrun -v -deffnm em"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Writting minimized structure pdb")
		print("/usr/bin/gmx_d editconf -f em.gro -o "+self.pdb[:-4]+"_min.pdb")
		os.system("gmx editconf -f em.gro -o "+self.pdb[:-4]+"_min.pdb")	
		
		print("=======================================================")
		print("Preparing gromacs NVT equilibration")
		print("gmx grompp -f nvt.mdp -c em.gro -p "+self.pdb[:-4]+".top -o nvt.tpr -maxwarn 50")
		os.system("gmx grompp -f nvt.mdp -c em.gro -p "+self.pdb[:-4]+".top -o nvt.tpr -maxwarn 50")
		
		print("=======================================================")
		print("Running nvt equilibration in gromacs.")
		text_to_run = "/usr/bin/gmx_d" + " mdrun -v -deffnm nvt"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Preparing gromacs NPT equilibration")
		print("gmx grompp -f npt.mdp -c nvt.gro -p "+self.pdb[:-4]+".top -o npt.tpr -maxwarn 50")
		os.system("gmx grompp -f npt.mdp -c nvt.gro -p "+self.pdb[:-4]+".top -o npt.tpr -maxwarn 50")
		
		print("=======================================================")
		print("Running npt equilibration in gromacs.")
		text_to_run = "/usr/bin/gmx_d" + " mdrun -v -deffnm npt"
		os.system(text_to_run)	
		
		print("=======================================================")
		print("Preparing gromacs NPT equilibration")
		print("gmx grompp -f md.mdp -c npt.gro -p "+self.pdb[:-4]+".top -o md.tpr -maxwarn 50")
		os.system("gmx grompp -f md.mdp -c npt.gro -p "+self.pdb[:-4]+".top -o md.tpr -maxwarn 50")
		
		
	def production(self):
		pass 
		
	def organize_files():
		pass
		
	def process_traj():
		pass
		
		

