#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

from pdb_class import *
import parmed as pmd
import glob
import sys,os 


#Cofactor list
cofac_list = ["ATP","ADN","atp"]
#ATP atoms list
atp_list = ["O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]

#names of bin of amber and gromacs
path_amber = "/home/igorchem/programs/amber18/bin"
Reduce     = path_amber + "/reduce -Trim "
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
	
#***********************************************************************
def gromacs_inp():
	'''Function Doc
	Move to gmx_module class file
	Put options to change basic DM runs.
	'''
		
	mdp_file2 =  "integrator    = steep  \n"
	mdp_file2 += "emtol         = 1000.0 \n"
	mdp_file2 += "emstep        = 0.01   \n"
	mdp_file2 += "nsteps        = 100000 \n"		
	mdp_file2 += "nstlist		= 15     \n"
	mdp_file2 += "cutoff-scheme = Verlet \n"
	mdp_file2 += "ns_type		= grid   \n"
	mdp_file2 += "coulombtype	= PME    \n"
	mdp_file2 += "rcoulomb	    = 1.0	 \n"
	mdp_file2 += "rvdw		    = 1.0	 \n"
	mdp_file2 += "pbc		    = xyz    "
	
	equi_file = open("nvt.mdp",'w')
		
	equi_text = " define		= -DPOSRES	; position restrain the protein\n"		
	#Run parameters"
	equi_text += "integrator	= md		; leap-frog integrator\n"
	equi_text += "nsteps		= 50000		; 2 * 50000 = 100 ps\n"
	equi_text += "dt		    = 0.002		; 2 fs\n"
	#Output control
	equi_text += "nstxout		= 1000		; save coordinates every 2.0 ps\n"
	equi_text += "nstvout		= 1000		; save velocities every 2.0 ps\n"
	equi_text += "nstenergy	    = 500		; save energies every 1.0 ps\n"
	equi_text += "nstlog		= 1000		; update log file every 2.0 ps\n"
	#Bond parameters"
	equi_text += "continuation	          = no		; first dynamics run\n"
	equi_text += "constraint_algorithm    = lincs	    ; holonomic constraints \n"
	equi_text += "constraints	          = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
	equi_text += "lincs_iter	          = 1		    ; accuracy of LINCS\n"
	equi_text += "lincs_order	          = 4		    ; also related to accuracy\n"
	#Neighborsearching"
	equi_text += "cutoff-scheme = Verlet\n"
	equi_text += "ns_type		= grid		; search neighboring grid cells\n"
	equi_text += "nstlist	    = 10		; 20 fs, largely irrelevant with Verlet\n"
	equi_text += "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
	equi_text += "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
	#Electrostatics"
	equi_text += "coulombtype    = PME	; Particle Mesh Ewald for long-range electrostatics\n"
	equi_text += "pme_order	     = 4	; cubic interpolation\n"
	equi_text += "fourierspacing = 0.16	; grid spacing for FFT\n"
	#Temperature coupling is on"
	equi_text += "tcoupl	= V-rescale	            ; modified Berendsen thermostat\n"
	equi_text += "tc-grps	= Protein Non-Protein	; two coupling groups - more accurate\n"
	equi_text += "tau_t		= 0.1	  0.1           ; time constant, in ps\n"
	equi_text += "ref_t		= 300 	  300           ; reference temperature, one for each group, in K\n"
	#Pressure coupling is off"
	equi_text += "pcoupl		= no 		; no pressure coupling in NVT\n"
	#Periodic boundary conditions"
	equi_text += "pbc		= xyz		    ; 3-D PBC\n"
	#Dispersion correction"
	equi_text += "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
	#Velocity generation\n"
	equi_text += "gen_vel		= yes	; assign velocities from Maxwell distribution\n"
	equi_text += "gen_temp	= 300		; temperature for Maxwell distribution\n"
	equi_text += "gen_seed	= -1		; generate a random seed\n"

	equi_file.write(equi_text)
	equi_file.close()
		
	mdp_inp = open("em.mdp",'w')
	mdp_inp.write(mdp_file2)
	mdp_inp.close()
	
	equi_file2 = open("npt.mdp",'w')
		
	equi_text2  = "		define		= -DPOSRES	; position restrain the protein\n"
	#Run parameters
	equi_text2 += "integrator	= md		; leap-frog integrator\n"
	equi_text2 += "nsteps		= 50000		; 2 * 50000 = 100 ps\n"
	equi_text2 += "dt		    = 0.002		; 2 fs\n"
	#Output control
	equi_text2 += "nstxout		= 1000		; save coordinates every 1.0 ps\n"
	equi_text2 += "nstvout		= 1000		; save velocities every 1.0 ps\n"
	equi_text2 += "nstenergy	= 500		; save energies every 1.0 ps\n"
	equi_text2 += "nstlog		= 500		; update log file every 1.0 ps\n"
	#Bond parameters
	equi_text2 += "continuation	           = yes		; Restarting after NVT \n"
	equi_text2 += "constraint_algorithm    = lincs	    ; holonomic constraints \n"
	equi_text2 += "constraints	           = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
	equi_text2 += "lincs_iter	           = 1		    ; accuracy of LINCS\n"
	equi_text2 += "lincs_order	           = 4		    ; also related to accuracy\n"
	#Neighborsearching
	equi_text2 += "cutoff-scheme   = Verlet\n"
	equi_text2 += "ns_type		    = grid		; search neighboring grid cells\n"
	equi_text2 += "nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n"
	equi_text2 += "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
	equi_text2 += "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
	#Electrostatics
	equi_text2 += "coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics\n"
	equi_text2 += "pme_order	    = 4		    ; cubic interpolation\n"
	equi_text2 += "fourierspacing	= 0.16		; grid spacing for FFT\n"
	equi_text2 += "; Temperature coupling is on\n"
	equi_text2 += "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
	equi_text2 += "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
	equi_text2 += "tau_t		= 0.1	  0.1	        ; time constant, in ps\n"
	equi_text2 += "ref_t		= 300 	  300	        ; reference temperature, one for each group, in K\n"
	#Pressure coupling is on
	equi_text2 += "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT\n"
	equi_text2 += "pcoupltype	        = isotropic	            ; uniform scaling of box vectors\n"
	equi_text2 += "tau_p		        = 2.0		            ; time constant, in ps\n"
	equi_text2 += "ref_p		        = 1.0		            ; reference pressure, in bar\n"
	equi_text2 += "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n"
	equi_text2 += "refcoord_scaling    = com\n"
	#Periodic boundary conditions
	equi_text2 += "pbc		= xyz		; 3-D PBC\n"
	#Dispersion correction
	equi_text2 += "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
	#Velocity generation
	equi_text2 += "gen_vel		= no		; Velocity generation is off \n"
	
	equi_file2.write(equi_text2)
	equi_file2.close()
	
	equi_file3 = open("md.mdp",'w')
				 
	#Run parameters
	equi_text3  =  "integrator	= md		; leap-frog integrator  \n"
	equi_text3 +=  "nsteps		= 5000000	; 2 * 5000000 = 10000 ps (10 ns)  \n"
	equi_text3 +=  "dt		    = 0.002		; 2 fs  \n"
	#Output control
	equi_text3 +=  "nstxout		        = 5000		; save coordinates every 10.0 ps  \n"
	equi_text3 +=  "nstvout		        = 5000		; save velocities every 10.0 ps  \n"
	equi_text3 +=  "nstenergy	        = 5000		; save energies every 10.0 ps  \n"
	equi_text3 +=  "nstlog		        = 5000		; update log file every 10.0 ps  \n"
	equi_text3 +=  "nstxout-compressed  = 5000      ; save compressed coordinates every 10.0 ps  \n"
	equi_text3 +=  "                                ; nstxout-compressed replaces nstxtcout  \n"
	equi_text3 +=  "compressed-x-grps   = System    ; replaces xtc-grps\n"
	#Bond parameters
	equi_text3 +=  "continuation	        = yes		; Restarting after NPT \n"
	equi_text3 +=  "constraint_algorithm    = lincs	    ; holonomic constraints \n"
	equi_text3 +=  "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
	equi_text3 +=  "lincs_iter	            = 1		    ; accuracy of LINCS\n"
	equi_text3 +=  "lincs_order	            = 4		    ; also related to accuracy\n"
	#Neighborsearching
	equi_text3 +=  "cutoff-scheme   = Verlet\n"
	equi_text3 +=  "ns_type		    = grid		; search neighboring grid cells\n"
	equi_text3 +=  "nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n"
	equi_text3 +=  "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
	equi_text3 +=  "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
	#Electrostatics
	equi_text3 +=  "coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics\n"
	equi_text3 +=  "pme_order	    = 4		    ; cubic interpolation\n"
	equi_text3 +=  "fourierspacing	= 0.16		; grid spacing for FFT\n"
	#Temperature coupling is on
	equi_text3 +=  "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
	equi_text3 +=  "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
	equi_text3 +=  "tau_t		= 0.1	  0.1	        ; time constant, in ps\n"
	equi_text3 +=  "ref_t		= 300 	  300	        ; reference temperature, one for each group, in K\n"
	#Pressure coupling is on
	equi_text3 +=  "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT\n"
	equi_text3 +=  "pcoupltype	        = isotropic	            ; uniform scaling of box vectors\n"
	equi_text3 +=  "tau_p		        = 2.0		            ; time constant, in ps\n"
	equi_text3 +=  "ref_p		        = 1.0		            ; reference pressure, in bar\n"
	equi_text3 +=  "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n"
	#Periodic boundary conditions
	equi_text3 +=  "pbc		= xyz		; 3-D PBC\n"
	#Dispersion correction
	equi_text3 +=  "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
	#Velocity generation
	equi_text3 +=  "gen_vel		= no		; Velocity generation is of\n"
			 
	equi_file3.write(equi_text)
	equi_file3.close()

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
		print("Removing any hydrogens in the ligand pdb") 
		print(Reduce + self.lig+".pdb > " + self.lig+"_h.pdb")
		os.system(Reduce + self.lig+".pdb > " + self.lig+"_h.pdb")
		
		if self.lig in cofac_list:
			self.lig = self.lig+"_h"
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
			if self.lig+".frcmod" in f:
				print("Found parameters for " + self.lig)
				print("FRCMOD file found for this ligand, antechamber parametrization will be skipped!") 
				par = True					
			if not par:
				print("===================================================")
				print("Run ANTECHAMBER:")
				print(antech+" -i "+self.lig+"_h.pdb -fi pdb -o "+self.lig+".mol2 -fo mol2 -c bcc -nc "+chg)
				os.system(antech + " -i " + self.lig+"_h.pdb -fi pdb -o " + self.lig+".mol2  -fo mol2 -c bcc -nc "+chg )
				os.system("rm ANTECHAMBER*")
		
				print("===================================================")
				print("Run Pamchek and generate frcmod")				
				print(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
				print(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
				os.system(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
		
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
		print(Reduce + self.current_pdb+  " > " +self.current_pdb[:-4]+"_p.pdb") 
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
		
		

