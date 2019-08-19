#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

from pdb_class import *
import parmed as pmd
import glob
import sys,os 


#names of bin of amber and gromacs
path_amber = "/home/igorchem/programs/amber18/bin"
Reduce = path_amber + "/reduce "
pdb4 = path_amber + "/pdb4amber "
antech = path_amber + "/antechamber "
parmchk = path_amber + "/parmchk2 "
tleap = path_amber + "/tleap "

def my_replace(fl,old,new):
	with open(fl, 'r+') as f:
			s = f.read()
			s = s.replace(old, new)
			f.write(s)



def pdb_cat(pdb1,pdb2):
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
	pdb_res.close()
	
def gromacs_inp():

	mdp_file =  "integrator = steep \n"
	mdp_file += "emtol = 1000.0 \n"
	mdp_file += "emstep = 0.01 \n"
	mdp_file += "nsteps = 50000 \n"		
	mdp_file += "nstlist		 = 1	\n "
	mdp_file += "cutoff-scheme   = Verlet \n"
	mdp_file += "ns_type		 = grid\n "
	mdp_file += "coulombtype	 = PME\n "
	mdp_file += "rcoulomb	    = 1.0	\n"
	mdp_file += "rvdw		    = 1.0	\n"
	mdp_file += "pbc		    = xyz\n "
	mdp_file =  "integrator = steep \n"
		
		
	mdp_file2 = "emtol = 1000.0 \n"
	mdp_file2 += "emstep = 0.0051 \n"
	mdp_file2 += "nsteps = 50000 \n"		
	mdp_file2 += "nstlist		 = 2	\n "
	mdp_file2 += "cutoff-scheme   = Verlet \n"
	mdp_file2 += "ns_type		 = grid\n "
	mdp_file2 += "coulombtype	 = PME\n "
	mdp_file2 += "rcoulomb	    = 1.0	\n"
	mdp_file2 += "rvdw		    = 1.0	\n"
	mdp_file2 += "pbc		    = xyz\n "



	ls = os.listdir(".")	
	for fl in ls:		
		if not fl == "ions.mdp":
			mdp_inp = open("ions.mdp",'w')
			mdp_inp.write(mdp_file)
			mdp_inp.close()
		if not fl == "em.mdp":
			mdp_inp = open("em.mdp",'w')
			mdp_inp.write(mdp_file2)
			mdp_inp.close()
			
class md_prep:
	
	def __init__(self,pdb):
		self.pdb = pdb
		self.current_pdb = pdb 
		self.lig = "none"
		self.net_charge = 0
	def prepare_lig(self,lign,chg=0,rwat=True):
		self.lig = lign
		pdb = protein(self.pdb)
		pdb.pdb_parse()
		if rwat:
			pdb.remove_waters()
		lig = []
	
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
			
		print("grep -v EFZ "+self.pdb+" > "+ self.pdb[:-4] +"_wl.pdb ")
		os.system("grep -v "+ self.lig +" "+self.pdb+" > "+ self.pdb[:-4] +"_wl.pdb ")
		
		os.system(Reduce + self.lig+".pdb > " + self.lig+"_h.pdb")
		
		os.system(antech + " -i " + self.lig+"_h.pdb -fi pdb -o " + self.lig+".mol2  -fo mol2 -c bcc -nc "+chg )
		os.system("rm ANTECHAMBER*")
		
		print(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
		os.system(parmchk + " -i "+ self.lig+".mol2 -f mol2 -o " + self.lig+".frcmod")
			
		tleap_in = "source leaprc.gaff \n"
		tleap_in += self.lig + " = loadmol2 " + self.lig +".mol2\n"
		tleap_in += "check " + self.lig + "\n"
		tleap_in += "loadamberparams " + self.lig + ".frcmod\n"
		tleap_in += "saveoff " + self.lig + " " + self.lig +".lib \n"		
		tleap_in += "saveamberparm prot prmtop inpcrd\n"
		tleap_in += "quit"
		
		tleap_file = open('tleap_in','w')
		tleap_file.write(tleap_in)
		tleap_file.close()

		os.system(tleap + " -f tleap_in" )
		
		self.current_pdb = self.pdb[:-4] +"_wl.pdb"
		
	def build_complex(self):
		
		os.system(Reduce +"-Trim "+ self.current_pdb+  " > " +self.current_pdb[:-4]+"_p.pdb")
		os.system(pdb4 +  self.current_pdb[:-4]+"_p.pdb"+ " > " +self.current_pdb[:-4]+"_c.pdb")
		
		pdb_cat(self.current_pdb[:-4]+"_c.pdb",self.lig+"_h.pdb")
		os.rename(self.current_pdb[:-4]+"_c_comp.pdb",self.pdb[:-4]+"_comp.pdb")
		self.current_pdb = self.pdb[:-4]+"_comp.pdb"
		
		tleap_in = "source leaprc.protein.ff14SB \n"
		tleap_in += "source leaprc.gaff \n"
		tleap_in += "loadamberparams " + self.lig + ".frcmod\n"
		tleap_in += "loadoff " + self.lig + ".lib\n"
		tleap_in += "complex = loadPdb " + self.current_pdb+ " \n"
		tleap_in += "source leaprc.water.tip3p \n"
		tleap_in += "solvatebox complex TIP3PBOX 10.0 \n"
		tleap_in += "savePdb complex "+self.current_pdb+"\n"
		tleap_in += "saveamberparm complex " +self.pdb[:-4]+".prmtop "+ self.pdb[:-4] +".inpcrd\n"
		tleap_in += "quit"

		tleap_file = open('tleap_in','w')
		tleap_file.write(tleap_in)
		tleap_file.close()

		os.system(tleap + " -f tleap_in" )
		
		try:
			ap = pmd.load_file(self.pdb[:-4]+".prmtop",self.pdb[:-4] +".inpcrd")
			ap.save(self.pdb[:-4]+".top")
			ap.save(self.pdb[:-4]+".gro")
		except:
			print("Topologies and coordinates conversions for gromacs was not done this time. ")

	def prepare_gromacs(self):
		
		topol_file = open(self.pdb[:-4]+".top",'r')
		for line in topol_file:
			line2 = line.split()			
			if len(line2) == 8:
				if line2[1] == "residue":
					if line2[7] == "1.0" or line2[7] == "+1.0":
						self.net_charge +=1
					elif line2[7] == "2.0" or line2[7] == "+2.0":
						self.net_charge +=2
					elif line2[7] == "-1.0":
						self.net_charge -=1
					elif line2[7] == "-2.0":
						self.net_charge -=2
					print(self.net_charge)
		NN = 0
		NP = 0

		if self.net_charge > 0:
			NN = self.net_charge
		elif self.net_charge < 0:
			NP = self.net_charge
		
		gromacs_inp()
		my_replace(self.current_pdb,"WAT","SOL")
		my_replace(self.pdb[:-4]+".top","WAT","SOL")
		my_replace(self.pdb[:-4]+".gro","WAT","SOL")
		
		
		text_to_run = "/usr/bin/gmx" + " grompp -f ions.mdp -c " + self.current_pdb+" -p "+ self.pdb[:-4] +".top -o ions.tpr -maxwarn 50<< EOF\n"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " genion -s ions.tpr -o " + self.current_pdb[:-4] + "_i.pdb -p " +self.pdb[:-4] +".top -pname NA -nname CL -nn {0} -np {1}".format(NN,NP)
		os.system(text_to_run)	
		
	
	def min_gromacs(self):
		text_to_run = "/usr/bin/gmx" +" grompp -f em.mdp -c "+self.current_pdb[:-4] + "_i.pdb -p "+ self.pdb[:-4] +".top -o em.tpr -maxwarn 50"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " mdrun -v -deffnm em"
		os.system(text_to_run)		
		
	def equilibration(self):
		pass
	
	def production(self):
		pass 
		
		
	
if __name__ == "__main__":
	a = md_prep(sys.argv[1])
	a.prepare_lig(sys.argv[2],sys.argv[3])
	a.build_complex()
	a.prepare_gromacs()
	a.min_gromacs()
	'''
	try:
		a.prepare_lig(sys.argv[2],sys.argv[3])
	except:
		print("Some Error!")
	'''
	
