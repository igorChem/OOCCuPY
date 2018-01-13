#/usr/bin/env python
# -*- coding: utf-8 -*-
# amber_module.py


import os 
from pdb_class import*
from gmx_module import*


class amber_mod:

	def __init__(self,
			pdbfile,
			H_opt=True,
			path="/home/barden/programs/amber16/bin"):

		self.name = pdbfile
		self.nameAf = self.name[:-4] + "_h.pdb"
		self.clear = self.name[:-4] +"_clear_h.pdb"
		self.H_opt = H_opt		
		self.amb_reduce = path + "/reduce"
		self.pdb4amber = path + "/pdb4amber"
		self.tleap = path + "/tleap"
		self.sander = path + "/sander" 
		self.ambpdb = path + "/ambpdb"


		print("========== starting reduce ===========")

		print(self.amb_reduce +" " + self.name + " > " + self.nameAf)

		os.system(self.amb_reduce +" " + self.name + " > " + self.nameAf)
		
		print("========== end of reduce =========")

		print(self.pdb4amber +" -i " + self.nameAf +" -o " + self.clear)
		
		os.system(self.pdb4amber +" -i " + self.nameAf +" -o " + self.clear)

	
	def tleap_call(self):

		tleap_in ="tleap \n" 
		tleap_in += "source leaprc.protein.ff14SB \n"
		tleap_in += "prot = loadPdb " + self.clear + "\n"
		tleap_in += "source leaprc.water.tip3p \n"
		tleap_in += "solvatebox prot TIP3PBOX 10.0 \n"
		tleap_in += "saveamberparm prot prmtop inpcrd\n"
		tleap_in += "quit"

		tleap_file = open('tleap_in','w')
		tleap_file.write(tleap_in)
		tleap_file.close()

		os.system(self.tleap + " -f tleap_in" )


	def sandro(self):

		if self.H_opt==True:
			mini_in =  "minimize hydrogens \n"
			mini_in += " &cntrl \n"
			mini_in += "  imin=1,maxcyc=5000,ncyc=500, \n"
			mini_in += "  ntb=1,ntmin=1,ntpr=100,ntr=1, \n"
			mini_in += "  restraintmask="+'"!@H="'+",igb=0, \n"
			mini_in += "  restraint_wt=500.0,cut=10.0 \n"
			mini_in += "/"

		mini_file = open("minimize.in",'w')
		mini_file.write(mini_in)
		mini_file.close()

	def run_sander(self):

		text_to_run = self.sander +" -O -i minimize.in -o minimize.out -p prmtop -c inpcrd -ref inpcrd -r minimize.rst -inf minimize.mdinfo"
		print(text_to_run)

		os.system(text_to_run)		
	
	def get_pdb(self):
	
		text_to_run = self.ambpdb +" -c inpcrd < minimize.rst > "+ self.name[:-4] +"_opt.pdb"
		
		os.system(text_to_run)

		pdb_opt = protein(self.name,amber=True)
		pdb_opt.pdb_parse(self.name[:-4]+"_opt.pdb")
		pdb_opt.write_pdb(self.name[:-4]+"_opt_final.pdb")

	def mopac_inp(self): 

		opt_pdb = protein(self.name[:-4]+"_opt_final.pdb")
		opt_pdb.pdb_parse()
		opt_pdb.mopac_mop()
	
'''	
a = amber_mod("1a2y_p1.pdb")
a.tleap_call()
a.sandro()
'''
