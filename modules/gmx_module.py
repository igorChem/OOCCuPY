#/usr/bin/env python
# -*- coding: utf-8 -*-
# gmx_module.py


#load modules

import os 
from pdb_class import*
from xyz_class import*



class min_prot:

	def __init__(self,
				 pdb_code):

		self.pdb_code = pdb_code
		self.net_charge = 0
		self.protein = pdb_code[:-4]

		

	def top_init(self):

		text_to_run = "/usr/bin/gmx" + " pdb2gmx -ignh -f " + self.protein + ".pdb -o " + self.protein + "_processed.pdb -water spce << EOF \n"
		text_to_run += "1 \n"
		text_to_run += "EOF"

		os.system(text_to_run)

	def solvate(self):

		text_to_run = "/usr/bin/gmx" + " editconf -f " + self.protein + "_processed.pdb -o " + self.protein +"_newbox.pdb -c -d 1.0 -bt cubic"

		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " solvate -cp " +  self.protein + "_newbox.pdb -cs spc216.gro -o " + self.protein + "_solv.pdb -p topol.top"

		os.system(text_to_run)

		topol_file = open("topol.top",'r')
		for line in topol_file:
			line2 = line.split()			
			if len(line2) == 8:
				if line2[1] == "residue":
					if line2[7] == "+1.0":
						self.net_charge +=1
					elif line2[7] == "+2.0":
						self.net_charge +=2
					elif line2[7] == "-1.0":
						self.net_charge -=1
					elif line2[7] == "-2.0":
						self.net_charge -=2		

	def add_ions(self):

		NN = 0
		NP = 0

		if self.net_charge > 0:
			NN = self.net_charge
		elif self.net_charge < 0:
			NP = self.net_charge

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

		mdp_inp = open("ions.mdp",'w')
		mdp_inp.write(mdp_file)
		mdp_inp.close()

		text_to_run = "/usr/bin/gmx" + " grompp -f ions.mdp -c " + self.protein + "_solv.pdb -p topol.top -o ions.tpr << EOF\n"
		text_to_run += "13 \n"
		text_to_run += "EOF"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " genion -s ions.tpr -o " + self.protein + "_solv_ions.pdb -p topol.top -pname NA -nname CL -nn {0} -np {1}".format(NN,NP)
		os.system(text_to_run)

	def Minimization(self):	

		mdp_file =  "integrator = steep \n"
		mdp_file += "emtol = 1000.0 \n"
		mdp_file += "emstep = 0.0051 \n"
		mdp_file += "nsteps = 50000 \n"		
		mdp_file += "nstlist		 = 1	\n "
		mdp_file += "cutoff-scheme   = Verlet \n"
		mdp_file += "ns_type		 = grid\n "
		mdp_file += "coulombtype	 = PME\n "
		mdp_file += "rcoulomb	    = 1.0	\n"
		mdp_file += "rvdw		    = 1.0	\n"
		mdp_file += "pbc		    = xyz\n "

		mdp_inp = open("minim.mdp",'w')
		mdp_inp.write(mdp_file)
		mdp_inp.close()

		text_to_run = "/usr/bin/gmx" +" grompp -f minim.mdp -c " + self.protein + "_solv_ions.pdb -p topol.top -o em.tpr"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " mdrun -v -deffnm em"
		os.system(text_to_run)		
		
	def write_minStruct(self):
		
		emgro = open("em.gro",'r')

		res_num = []
		res_name = []
		atom_name = [] 
		atomN = [] 
		atom_sym = []
		xcoord = []
		ycoord = []
		zcoord = []

		word1 = "NUCLEAR"
		word2 = "Protein"

		couter = 0
		for line in emgro:
			print(line)
			line2 = line.split()
			if line2[0] == word1:
				continue			
			if len(line2) == 6:
				if not line2[0][-3:] == "SOL":				
					res_name.append(line2[0][-3:])
					res_num.append(line2[0][:-3])
					atom_name.append(line2[1])
					atomN.append(line2[2])
					xxx = round(10*(float(line2[3])),5)
					yyy = round(10*(float(line2[4])),5)
					zzz = round(10*(float(line2[5])),5)
					xcoord.append(xxx)
					ycoord.append(yyy)
					zcoord.append(zzz)
					atom_sym.append(line2[1][0])

		pdb_file = open(self.protein+"_min.pdb",'w')
		
		xyz_file = open(self.protein+"_min.xyz",'w')

		pdb_text = ""
		xyz_text = "{0} \n\n".format(len(res_num))

		for i in range(len(res_num)):
			pdb_text += "ATOM      {0:>4}  {1:<5} {2:>4} A {3:>4} {4:>5} {5:>5} {6:>5}  1.00  0.00  {7:>8} \n".format(atomN[i],atom_name[i],res_name[i],res_num[i],xcoord[i],ycoord[i],zcoord[i],atom_sym[i])

		for i in range(len(res_num)):
			xyz_text += "{0} {1} {2} {3} \n".format(atom_sym[i],xcoord[i],ycoord[i],zcoord[i])

		pdb_file.write(pdb_text)
		pdb_file.close()

		xyz_file.write(xyz_text)
		xyz_file.close()

	def rewrite_pdb(self):

		a = protein(self.protein)
		a.pdb_parse(self.protein+"_min.pdb")
		a.residue_def(reorg=True)
		a.write_pdb(self.protein+"_.pdb")		

	def run(self):
		self.top_init()
		self.solvate()
		self.add_ions()
		self.Minimization()
		self.write_minStruct()
		self.rewrite_pdb()	
		
		os.system("/usr/bin/pymol "+self.protein+"_min.xyz")




'''
a = min_prot("1AKI.pdb")
#a.top_init()
#a.solvate()
#a.add_ions()
#a.Minimization()
#a.write_minStruct()
a.run()
'''

