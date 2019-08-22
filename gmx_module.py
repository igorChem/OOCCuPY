#/usr/bin/env python
# -*- coding: utf-8 -*-
# gmx_module.py


#load modules

import os 
from pdb_class import*
from xyz_class import*

def gromacs_inp():
	'''Function Doc
	Put options to change basic DM runs.
		time in ns 
		temperature in K
	'''
	
	em_file    = open("em.mdp",'w')		
	nvt_file   = open("nvt.mdp",'w')
	npt_file   = open("npt.mdp",'w')
	md_file    = open("md.mdp",'w')
		
	em_text =  "integrator    = steep  \n"
	em_text += "emtol         = 1000.0 \n"
	em_text += "emstep        = 0.01   \n"
	em_text += "nsteps        = 100000 \n"		
	em_text += "nstlist		= 15     \n"
	em_text += "cutoff-scheme = Verlet \n"
	em_text += "ns_type		= grid   \n"
	em_text += "coulombtype	= PME    \n"
	em_text += "rcoulomb	    = 1.0	 \n"
	em_text += "rvdw		    = 1.0	 \n"
	em_text += "pbc		    = xyz    "
	
		
	nvt_text = " define		= -DPOSRES	; position restrain the protein\n"		
	#Run parameters"
	nvt_text += "integrator	= md		; leap-frog integrator\n"
	nvt_text += "nsteps		= 50000		; 2 * 50000 = 100 ps\n"
	nvt_text += "dt		    = 0.002		; 2 fs\n"
	#Output control
	nvt_text += "nstxout		= 1000		; save coordinates every 2.0 ps\n"
	nvt_text += "nstvout		= 1000		; save velocities every 2.0 ps\n"
	nvt_text += "nstenergy	    = 500		; save energies every 1.0 ps\n"
	nvt_text += "nstlog		= 1000		; update log file every 2.0 ps\n"
	#Bond parameters"
	nvt_text += "continuation	          = no		; first dynamics run\n"
	nvt_text += "constraint_algorithm    = lincs	    ; holonomic constraints \n"
	nvt_text += "constraints	          = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
	nvt_text += "lincs_iter	          = 1		    ; accuracy of LINCS\n"
	nvt_text += "lincs_order	          = 4		    ; also related to accuracy\n"
	#Neighborsearching"
	nvt_text += "cutoff-scheme = Verlet\n"
	nvt_text += "ns_type		= grid		; search neighboring grid cells\n"
	nvt_text += "nstlist	    = 10		; 20 fs, largely irrelevant with Verlet\n"
	nvt_text += "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
	nvt_text += "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
	#Electrostatics"
	nvt_text += "coulombtype    = PME	; Particle Mesh Ewald for long-range electrostatics\n"
	nvt_text += "pme_order	     = 4	; cubic interpolation\n"
	nvt_text += "fourierspacing = 0.16	; grid spacing for FFT\n"
	#Temperature coupling is on"
	nvt_text += "tcoupl	= V-rescale	            ; modified Berendsen thermostat\n"
	nvt_text += "tc-grps	= Protein Non-Protein	; two coupling groups - more accurate\n"
	nvt_text += "tau_t		= 0.1	  0.1           ; time constant, in ps\n"
	nvt_text += "ref_t		= 300 	  300           ; reference temperature, one for each group, in K\n"
	#Pressure coupling is off"
	nvt_text += "pcoupl		= no 		; no pressure coupling in NVT\n"
	#Periodic boundary conditions"
	nvt_text += "pbc		= xyz		    ; 3-D PBC\n"
	#Dispersion correction"
	nvt_text += "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
	#Velocity generation\n"
	nvt_text += "gen_vel		= yes	; assign velocities from Maxwell distribution\n"
	nvt_text += "gen_temp	= 300		; temperature for Maxwell distribution\n"
	nvt_text += "gen_seed	= -1		; generate a random seed\n"

				
	npt_text2  = "		define		= -DPOSRES	; position restrain the protein\n"
	#Run parameters
	npt_text2 += "integrator	= md		; leap-frog integrator\n"
	npt_text2 += "nsteps		= 50000		; 2 * 50000 = 100 ps\n"
	npt_text2 += "dt		    = 0.002		; 2 fs\n"
	#Output control
	npt_text2 += "nstxout		= 1000		; save coordinates every 1.0 ps\n"
	npt_text2 += "nstvout		= 1000		; save velocities every 1.0 ps\n"
	npt_text2 += "nstenergy	= 500		; save energies every 1.0 ps\n"
	npt_text2 += "nstlog		= 500		; update log file every 1.0 ps\n"
	#Bond parameters
	npt_text2 += "continuation	           = yes		; Restarting after NVT \n"
	npt_text2 += "constraint_algorithm    = lincs	    ; holonomic constraints \n"
	npt_text2 += "constraints	           = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
	npt_text2 += "lincs_iter	           = 1		    ; accuracy of LINCS\n"
	npt_text2 += "lincs_order	           = 4		    ; also related to accuracy\n"
	#Neighborsearching
	npt_text2 += "cutoff-scheme   = Verlet\n"
	npt_text2 += "ns_type		    = grid		; search neighboring grid cells\n"
	npt_text2 += "nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n"
	npt_text2 += "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
	npt_text2 += "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
	#Electrostatics
	npt_text2 += "coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics\n"
	npt_text2 += "pme_order	    = 4		    ; cubic interpolation\n"
	npt_text2 += "fourierspacing	= 0.16		; grid spacing for FFT\n"
	npt_text2 += "; Temperature coupling is on\n"
	npt_text2 += "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
	npt_text2 += "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
	npt_text2 += "tau_t		= 0.1	  0.1	        ; time constant, in ps\n"
	npt_text2 += "ref_t		= 310 	  310	        ; reference temperature, one for each group, in K\n"
	#Pressure coupling is on
	npt_text2 += "pcoupl		      = Parrinello-Rahman	    ; Pressure coupling on in NPT\n"
	npt_text2 += "pcoupltype	      = isotropic	            ; uniform scaling of box vectors\n"
	npt_text2 += "tau_p		          = 2.0		            ; time constant, in ps\n"
	npt_text2 += "ref_p		          = 1.0		            ; reference pressure, in bar\n"
	npt_text2 += "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n"
	npt_text2 += "refcoord_scaling    = com\n"
	#Periodic boundary conditions
	npt_text2 += "pbc		= xyz		; 3-D PBC\n"
	#Dispersion correction
	npt_text2 += "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
	#Velocity generation
	npt_text2 += "gen_vel		= no		; Velocity generation is off \n"
	
			 
	#Run parameters
	md_text  =  "integrator	= md		; leap-frog integrator  \n"
	md_text +=  "nsteps		= 5000000	; 2 * 5000000 = 10000 ps (10 ns)  \n"
	md_text +=  "dt		    = 0.002		; 2 fs  \n"
	#Output control
	md_text +=  "nstxout		        = 5000		; save coordinates every 10.0 ps  \n"
	md_text +=  "nstvout		        = 5000		; save velocities every 10.0 ps  \n"
	md_text +=  "nstenergy	        = 5000		; save energies every 10.0 ps  \n"
	md_text +=  "nstlog		        = 5000		; update log file every 10.0 ps  \n"
	md_text +=  "nstxout-compressed  = 5000      ; save compressed coordinates every 10.0 ps  \n"
	md_text +=  "                                ; nstxout-compressed replaces nstxtcout  \n"
	md_text +=  "compressed-x-grps   = System    ; replaces xtc-grps\n"
	#Bond parameters
	md_text +=  "continuation	        = yes		; Restarting after NPT \n"
	md_text +=  "constraint_algorithm    = lincs	    ; holonomic constraints \n"
	md_text +=  "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
	md_text +=  "lincs_iter	            = 1		    ; accuracy of LINCS\n"
	md_text +=  "lincs_order	            = 4		    ; also related to accuracy\n"
	#Neighborsearching
	md_text +=  "cutoff-scheme   = Verlet\n"
	md_text +=  "ns_type		    = grid		; search neighboring grid cells\n"
	md_text +=  "nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n"
	md_text +=  "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
	md_text +=  "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
	#Electrostatics
	md_text +=  "coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics\n"
	md_text +=  "pme_order	    = 4		    ; cubic interpolation\n"
	md_text +=  "fourierspacing	= 0.16		; grid spacing for FFT\n"
	#Temperature coupling is on
	md_text +=  "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
	md_text +=  "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
	md_text +=  "tau_t		= 0.1	  0.1	        ; time constant, in ps\n"
	md_text +=  "ref_t		= 310 	  310	        ; reference temperature, one for each group, in K\n"
	#Pressure coupling is on
	md_text +=  "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT\n"
	md_text +=  "pcoupltype	        = isotropic	            ; uniform scaling of box vectors\n"
	md_text +=  "tau_p		        = 2.0		            ; time constant, in ps\n"
	md_text +=  "ref_p		        = 1.0		            ; reference pressure, in bar\n"
	md_text +=  "compressibility    = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n"
	#Periodic boundary conditions
	md_text +=  "pbc		= xyz		; 3-D PBC\n"
	#Dispersion correction
	md_text +=  "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
	#Velocity generation
	md_text +=  "gen_vel		= no		; Velocity generation is of\n"

	em_file.write(em_text)
	nvt_file.write(nvt_text)
	npt_file.write(npt_text)
	md_file.write(md_text)
	
	em_file.close()
	nvt_file.close()
	npt_file.close()
	md_file.close()
	

class min_prot:

	def __init__(self,
				 pdb_code):

		self.pdb_code   = pdb_code
		self.net_charge = 0
		self.protein    = pdb_code[:-4]		

	def top_init(self):

		text_to_run = "/usr/bin/gmx" + " pdb2gmx -ignh -f " + self.protein + ".pdb -o " + self.protein + "_P.pdb -water tip3p << EOF \n"
		text_to_run += "5 \n"
		text_to_run += "EOF"

		os.system(text_to_run)

	def solvate(self):

		text_to_run = "/usr/bin/gmx" + " editconf -f " + self.protein + "_P.pdb -o " + self.protein +"_NB.pdb -c -d 1.0 -bt cubic"

		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " solvate -cp " +  self.protein + "_NB.pdb -cs spc216.gro -o " + self.protein + "_solv.pdb -p topol.top"

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

		text_to_run = "/usr/bin/gmx" + " grompp -f ions.mdp -c " + self.protein + "_solv.pdb -p topol.top -o ions.tpr -maxwarn 50<< EOF\n"
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
		mdp_file += "nstlist		 = 2	\n "
		mdp_file += "cutoff-scheme   = Verlet \n"
		mdp_file += "ns_type		 = grid\n "
		mdp_file += "coulombtype	 = PME\n "
		mdp_file += "rcoulomb	    = 1.0	\n"
		mdp_file += "rvdw		    = 1.0	\n"
		mdp_file += "pbc		    = xyz\n "

		mdp_inp = open("minim.mdp",'w')
		mdp_inp.write(mdp_file)
		mdp_inp.close()

		text_to_run = "/usr/bin/gmx" +" grompp -f minim.mdp -c " + self.protein + "_solv_ions.pdb -p topol.top -o em.tpr -maxwarn 50"
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

	def equilibration(self):

		equi_file = open("nvt.mdp",'w')
		
		equi_text = " define		= -DPOSRES	; position restrain the protein\n"
		
		equi_text += "		; Run parameters \n"
		equi_text += "integrator	= md		; leap-frog integrator\n"
		equi_text += "nsteps		= 50000		; 2 * 50000 = 100 ps\n"
		equi_text += "dt		    = 0.002		; 2 fs\n"
		equi_text += "; Output control\n"
		equi_text += "nstxout		= 500		; save coordinates every 1.0 ps\n"
		equi_text += "nstvout		= 500		; save velocities every 1.0 ps\n"
		equi_text += "nstenergy	= 500		; save energies every 1.0 ps\n"
		equi_text += "nstlog		= 500		; update log file every 1.0 ps\n"
		equi_text += "; Bond parameters\n"
		equi_text += "continuation	        = no		; first dynamics run\n"
		equi_text += "constraint_algorithm    = lincs	    ; holonomic constraints \n"
		equi_text += "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
		equi_text += "lincs_iter	            = 1		    ; accuracy of LINCS\n"
		equi_text += "lincs_order	            = 4		    ; also related to accuracy\n"
		equi_text += "; Neighborsearching\n"
		equi_text += "cutoff-scheme   = Verlet\n"
		equi_text += "ns_type		    = grid		; search neighboring grid cells\n"
		equi_text += "nstlist		    = 10		; 20 fs, largely irrelevant with Verlet\n"
		equi_text += "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
		equi_text += "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
		equi_text += "; Electrostatics\n"
		equi_text += "coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics\n"
		equi_text += "pme_order	    = 4		; cubic interpolation\n"
		equi_text += "fourierspacing	= 0.16	; grid spacing for FFT\n"
		equi_text += "; Temperature coupling is on\n"
		equi_text += "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
		equi_text += "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
		equi_text += "tau_t		= 0.1	  0.1           ; time constant, in ps\n"
		equi_text += "ref_t		= 300 	  300           ; reference temperature, one for each group, in K\n"
		equi_text += "; Pressure coupling is off\n"
		equi_text += "pcoupl		= no 		; no pressure coupling in NVT\n"
		equi_text += "; Periodic boundary conditions\n"
		equi_text += "pbc		= xyz		    ; 3-D PBC\n"
		equi_text += "; Dispersion correction\n"
		equi_text += "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
		equi_text += "; Velocity generation\n"
		equi_text += "gen_vel		= yes		; assign velocities from Maxwell distribution\n"
		equi_text += "gen_temp	= 300		; temperature for Maxwell distribution\n"
		equi_text += "gen_seed	= -1		; generate a random seed\n"

		equi_file.write(equi_text)
		equi_file.close()
		
		text_to_run = "/usr/bin/gmx" +" grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " mdrun -deffnm nvt"
		os.system(text_to_run)

		equi_file = open("npt.mdp",'w')
		
		equi_text  = " title		= OPLS Lysozyme NPT equilibration \n"
		equi_text += "		define		= -DPOSRES	; position restrain the protein\n"
		equi_text += "; Run parameters\n"
		equi_text += "integrator	= md		; leap-frog integrator\n"
		equi_text += "nsteps		= 50000		; 2 * 50000 = 100 ps\n"
		equi_text += "dt		    = 0.002		; 2 fs\n"
		equi_text += "; Output control\n"
		equi_text += "nstxout		= 500		; save coordinates every 1.0 ps\n"
		equi_text += "nstvout		= 500		; save velocities every 1.0 ps\n"
		equi_text += "nstenergy	= 500		; save energies every 1.0 ps\n"
		equi_text += "nstlog		= 500		; update log file every 1.0 ps\n"
		equi_text += "; Bond parameters\n"
		equi_text += "continuation	        = yes		; Restarting after NVT \n"
		equi_text += "constraint_algorithm    = lincs	    ; holonomic constraints \n"
		equi_text += "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
		equi_text += "lincs_iter	            = 1		    ; accuracy of LINCS\n"
		equi_text += "lincs_order	            = 4		    ; also related to accuracy\n"
		equi_text += "; Neighborsearching\n"
		equi_text += "cutoff-scheme   = Verlet\n"
		equi_text += "ns_type		    = grid		; search neighboring grid cells\n"
		equi_text += "nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n"
		equi_text += "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
		equi_text += "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
		equi_text += "; Electrostatics\n"
		equi_text += "coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics\n"
		equi_text += "pme_order	    = 4		    ; cubic interpolation\n"
		equi_text += "fourierspacing	= 0.16		; grid spacing for FFT\n"
		equi_text += "; Temperature coupling is on\n"
		equi_text += "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
		equi_text += "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
		equi_text += "tau_t		= 0.1	  0.1	        ; time constant, in ps\n"
		equi_text += "ref_t		= 300 	  300	        ; reference temperature, one for each group, in K\n"
		equi_text += "; Pressure coupling is on\n"
		equi_text += "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT\n"
		equi_text += "pcoupltype	        = isotropic	            ; uniform scaling of box vectors\n"
		equi_text += "tau_p		        = 2.0		            ; time constant, in ps\n"
		equi_text += "ref_p		        = 1.0		            ; reference pressure, in bar\n"
		equi_text += "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n"
		equi_text += "refcoord_scaling    = com\n"
		equi_text += "; Periodic boundary conditions\n"
		equi_text += "pbc		= xyz		; 3-D PBC\n"
		equi_text += "; Dispersion correction\n"
		equi_text += "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
		equi_text += "; Velocity generation\n"
		equi_text += "gen_vel		= no		; Velocity generation is off \n"
		
		equi_file.write(equi_text)
		equi_file.close()
		
		text_to_run = "/usr/bin/gmx" +" grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " mdrun  -deffnm npt"
		os.system(text_to_run)
		
	def production(self):
		
		equi_file3 = open("md.mdp",'w')
				 
		equi_text3 =  "title		= OPLS Lysozyme MD simulation  \n"
		equi_text3 +=  "		; Run parameters  \n"
		equi_text3 +=  "integrator	= md		; leap-frog integrator  \n"
		equi_text3 +=  "nsteps		= 5000000	; 2 * 5000000 = 10000 ps (10 ns)  \n"
		equi_text3 +=  "dt		    = 0.002		; 2 fs  \n"
		equi_text3 +=  "; Output control  \n"
		equi_text3 +=  "nstxout		        = 5000		; save coordinates every 10.0 ps  \n"
		equi_text3 +=  "nstvout		        = 5000		; save velocities every 10.0 ps  \n"
		equi_text3 +=  "nstenergy	        = 5000		; save energies every 10.0 ps  \n"
		equi_text3 +=  "nstlog		        = 5000		; update log file every 10.0 ps  \n"
		equi_text3 +=  "nstxout-compressed  = 5000      ; save compressed coordinates every 10.0 ps  \n"
		equi_text3 +=  "                                ; nstxout-compressed replaces nstxtcout  \n"
		equi_text3 +=  "compressed-x-grps   = System    ; replaces xtc-grps\n"
		equi_text3 +=  "; Bond parameters\n"
		equi_text3 +=  "continuation	        = yes		; Restarting after NPT \n"
		equi_text3 +=  "constraint_algorithm    = lincs	    ; holonomic constraints \n"
		equi_text3 +=  "constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained\n"
		equi_text3 +=  "lincs_iter	            = 1		    ; accuracy of LINCS\n"
		equi_text3 +=  "lincs_order	            = 4		    ; also related to accuracy\n"
		equi_text3 +=  "; Neighborsearching\n"
		equi_text3 +=  "cutoff-scheme   = Verlet\n"
		equi_text3 +=  "ns_type		    = grid		; search neighboring grid cells\n"
		equi_text3 +=  "nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme\n"
		equi_text3 +=  "rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)\n"
		equi_text3 +=  "rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)\n"
		equi_text3 +=  "; Electrostatics\n"
		equi_text3 +=  "coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics\n"
		equi_text3 +=  "pme_order	    = 4		    ; cubic interpolation\n"
		equi_text3 +=  "fourierspacing	= 0.16		; grid spacing for FFT\n"
		equi_text3 +=  "; Temperature coupling is on\n"
		equi_text3 +=  "tcoupl		= V-rescale	            ; modified Berendsen thermostat\n"
		equi_text3 +=  "tc-grps		= Protein Non-Protein	; two coupling groups - more accurate\n"
		equi_text3 +=  "tau_t		= 0.1	  0.1	        ; time constant, in ps\n"
		equi_text3 +=  "ref_t		= 300 	  300	        ; reference temperature, one for each group, in K\n"
		equi_text3 +=  "; Pressure coupling is on\n"
		equi_text3 +=  "pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT\n"
		equi_text3 +=  "pcoupltype	        = isotropic	            ; uniform scaling of box vectors\n"
		equi_text3 +=  "tau_p		        = 2.0		            ; time constant, in ps\n"
		equi_text3 +=  "ref_p		        = 1.0		            ; reference pressure, in bar\n"
		equi_text3 +=  "compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1\n"
		equi_text3 +=  "; Periodic boundary conditions\n"
		equi_text3 +=  "pbc		= xyz		; 3-D PBC\n"
		equi_text3 +=  "; Dispersion correction\n"
		equi_text3 +=  "DispCorr	= EnerPres	; account for cut-off vdW scheme\n"
		equi_text3 +=  "; Velocity generation\n"
		equi_text3 +=  "gen_vel		= no		; Velocity generation is of\n"
				 
		equi_file3.write(equi_text)
		equi_file3.close()
		
		text_to_run = "/usr/bin/gmx" +" grompp -f md.mdp -c npt.gro -p topol.top -o md_0_1.tpr"
		os.system(text_to_run)

		text_to_run = "/usr/bin/gmx" + " mdrun  -deffnm md_0_1"
		os.system(text_to_run)
		
		text_to_run = "/usr/bin/gmx" + " trjconv -f  md_0_1.trr -s md_0_1.tpr -pbc mol -dt 10 -o output.pdb"
		os.system(text_to_run)
		
	def run(self):
		self.top_init()
		self.solvate()
		self.add_ions()
		self.Minimization()
		self.equilibration()
		self.production()
		#self.write_minStruct()
		#self.rewrite_pdb()	
		
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

