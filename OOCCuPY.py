#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

from pdb_class import *
from mopac_module import *
import glob
import sys,os 

'''
1. Make input for mopac
2. Make input for gamess
3. Make input for orca
4. make input for FMO
4. make input for primordia
5. Make protein minimizations with gromacs
6. Make simple molecular dynamics with gromacs
7. Find the most probable Structure in molecular dynamics
'''

sem = "PM7"
moz = ""
mgf = False

def inp_mopac_from_pdb(pdb):
	a = protein(pdb)
	a.pdb_parse()
	a.write_xyz()
	a = mopac_inp(xyzfile=pdb[:-4]+".xyz",charge=0,multi=0,inpnam=pdb[:-4]+".mop",method="PM7",mgf=mgf)
	a.write_mop()

def primordia_inp(option=3,program="mopac",lh="potential_fukui",gridn=0,eband=5,bandm="BD",norb=100):

	f = open("input_pri",'w')
	text ="eband {0} pymols\n".format(eband)


	if 	 option == '1':
		if program   == "mopac":
			aux = glob.glob("*aux")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)
		elif program == "gamess":
			aux = glob.glob("*log")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)
		elif program == "orca":
			aux = glob.glob("*out")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)
		elif program == "gaussian":
			aux = glob.glob("*fchk")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)		
	elif option == '2':
		if program   == "mopac":
			aux = glob.glob("*aux")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		elif program == "gamess":
			aux = glob.glob("*log")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		elif program == "orca":
			aux = glob.glob("*out")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		elif program == "gaussian":
			aux = glob.glob("*fchk")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
	elif option == '3':
		pdb = glob.glob("*aux")		
		for i in range(len(pdb)):
			text+="3 {0} {1} {2} {3} {4} {5} 0 0 0 0 {6}\n".format(pdb[i],lh,gridn,norb,pdb[i][:-11]+".pdb",program,bandm)
	elif option == '4':
		aux = glob.glob("*aux")
		for i in sorted(aux,reverse=True):
			if aux[i][-6:-4] == "cat":
				del aux[i]
			elif aux[i][-7:-4] == "an":
				del aux[i]
		log = glob.glob("*log")
		for i in sorted(log,reverse=True):
			if aux[i][-6:-4] == "an":
				del aux[i]
			elif log[i][-7:-4] == "cat":
				del log[i]
		for i in range(len(aux)):
			text+="1 {0} {1} {2} {3} \n".format(aux[i],lh,gridn,"mopac")
		for i in range(len(log)):
			text+="1 {0} {1} {2} {3} \n".format(log[i],lh,gridn,"gamess")
		for i in range(len(aux)):
			text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		for i in range(len(log)):
			text+="2 {0} {1} {2} {3} {4} {5}\n".format(log[i],log[i][:-4]+"_cat.log",log[i][:-4]+"_an.log",lh,gridn,2,"gamess")


	f.write(text)
	f.close() 

def inp_mopac_from_all_pdbs():

	list = glob.glob("*.pdb")

	for pdb in list:
		a =  protein(pdb)
		a.pdb_parse()
		a.remove_waters()
		a.write_xyz()
		list2 = glob.glob("*.xyz")

	lmo = False
	if ( moz == "_zy" ):
		lmo = True

	for xyz in list2:
		a = mopac_inp(xyzfile=xyz,charge=0,multi=0,mozyme=lmo,inpnam=xyz[:-4]+"_"+sem+moz+".mop",method=sem,mgf=mgf)
		if not lmo:
			b = mopac_inp(xyzfile=xyz,charge=2,multi=0,mozyme=lmo,inpnam=xyz[:-4]+"_"+sem+moz+"_cat.mop",method=sem,mgf=mgf)
			c = mopac_inp(xyzfile=xyz,charge=-2,multi=0,mozyme=lmo,inpnam=xyz[:-4]+"_"+sem+moz+"_an.mop",method=sem,mgf=mgf)
			b.write_mop()
			c.write_mop()
		a.write_mop()


def script_shell_from_mop():

	list = glob.glob("*.mop")
	txt = "#!/bin/sh \n"
	f = open("mopac_shell.sh","w")
	for i in range(len(list)):
		ee= "&&"
		if i == len(list)-1:
			ee =" "
		txt += "/opt/mopac/MOPAC2016.exe {0} {1}\n".format(list[i],ee)
	f.write(txt)
	f.close()


if __name__ == "__main__":
	if ( sys.argv[1] == "-imop" ):
		sem = sys.argv[2]
		for i in range(len(sys.argv)):
			if  sys.argv[i] == "-mozyme":
				moz = "_zy"
			elif  sys.argv[i] == "-mgf":
				mgf = True
		inp_mopac_from_all_pdbs()

	elif ( sys.argv[1] == "-sh" ):
		script_shell_from_mop()

	elif ( sys.argv[1] == "-pri" ):
		LH      = "none"
		grid    = 40 
		prog    = "mopac"
		eb      = 5
		bmtd    = "BD" 
		nrb     = 100
		for i in range(len(sys.argv)):
			if ( sys.argv[i] == "-prog"):
				prog = sys.argv[i+1]
			elif ( sys.argv[i] == "-lh"):
				LH = sys.argv[i+1]
			elif ( sys.argv[i] == "-grid"):
				grid = int(sys.argv[i+1])	
		primordia_inp(option=sys.argv[2],program=prog,lh=LH,gridn=grid,eband=eb,bandm=bmtd,norb=nrb)


