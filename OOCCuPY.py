#/usr/bin/env python
# -*- coding: utf-8 -*-
#OOCCuPY.py

from pdb_class import *
from mopac_module import *
from primordia import *
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
			elif ( sys.argv[i] == "-norb"):
				nrb = int(sys.argv[i+1])	
		primordia_inp(option=sys.argv[2],program=prog,lh=LH,gridn=grid,eband=eb,bandm=bmtd,norb=nrb)

	elif ( sys.argv[1] == "-prd"):
		a = pair_RD()
		a.write()
		a.r_scripts()


