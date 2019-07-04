#!/usr/bin/env python
# -*- coding: utf-8 -*-
# script_input.py

from pdb_class import *
from mopac_module import *
import glob

list = glob.glob("*.pdb")

for pdb in list:
	a =  protein(pdb)
	a.pdb_parse(pdb)
	a.write_xyz()

list2 = glob.glob("*.xyz")

for xyz in list2:
	a = mopac_inp(xyzfile=xyz,charge=0,multi=0,inpnam=xyz[:-4]+".mop",method="PM7")
	a.write_mop()


def inp(pdb):
	a = protein(pdb)
	a.pdb_parse()
	a.write_xyz()
	a = mopac_inp(xyzfile=pdb[:-4]+".xyz",charge=0,multi=0,inpnam=pdb[:-4]+".mop",method="PM7")
	a.write_mop()
 
