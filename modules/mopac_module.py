#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mopac_module.py

#=======================================================================

#load modules

import os
from xyz_class import*

#========================================================================

class mopac_inp:

	def __init__(self,xyzfile):
		self.name = xyzfile
		self.charge = 0
		self.multi = 1
		self.mult = "Singlet"
		self.solvent = True
		self.hamilt = "RM1"
		self.xyz = None

		if not self.multi == 1:
			self.mult = "Duplet"

		self.xyz = xyz_parser(self.name)
		self.xyz.parse_xyz()

	def write_mop(self):

		mop_inp = open(self.name[:-4]+".mop",'w')

		mop_text = '' 
		mop_text += '{0} 1SCF  PL T=1D TIMES MOZYME charge={1} {2} AUX \n\n'.format(self.hamilt,self.charge,self.mult)

		for i in range(self.xyz.Natoms):
			mop_text +="{0}  {1}  1 {2}  1 {3} \n".format(self.xyz.AtomLabels[i],self.xyz.xCoord[i],self.xyz.yCoord[i],self.xyz.zCoord[i])

		mop_inp.write(mop_text)
		mop_inp.close()


test = mopac_inp("1l2y_min.xyz")
test.write_mop()
