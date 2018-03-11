#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mopac_module.py

#=======================================================================

#load modules

import os
from xyz_class import*

#========================================================================

class mopac_inp:

	def __init__(self  ,
				xyzfile,
				charge ,
				multi  , 
				inpnam ,
				method):					
				
		self.name = xyzfile
		self.inpnam = inpnam
		self.charge = charge 
		self.multi = multi 
		self.mult = "Singlet"
		self.solvent = True
		self.hamilt = method
		self.xyz = None

		if not self.multi == 1:
			self.mult = "Doublet"

		self.xyz = xyz_parser(self.name)
		self.xyz.parse_xyz()

	def write_mop(self):

		mop_inp = open(self.inpnam,'w')

		mop_text = '' 
		mop_text += '{0} 1SCF  PL T=1D TIMES charge={1} {2} AUX VECTORS ALLVECS \n\n\n'.format(self.hamilt,self.charge,self.mult)

		for i in range(self.xyz.Natoms):
			mop_text +="{0}  {1}  1 {2}  1 {3} \n".format(self.xyz.AtomLabels[i],self.xyz.xCoord[i],self.xyz.yCoord[i],self.xyz.zCoord[i])

		mop_inp.write(mop_text)
		mop_inp.close()

