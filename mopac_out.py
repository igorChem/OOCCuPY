#!/usr/bin/env python
# -*- coding: utf-8 -*-
# mopac_out.py

#class to parse mopac output info 


class mopac_out:

	def __init__(self,outfile, method="RM1"):
		self.name = outfile
		self.aux_name = outfile[:-4] + "aux"
		self.energy = 0
		self.heat = 0
		self.energy = 0
		self.atoms = []
		self.homo_en = 0
		self.lumo_en = 0
		self.method = method

	def parse_out(self):

		out_file = open("outfile",'r')

		for line in out_file:
			line2 = line.split()
			if len(line2) == 9:
				if line2[1] == "HEAT" amd line2[3] == "FORMATION":
					self.heat = line2[7]
			elif len(line2) == 5:
				if line2[0] == "TOTAL" and line2[1] == "ENERGY":
					self.energy = line2[3]
			elif len(line2) == 7:
				elif line2[0] == "HOMO" and line2[1] == "LUMO":
					self.homo_en = line2[5]
					self.lumo_en = line2[6]




