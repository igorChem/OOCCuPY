#!/usr/bin/env python
# -*- coding: utf-8 -*-
# cubeclass.py


import os 
import glob 
import numpy as np

#=======================================================================
# Cube Extraction and wrtting methods 
#=======================================================================

class Cube:
	
	'''
	Class to read/write and modify .cube files from QM runs
	'''
          
	def __init__(self,Typ="MO",Prog="Gaussian",grid = 40):
		
		'''
		Initialization method for the Elec_Cube class
		'''
		
		self.name = ''
		self.header = ''
		self.Prog = Prog
		self.Typ  = Typ
		self.natoms = []
		self.threeLI = ''
		self.grid  = grid
		self.atoms  = ''
		self.scalar3d = []
		self.MO = []
		self.origin = []   
		self.gridx = []
		self.gridy = []
		self.gridz = []
        
	def Read_Elec_Cube (self,filename):
		
		'''
		Method to parse electronic data cube and fill the atributes
		of Elec_cube class 
		'''
		
		cubefile = open(filename,'r')	
		cubefile2 =	open(filename,'r')			
	
		for i in range(2):
			self.header += cubefile.readline()
			
		for i in range(1,2):	
			tmp = cubefile.readline().split()
			self.natoms = (int(tmp[0]))
			for j in range(3): 
				self.origin.append(tmp[j+1]) 	
			
		NA = int(abs(float(self.natoms)) + 6)


		if 	self.Typ == 'MO':
			self.MO = self.header.split()[7]
			NA = int(NA + 1)
			print(NA)
			
		for i in range(2,5):
			self.threeLI += cubefile.readline()
		
		for i in range(6):
			tmp = cubefile2.readline()
			if i == 3:
				self.gridx=int(tmp.split()[0])
			elif i == 4:
				 self.gridy=int(tmp.split()[0])
			elif i == 5:
				self.gridz=int(tmp.split()[0])
		
		    
		num = abs(self.natoms) + 5
		for i in range(5,num):
			self.atoms += cubefile.readline()
			
		cubefile2 = open(filename,'r')		
		str = ' '.join(cubefile2.readlines()[NA:])
		self.scalar3d = np.fromstring(str,sep = ' ')
		self.scalar3d.shape =(self.gridx, self.gridy, self.gridz)


	def write_cubeElec (self,filename2):
		
		'''
		Method to write the cube data in gaussian like style		
		'''
		
		text_to_write =''
		text_to_write += self.header 
		text_to_write += '   ' + str(self.natoms) +'   '+ str(self.origin[0]) + '  ' + str(self.origin[1]) + '  '+ str(self.origin[2]) + '\n'
		text_to_write += self.threeLI
		text_to_write += self.atoms	
		
		if 	self.Typ == 'MO':
			text_to_write += '    1   ' + str(self.MO) + ' \n'

		for i in range(self.gridx):
			for j in range(self.gridy):
				for k in range(self.gridz):
					text_to_write += '{0:.5e}'.format(self.scalar3d[i][j][k]) + '   '
					if (k % 6 == 5):
						text_to_write += ' \n'
				text_to_write +=  ' \n'
		
			
		file_write = open(filename2+'.cube','w')
		file_write.write(text_to_write)
        
	def __sub__ (self,other):
		
		'''
		Overload method to operte subtraction of scalar data from different
		Elec_Cube object
		'''		
		
		z = Cube()
		
		z.name = self.name + 'mod'
		z.header = self.header
		z.Prog = self.Prog
		z.Typ = self.Typ 
		z.natoms = self.natoms
		z.threeLI = self.threeLI
		z.grid = self.grid 
		z.atoms = self.atoms  
		z.scalar3d = self.scalar3d
		z.MO = self.MO
		z.origin = self.origin
		 
		z.scalar3d = np.absolute(self.scalar3d - other.scalar3d)

		return(z)
		

				
		
		
		
'''
a = Cube(Typ="Elecdensity")
a.Read_Elec_Cube('rta_Fukui0.cube')


os.system('sudo cp ~/Dropbox/RDscripts/orcamodule/cubeclass.py /usr/lib/python2.7')
'''	
		
		
		
		
		
		
		
		
		
		
		
