#!/usr/bin/env python
# -*- coding: utf-8 -*-
# xyz_class.py


#load modules

import os 
import glob 

	
#=======================================================================

#=======================================================================
# A class to storage xyz molecular information
#=======================================================================

atomnumber = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':36,'BR':36,'CL':17}

class xyz_parser:
	
	'''
	Class to read/write and storage xyz molecula coordinates file information
	'''
	def __init__(self,name):
		self.name = name 
		self.AtomLabels = []
		self.AtomicN = []
		self.xCoord = []
		self.yCoord = []
		self.zCoord = []
		self.Natoms = []
		self.origin = [0,0,0]
		self.xvec = [0,0,0]
		self.yvec = [0,0,0]
		self.zvec = [0,0,0]
		self.ElNumb = []
		self.elsden_text =''
		
	
	def parse_xyz(self):
		
		logfile = self.name
		log = open(logfile,'r')	
		
		
		for line in log:
			line2 = line.split()
			if len(line2) == 4:
				self.AtomLabels.append(line2[0])
				self.xCoord.append(float(line2[1]))
				self.yCoord.append(float(line2[2]))
				self.zCoord.append(float(line2[3]))
				
		self.Natoms = len(self.AtomLabels)		
					
					
	def get_atomnumber(self):
		
		for i in range(len(self.AtomLabels)):
			self.AtomicN.append(atomnumber[self.AtomLabels[i]])
	
	def parse_log(self):
				
		logfile = self.name
		
		phrase1 = '***** EQUILIBRIUM GEOMETRY LOCATED *****'
		phrase2 = 'INTERNUCLEAR DISTANCES (ANGS.)' 
		
		xyz_init = []
		xyz_fin  = []
		
		with open(logfile,'r') as text:
			for (i, line) in enumerate(text):
				if phrase1 in line:
					xyz_init.append(i)
				elif phrase2 in line:
					xyz_fin.append(i)
					
	
		with open(logfile,'r') as text:
			for(i,line) in enumerate(text):
				if i >= (xyz_init[0]+4) and i <= (xyz_fin[-1]-2):
					line2 = line.split()
					if len(line2) == 5:
						self.AtomLabels.append(line2[0])
						self.AtomicN.append(line2[1])
						self.xCoord.append(line2[2])
						self.yCoord.append(line2[3])
						self.zCoord.append(line2[4])	
	
	def mol2(self):	
						
		logfile = self.name
		log = open(logfile,'r')	
		
		
		for line in log:
			line2 = line.split()
			if len(line2) == 9:
				self.AtomLabels.append(line2[1])
				self.xCoord.append(line2[2])
				self.yCoord.append(line2[3])
				self.zCoord.append(line2[4])
				
		self.Natoms = len(self.AtomLabels)						
				
	def write_text(self,filename = '',writefile = False):
		
		if writefile == True:
			text_to_return = str(len(self.AtomLabels)) + '\n \n'
			for i in range(len(self.AtomLabels)):	
				text_to_return += str(self.AtomLabels[i]) +' '+ str(self.xCoord[i]) + ' '+str(self.yCoord[i]) + ' '+str(self.zCoord[i]) + '\n'
		else:
			text_to_return = ''		
			for i in range(len(self.AtomLabels)):
				text_to_return += str(self.AtomLabels[i]) +' '+ str(self.AtomicN[i]) +' '+ str(self.xCoord[i]) + ' '+str(self.yCoord[i]) + ' '+str(self.zCoord[i]) + '\n'
		
		
		if writefile == True:
			fileout = open(filename,'w')
			fileout.write(text_to_return)
			fileout.close()
			
		return(text_to_return)

		
	def get_origin(self): 
		
		self.origin[0] = min(self.xCoord)-3
		self.origin[1] = min(self.yCoord)-3
		self.origin[2] = min(self.zCoord)-3
		
		self.xvec[0] = max(self.xCoord) + 3
		self.xvec[1] = self.origin[1]
		self.xvec[2] = self.origin[2]
		
		self.yvec[0] = self.origin[0]
		self.yvec[1] = max(self.yCoord) + 3
		self.yvec[2] = self.origin[2]
		
		self.zvec[0] = self.origin[0]
		self.zvec[1] = self.origin[1]
		self.zvec[2] = max(self.zCoord) + 3
		
		text_to_print = ''
		text_to_print += ' $grid origin(1) = ' +str(self.origin[0])+','+str(self.origin[1])+','+str(self.origin[2])+' $end\n'
		text_to_print += ' $grid xvec(1)= '+str(self.xvec[0])+','+str(self.xvec[1])+','+str(self.xvec[2])+' $end\n'
		text_to_print += ' $grid yvec(1)= '+str(self.yvec[0])+','+str(self.yvec[1])+','+str(self.yvec[2])+' $end\n'
		text_to_print += ' $grid zvec(1)= '+str(self.zvec[0])+','+str(self.zvec[1])+','+str(self.zvec[2])+' $end \n'
		
		self.elsden_text = text_to_print
		
		return(self.origin,self.xvec,self.yvec,self.zvec)
		
	def get_elecNumber(self): 
		
		self.ElNumb = sum(self.AtomicN)
		
		print(self.ElNumb)
		return(self.ElNumb)
	
	def read_mop(self):		
		mopfile = self.name
		mop = open(mopfile,'r')	
		i=0
		for line in mop:
			if i>3:
				line2 = line.split()
				if len(line2) == 7:
					self.AtomLabels.append(line2[0])
					self.xCoord.append(line2[1])
					self.yCoord.append(line2[3])
					self.zCoord.append(line2[5])
			i+=1
					
					
