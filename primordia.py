#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pdb_Reader.py

import os,glob


def ls_gen():
	a = []
	sc = ",";
	while True:
		sc = input()
		if sc == ".":
			return(a)
		else:
			a.append(sc)


class primordia_inp:
	pdb = glob.glob("*pdb")
	f = open("input_pri",'w')
	text ="eband 5\n"
	for i in range(len(pdb)):
		text+="3 {0} potential_fukui 0 100 {1} mopac 0 0 0 0 EW elecdens mep\n".format(pdb[i][:-4]+".aux",pdb[i])
			
	f.write(text)
	f.close()

class primordia:
	def __init__(self,comp,pro,listr):
		self.name = ""
		self.eas_c    = []
		self.nas_c   = []
		self.ras_c   = []
		self.dual_c  = []
		self.pot_c   = []
		self.eas_p    = []
		self.nas_p   = []
		self.ras_p   = []
		self.dual_p  = []
		self.pot_p   = []
		self.eas_ch   = []
		self.nas_ch  = []
		self.ras_ch  = []
		self.dual_ch = []
		self.pot_ch  = []
	
		f2 = open(pro,'r')
		f1 = open(comp,'r')
		c = 0
		i = 0
		self.name = pro[:4]
		for line in f1:
			if len(listr) >=c: 
				if i == listr[c]:
					line2 = line.split()
					self.eas_c.append(float(line2[1]))
					self.nas_c.append(float(line2[2]))
					self.ras_c.append(float(line2[3]))
					self.dual_c.append(float(line2[4]))
					self.pot_c.append(float(line2[6]))
					if len(listr) == c+1:
						break
					else:
						c +=1
			i +=1
			
		c = 0
		i=0		
		for line in f2:
			if len(listr) >=c: 
				if i == listr[c]:
					line2 = line.split()
					self.eas_p.append(float(line2[1]))
					self.nas_p.append(float(line2[2]))
					self.ras_p.append(float(line2[3]))
					self.dual_p.append(float(line2[4]))
					self.pot_p.append(float(line2[6]))
					if len(listr) == c+1:
						break
					else:
						c +=1
			i +=1	
		
		f1.close()
		f2.close()
		
		for i in range(len(self.eas_p)):
			self.eas_ch.append(self.eas_c[i]-self.eas_p[i])
			self.nas_ch.append(self.nas_c[i]-self.nas_p[i])
			self.ras_ch.append(self.ras_c[i]-self.ras_p[i])
			self.dual_ch.append(self.dual_c[i]-self.dual_p[i])
			self.pot_ch.append(self.pot_c[i]-self.pot_p[i])
		
		self.eas_c.append(sum(self.eas_c)) 
		self.nas_c.append(sum(self.nas_c)) 
		self.ras_c.append(sum(self.ras_c)) 
		self.dual_c.append(sum(self.dual_c)) 
		self.pot_c.append(sum(self.pot_c)) 
		self.eas_p.append(sum(self.eas_p))  
		self.nas_p.append(sum(self.nas_p))     
		self.ras_p.append(sum(self.ras_p))     
		self.dual_p.append(sum(self.dual_p))   
		self.pot_p.append(sum(self.pot_p))     
		self.eas_ch.append(sum(self.eas_ch))    
		self.nas_ch.append(sum(self.nas_ch))   
		self.ras_ch.append(sum(self.ras_ch))   
		self.dual_ch.append(sum(self.dual_ch)) 
		self.pot_ch.append(sum(self.pot_ch))
		
			
			
		f = open(self.name+".pri",'w')
		tesxt = ""
		for i in range(len(self.eas_ch)):
			tesxt += "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(self.eas_c[i],self.nas_c[i],self.ras_c[i],self.dual_c[i],self.pot_c[i],self.eas_p[i],self.nas_p[i],self.ras_p[i],self.dual_p[i],self.pot_p[i],self.eas_ch[i],self.nas_ch[i],self.ras_ch[i],self.dual_ch[i],self.pot_ch[i])
		f.write(tesxt)
		f.close()

    
     
