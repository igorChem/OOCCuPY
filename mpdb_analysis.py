#/usr/bin/env python
# -*- coding: utf-8 -*-
#mpdb_analuysis.py

import os, math, sys

class m_pdb:
	
	def __init__(self,name):
		self.file   = name 
		self.cx1 = [] 
		self.cx2 = [] 
		self.cx3 = [] 
		self.cx4 = [] 
		self.cy1 = [] 
		self.cy2 = [] 
		self.cy3 = [] 
		self.cy4 = [] 
		self.cz1 = [] 
		self.cz2 = [] 
		self.cz3 = [] 
		self.cz4 = [] 
		self.r1  = []
		self.r2  = []
		self.npairs = 1

	def set_pairs(self,index):
		
		fl = open(self.file,'r')
		for line in fl:
			ln= line.split()
			if len(ln) > 1:
				if ln[0] == "ATOM" or ln[0] == "HETATM":
					for i in range(len(index)):
						if ln[1] == index[i]:
							if i == 0:
								self.cx1.append(float(ln[6]))
								self.cy1.append(float(ln[7]))
								self.cz1.append(float(ln[8]))
							elif i == 1:
								self.cx2.append(float(ln[6]))
								self.cy2.append(float(ln[7]))
								self.cz2.append(float(ln[8]))
							elif i == 2:
								self.cx3.append(float(ln[6]))
								self.cy3.append(float(ln[7]))
								self.cz3.append(float(ln[8]))
							elif i == 3:
								self.cx4.append(float(ln[6]))
								self.cy4.append(float(ln[7]))
								self.cz4.append(float(ln[8]))
		if len(self.cx3) > 1:
			self.npairs = 2
								
	def calc_distances(self):
		
		for i in range(len(self.cx1)):
			xx = (self.cx1[i] - self.cx2[i])**2
			yy = (self.cy1[i] - self.cy2[i])**2
			zz = (self.cz1[i] - self.cz2[i])**2
			self.r1.append(math.sqrt(xx+yy+zz))
			if self.npairs == 2:
				xx = (self.cx3[i] - self.cx4[i])**2
				yy = (self.cy3[i] - self.cy4[i])**2
				zz = (self.cz3[i] - self.cz4[i])**2
				self.r2.append(math.sqrt(xx+yy+zz))
		
		dist = open("dist_analysis",'w')
		if self.npairs == 2:
			dist_text = "r1 r2\n"
			for j in range(len(self.r1)):
				dist_text +="{} {}\n".format(self.r1[j],self.r2[j])
		else:
			dist_text = "r1 r2\n"
			for j in range(len(self.r1)):
				dist_text +="{}\n".format(self.r1[j])
		dist.write(dist_text)
		dist.close()		
		
		r_s = open("r_script.R",'w')
		
		r_s_text =  "library(hexbin)\n"
		r_s_text += "library(RColorBrewer)\n"
		r_s_text += "rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))\n"
		r_s_text += "r <-rf(32)\n"
		r_s_text += "dat = read.table('dist_analysis',header=T)\n"
		r_s_text +=	"h<-hexbin(dat)\n"
		r_s_text += "png('dist2D.png',units='in',res=400,width=5,height=4)\n"
		r_s_text +=	"plot(h,colramp=rf,xlab='Pair 1',ylab='Pair 2')\n"
		r_s_text += "dev.off()\n"
		r_s_text += "summary(dat)\n"
		
		r_s.write(r_s_text)
	
	def find_frame(self,r1,r2):
		a1 = []
		
		for i in range(len(self.r1)):
			if self.r1[i] >= (r1-0.02) and self.r1[i] <= (r1+0.02) and self.r2[i] >= (r2-0.02) and self.r2[i] <= (r2+0.02):
				a1.append(i)
		print(a1)		
				
				
	
				
					
	
	
	
