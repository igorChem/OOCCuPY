#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pdb_Reader.py

import os,glob
import datetime
import numpy as np

x = datetime.datetime.now()



class pair_RD:
	def __init__(self,mode="1d"):
		self.prds = glob.glob("*.prd")
		self.pairs = 1
		self.eas_a1 = []
		self.eas_a2 = []
		self.eas_a3 = []
		self.eas_a4 = []
		self.chg_a1 = []
		self.chg_a2 = []
		self.chg_a3 = []
		self.chg_a4 = []		
		self.nas_a1 = []
		self.nas_a2 = []
		self.nas_a3 = []
		self.nas_a4 = []
		self.hardness_a1 = []
		self.hardness_a2 = []
		self.hardness_a3 = []
		self.hardness_a4 = []
		self.CT_p1 = []
		self.CT_p2 = []
		self.SPI_p1 = []
		self.SPI_p2 = []
		self.HPI_p1 = []
		self.HPI_p2 = []
		self.HOF = []
		self.Elec_en = []
		self.hardness = []
		self.ECP = []
		self.electrophilicity = []
		self.gstep  = []
		self.gstep2  = []
		self.lstep  = []
		self.lstep2 = []
		self.CT_p1_p2 = []
		self.SPI_p1_p2 = []
		self.HPI_p1_p2 = []
		self.mode = mode 
		self.dist1   = []
		self.dist2   = []
		self.hamilt  = []
		self.lhamilt = []
		

		globaln = sorted(glob.glob("*global"))
		fgl = open(globaln[0],'r')
	
		#read global descriptors 
		#------------------------------------------------------------------#
		i=0
		
		for line in fgl:
			if i> 1:
				line2 = line.split()
				if mode == "2d":
					if len(line2[0]) == 16:
						self.gstep.append(int(line2[0][7:9]))
						self.gstep2.append(int(line2[0][10:12]))
						self.hamilt.append(line2[0][13:16])
					elif len(line2[0]) == 15:
						g1 = 0
						g2 = 0
						try:
							g1 = int(line2[0][7:8])
							g2 = int(line2[0][9:11])							
						except:
							try:
								g1 = int(line2[0][7:9])
								g2 = int(line2[0][10:11])
							except:
								print(line2)
								print(line2[0][7:9])
								print(line2[0][10:11])
						
						self.gstep.append(g1)
						self.gstep2.append(g2)
						self.hamilt.append(line2[0][12:15])
					elif len(line2[0]) == 14:
						self.gstep.append(int(line2[0][7:8]))
						self.gstep2.append(int(line2[0][9:10]))
						self.hamilt.append(line2[0][11:14])
				else:
					line2 = line.split()
					if len(line2[0]) > 10:
						#self.gstep.append(int(line2[0][6:-7]))
						try:
							self.gstep.append(int(line2[0][11:-9]))
						except:
							self.gstep.append(int(line2[0][11:-10]))
						self.gstep2.append(0)
						#self.hamilt.append(line2[0][-6:-3])
						self.hamilt.append("am1")
					else:
						self.gstep.append(int(line2[0][6:]))
						self.hamilt.append("pm7")
				self.HOF.append(float(line2[2]))
				self.Elec_en.append(float(line2[1]))
				self.hardness.append(float(line2[6]))
				self.ECP.append(float(line2[5]))
				self.electrophilicity.append(float(line2[7]))
			i +=1
		fgl.close()
		
		for j in range(len(self.prds)):
			prd = open(self.prds[j],'r')
			if len(self.prds[j]) > 16:
				if mode == "1d":
					#self.lstep.append(int(self.prds[j][6:-15]))
					try:
						self.lstep.append(int(self.prds[j][9:-17]))
					except:
						self.lstep.append(int(self.prds[j][11:-18]))
					self.lstep2.append(0)
					#self.lhamilt.append(self.prds[j][-14:-11])
					self.lhamilt.append("am1")
				elif mode == "2d":
					if len(self.prds[j]) == 24:
						self.lstep.append(int(self.prds[j][7:9]))
						self.lstep2.append(int(self.prds[j][10:12]))
						self.lhamilt.append(self.prds[j][13:16])
					elif len(self.prds[j]) == 23:
						g1 = 0
						g2 = 0
						try:
							g1 = int(self.prds[j][7:8])
							g2 = int(self.prds[j][9:11])
							self.lstep.append(g1)
							self.lstep2.append(g2) 
						except:							
							g1 = int(self.prds[j][7:9])
							g2 = int(self.prds[j][10:11])
							self.lstep.append(g1)
							self.lstep2.append(g2)
						self.lhamilt.append(self.prds[j][12:15])
					elif len(self.prds[j]) == 22:
						self.lstep.append(int(self.prds[j][7:8]))
						self.lstep2.append(int(self.prds[j][9:10]))
						self.lhamilt.append(self.prds[j][11:14])
			else:
				self.lstep.append(int(self.prds[j][6:-8]))
				self.lhamilt.append("PM7")
			i = 0
			for line in prd:
				line2 = line.split()
				if i == 0:
					if len(line2) == 7:
						self.pairs = 2
				if i > 0:
					if line2[0] == "dist":
						self.dist1.append(float(line2[1]))
						if self.pairs == 2:
							self.dist2.append(float(line2[2]))
					elif line2[0] == "EAS":
						self.eas_a1.append(float(line2[1]))
						self.eas_a2.append(float(line2[2]))
						if self.pairs == 2:
							self.eas_a3.append(float(line2[3]))
							self.eas_a4.append(float(line2[4]))
					elif line2[0] == "NAS":
						self.nas_a1.append(float(line2[1]))
						self.nas_a2.append(float(line2[2]))
						if self.pairs == 2:
							self.nas_a3.append(float(line2[3]))
							self.nas_a4.append(float(line2[4]))
					elif line2[0] == "Hardness":
						self.hardness_a1.append(float(line2[1]))
						self.hardness_a2.append(float(line2[2]))
						if self.pairs == 2:
							self.hardness_a3.append(float(line2[3]))
							self.hardness_a4.append(float(line2[4]))
					elif line2[0] == "Charge":
						self.chg_a1.append(float(line2[1]))
						self.chg_a2.append(float(line2[2]))
						if self.pairs == 2:
							self.chg_a3.append(float(line2[3]))
							self.chg_a4.append(float(line2[4]))
					elif line2[0] == "CT":
						self.CT_p1.append(float(line2[5]))
						if self.pairs == 2:
							self.CT_p2.append(float(line2[6]))
							self.CT_p1_p2.append(self.CT_p1[j] + self.CT_p2[j])
					elif line2[0] == "SPI":
						self.SPI_p1.append(float(line2[5]))
						if self.pairs == 2:
							self.SPI_p2.append(float(line2[6]))
							self.SPI_p1_p2.append(self.SPI_p1[j] + self.SPI_p2[j])
					elif line2[0] == "HPI":
						self.HPI_p1.append(float(line2[5]))
						if self.pairs == 2:
							self.HPI_p2.append(float(line2[6]))
							self.HPI_p1_p2.append(self.HPI_p1[j] + self.HPI_p2[j])

				i+=1
			prd.close()
		
	def write(self):
		
		'''
		shof = min(self.HOF)
		selec = min(self.Elec_en)
		for i in range(len(self.gstep)):
			self.HOF[i] = self.HOF[i] - shof
			self.Elec_\en[i] = self.Elec_en[i] - selec
		'''
		
		half_dist1 = 0
		half_dist2 = 0
		for i in range(len(self.lstep)):
			if len(self.lstep) %2 == 0:
				if self.lstep[i] == len(self.lstep)/2:
					half_dist1 = self.dist1[i]
					if self.pairs == 2:
						half_dist2 = self.dist2[i]
			elif not len(self.lstep) %2 == 0:
				if self.lstep[i] == (len(self.lstep)-1)/2:
					half_dist1 = self.dist1[i]
					if self.pairs == 2:
						half_dist2 = self.dist2[i]
					
		for i in range(len(self.dist1)):
			if self.dist1[i] > half_dist1:
				self.dist1[i] = half_dist1 - self.dist1[i]
			elif self.dist1[i] < half_dist1:
				self.dist1[i] =  half_dist1 - self.dist1[i]
			else:
				self.dist1[i] = 0
			if self.pairs == 2:
				if self.dist2[i] > half_dist2:
					self.dist2[i] = half_dist2 - self.dist2[i]
				elif self.dist2[i] < half_dist2:
					self.dist2[i] = self.dist2[i] - half_dist2
				else:
					self.dist1[i] = 0	
					
		
		fgr_text = "n HOF Energy Hardness ECP Electrophilicity method\n"
		fgr = open("global_resume"+str(x.hour)+"_"+str(x.minute),'w')
		print(len(self.gstep),len(self.gstep2),len(self.HOF),len(self.Elec_en),len(self.hardness),len(self.ECP),len(self.electrophilicity),len(self.hamilt))
		for i in range(len(self.gstep)):
			fgr_text += "{} {} {} {} {} {} {} {}\n".format(self.gstep[i],self.HOF[i],self.Elec_en[i],self.hardness[i],self.ECP[i],self.electrophilicity[i],self.hamilt[i],self.gstep2[i])
		fgr.write(fgr_text)
		fgr.close()
		
		flr_text = ""
		flr = open("local_resume"+str(x.hour)+"_"+str(x.minute),'w')
		if self.pairs == 1:
			flr_text = "n eas_a1 nas_a1 hardness_a1 chg_a1 chg_a2 eas_a2 nas_a2 hardness_a2 CT SPI HPI method\n"
			for i in range(len(self.lstep)):
				flr_text += "{} {} {} {} ".format(self.lstep[i],self.eas_a1[i],self.nas_a1[i],self.hardness_a1[i])
				flr_text += "{} {} ".format(self.chg_a1[i],self.chg_a2[i])
				flr_text += "{} {} {} {} {} {} {}\n".format(self.eas_a2[i],self.nas_a2[i],self.hardness_a2[i],self.CT_p1[i],self.SPI_p1[i],self.HPI_p1[i],self.lhamilt[i])
		elif self.pairs == 2:
			flr_text = "n eas_a1 nas_a1 hardness_a1 eas_a2 nas_a2 hardness_a2 eas_a3 nas_a3 hardness_a3 eas_a4 nas_a4 hardness_a4 chg_a1 chg_a2 chg_a3 chg_a4 CT1 SPI1 HPI1 CT2 SPI2 HPI2 method CT12 SPI12 HPI12 m\n"
			for i in range(len(self.lstep)):
				flr_text += "{} {} {} {} ".format(self.lstep[i],self.eas_a1[i],self.nas_a1[i],self.hardness_a1[i])
				flr_text += "{} {} {} ".format(self.eas_a2[i],self.nas_a2[i],self.hardness_a2[i])
				flr_text += "{} {} {} ".format(self.eas_a3[i],self.nas_a3[i],self.hardness_a3[i])
				flr_text += "{} {} {} ".format(self.eas_a4[i],self.nas_a4[i],self.hardness_a4[i])
				flr_text += "{} {} {} {} ".format(self.chg_a1[i],self.chg_a2[i],self.chg_a3[i],self.chg_a4[i])
				flr_text += "{} {} {} {} {} {} {} ".format(self.CT_p1[i],self.SPI_p1[i],self.HPI_p1[i],self.CT_p2[i],self.SPI_p2[i],self.HPI_p2[i],self.lhamilt[i])
				flr_text += "{} {} {} {}\n".format(self.CT_p1_p2[i],self.SPI_p1_p2[i],self.HPI_p1_p2[i],self.lstep2[i])
		flr.write(flr_text)
		flr.close()

	def r_scripts(self):
		
		r_text = "#r script generated for OOCCuPY for PRIMoRDiA software\n"
		r_text += "gldat = read.table('global_resume_data',header=T)\n"
		r_text += "lldat = read.table('local_resume_data',header=T)\n"
		r_text += "attach(gldat)\n"
		r_text += "tiff('hof',units='in',width=4,height=4,res=400)\n"
		r_text += "plot(HOF~n,data=gldat,type='p',col='black',xlab='Step',ylab='Heat of Formation(kcalmol)')\n"
		r_text += "dev.off()\n"
		r_text += "tiff('hardness',units='in',width=4,height=4,res=400)\n"
		r_text += "plot(Hardness~n,data=gldat,type='p',col='green',xlab='Step',ylab='Hardness(eV)')\n"
		r_text += "dev.off()\n"
		r_text += "detach(gldat)\n"
		r_text += "attach(lldat)\n"
		
		if self.pairs == 1:			
			r_text += "tiff('CT',units='in',width=4,height=4,res=400)\n"
			r_text += "plot(CT~n,data=lldat,type='p',col='blue',xlab='Step',ylab='Charge Transfer')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('SPI',units='in',width=4,height=4,res=400)\n"
			r_text += "plot(SPI~n,data=lldat,type='p',col='orange',xlab='Step',ylab='Softness Pair Interaction')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('HPI',units='in',width=4,height=4,res=400)\n"
			r_text += "plot(HPI~n,data=lldat,type='p',col='red',xlab='Step',ylab='Hardness Pair Interaction')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('atom1',units='in',width=6,height=6,res=400)\n"
			r_text += "par(mfrow=c(2,2))\n"
			r_text += "plot(eas_a1~n,data=lldat,type='p',col='blue',xlab='Step',ylab='EAS')\n"
			r_text += "plot(nas_a1~n,data=lldat,type='p',col='red',xlab='Step',ylab='NAS')\n"
			r_text += "plot(hardness_a1~n,data=lldat,type='p',col='green',xlab='Step',ylab='Local Hardness')\n"
			r_text += "plot(chg_a1~n,data=lldat,type='p',col='black',xlab='Step',ylab='Partial Charge')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('atom2',units='in',width=6,height=6,res=400)\n"
			r_text += "par(mfrow=c(2,2))\n"
			r_text += "plot(eas_a2~n,data=lldat,type='p',col='blue',xlab='Step',ylab='EAS')\n"
			r_text += "plot(nas_a2~n,data=lldat,type='p',col='red',xlab='Step',ylab='NAS')\n"
			r_text += "plot(hardness_a2~n,data=lldat,type='p',col='green',xlab='Step',ylab='Local Hardness')\n"
			r_text += "plot(chg_a2~n,data=lldat,type='p',col='black',xlab='Step',ylab='Partial Charge')\n"
			r_text += "dev.off()\n"
		elif self.pairs == 2:			
			r_text += "tiff('CT_p1',units='in',width=4,height=5,res=400)\n"
			r_text += "plot(CT1~n,data=lldat,type='p',col='blue',xlab='Step',ylab='Charge Transfer')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('SPI1',units='in',width=4,height=5,res=400)\n"
			r_text += "plot(SPI1~n,data=lldat,type='p',col='orange',xlab='Step',ylab='Softness Pair Interaction')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('HPI_p1',units='in',width=4,height=5,res=400)\n"
			r_text += "plot(HPI1~n,data=lldat,type='p',col='green',xlab='Step',ylab='Hardness Pair Interaction')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('CT_p2',units='in',width=4,height=5,res=400)\n"
			r_text += "plot(CT2~n,data=lldat,type='p',col='blue',xlab='Step',ylab='Charge Transfer')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('SPI2',units='in',width=4,height=5,res=400)\n"
			r_text += "plot(SPI2~n,data=lldat,type='p',col='orange',xlab='Step',ylab='Softness Pair Interaction')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('HPI_p2',units='in',width=4,height=5,res=400)\n"
			r_text += "plot(HPI2~n,data=lldat,type='p',col='green',xlab='Step',ylab='Hardness Pair Interaction')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('atom1',units='in',width=4,height=5,res=400)\n"
			r_text += "par(mfrow=c(2,2))\n"
			r_text += "plot(eas_a1~n,data=lldat,type='p',col='blue',xlab='Step',ylab='EAS')\n"
			r_text += "plot(nas_a1~n,data=lldat,type='p',col='red',xlab='Step',ylab='NAS')\n"
			r_text += "plot(hardness_a1~n,data=lldat,type='p',col='green',xlab='Step',ylab='Local Hardness')\n"
			r_text += "plot(chg_a1~n,data=lldat,type='p',col='black',xlab='Step',ylab='Partial Charge')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('atom2',units='in',width=4,height=5,res=400)\n"
			r_text += "par(mfrow=c(2,2))\n"
			r_text += "plot(eas_a2~n,data=lldat,type='p',col='blue',xlab='Step',ylab='EAS')\n"
			r_text += "plot(nas_a2~n,data=lldat,type='p',col='red',xlab='Step',ylab='EAS')\n"
			r_text += "plot(hardness_a2~n,data=lldat,type='p',col='green',xlab='Step',ylab='Local Hardness')\n"
			r_text += "plot(chg_a2~n,data=lldat,type='p',col='black',xlab='Step',ylab='Partial Charge')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('atom3',units='in',width=4,height=5,res=400)\n"
			r_text += "par(mfrow=c(2,2))\n"
			r_text += "plot(eas_a3~n,data=lldat,type='p',col='blue',xlab='Step',ylab='EAS')\n"
			r_text += "plot(nas_a3~n,data=lldat,type='p',col='red',xlab='Step',ylab='NAS')\n"
			r_text += "plot(hardness_a3~n,data=lldat,type='p',col='green',xlab='Step',ylab='Local Hardness')\n"
			r_text += "plot(chg_a3~n,data=lldat,type='p',col='black',xlab='Step',ylab='Partial Charge')\n"
			r_text += "dev.off()\n"
			r_text += "tiff('atom4',units='in',width=4,height=5,res=400)\n"
			r_text += "par(mfrow=c(2,2))\n"
			r_text += "plot(eas_a4~n,data=lldat,type='p',col='blue',xlab='Step',ylab='EAS')\n"
			r_text += "plot(nas_a4~n,data=lldat,type='p',col='red',xlab='Step',ylab='EAS')\n"
			r_text += "plot(hardness_a4~n,data=lldat,type='p',col='green',xlab='Step',ylab='Local Hardness')\n"
			r_text += "plot(chg_a4~n,data=lldat,type='p',col='black',xlab='Step',ylab='Partial Charge')\n"
			r_text += "dev.off()\n"
			
			
		r_scr = open("pair_rd.R",'w')
		r_scr.write(r_text)
		r_scr.close()


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

#***********************************************************************    
class primordia_traj:
	def __init__(self,avg,residue_lst):
		self.rds_avg    = None
		self.rds_sd     = None
		self.labels_rd  = "EAS NAS RAS Netphilicity Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic Electrophilic Fukushima Electron_Density Softness_dual MEP"
		self.labels_res = []
		self.res_list   = residue_lst
		self.avg        = avg 
		self.trajs      = glob.glob("*residues.lrd")
		self.rds_values = None
		self.rds        = None
		
		self.traj_len = len(self.trajs)
		self.list_len = len(self.res_list)
		
		self.rds_avg    = np.zeros( (self.list_len,15),dtype=float )
		self.rds_sd     = np.zeros( (self.list_len,15),dtype=float ) 
		self.rds_values = np.zeros( (self.traj_len,self.list_len,15),dtype=float )
		self.rds        = np.zeros( (self.traj_len,300,15),dtype=float )
		
	#===================================================================
	
	def fill_arrays(self):
	
		#---------------------------------------------------------------	
		i = 0
		for fl in self.trajs:
			res_rd = open(fl,'r')
			j = 0
			for line in res_rd:
				line2 = line.split()
				if j > 0:
					self.labels_res.append(line2[0])
					for k in range(14):
						#print(line2[k+1])
						#input()
						self.rds[i][j-1][k] = line2[k+1]					
					#print(self.rds[i][j-1][k])
					#input()				
				j += 1
			i += 1
	

		#---------------------------------------------------------------
		m = 0
		for x in range(self.traj_len):
			for y in range(j):
				if y in self.res_list:
					if m < self.list_len:
						for z in range(15):
							if m < self.list_len:
								print(self.rds[x][y][z])
								self.rds_values[x][m][z] = self.rds[x][y][z]
								print(self.rds_values[x][y][z])
								input()
								
					m+=1
					
			print(self.rds[x][y][z])	
			input()				
					
						
						
		#---------------------------------------------------------------
		
		sumtmp = 0.0		
		for x in range(15):
			for y in range(self.list_len):
				sumtmp = 0.0;
				for z in range(self.traj_len):
					sumtmp += self.rds_values[z][y][x]
					
				self.rds_avg[y][x] = sumtmp/self.traj_len

		devtmp = 0.0
		for x in range(15):
			for y in range(self.list_len):
				devtmp = 0.0
				for z in range(self.traj_len):
					devtmp += (abs(self.rds_values[z][y][x] - self.rds_avg[y][x]))**2
				
				self.rds_sd[y][x] = np.sqrt(devtmp/self.traj_len)
		
	
	#===================================================================
	
	def write_data(self):
		
			
		if not self.avg:
			k = 0
			for res in self.res_list:
				report_fl = open(self.labels_res[int(res)]+"_data_residue",'w')
				report_txt = self.labels_rd
				report_txt += "\n"; 
				for j in range(self.traj_len):
					for i in range(15):
						print(self.rds_values[j][k][i])
						print(j,k,i)
						input()
						report_txt += str(self.rds_values[j][k][i])
						report_txt += " "
					report_txt+="\n"
				report_fl.write(report_txt)
				report_fl.close()						
				k+=1
						

			
		pass
		
	#===================================================================	
	
	def write_R(self):
		pass
		
		
	
	
	
