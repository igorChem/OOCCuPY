#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pdb_Reader.py

import os,glob

def primordia_inp(option=3,program="mopac",lh="potential_fukui",gridn=0,eband=5,bandm="BD",norb=100):

	f = open("input_pri",'w')
	text ="eband {0} pymols\n".format(eband)


	if 	 option == '1':
		if program   == "mopac":
			aux = glob.glob("*aux")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)
		elif program == "gamess":
			aux = glob.glob("*log")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)
		elif program == "orca":
			aux = glob.glob("*out")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)
		elif program == "gaussian":
			aux = glob.glob("*fchk")
			for ax in aux:
				text+="1 {0} {1} {2} {3} \n".format(ax,lh,gridn,program)		
	elif option == '2':
		if program   == "mopac":
			aux = glob.glob("*aux")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		elif program == "gamess":
			aux = glob.glob("*log")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		elif program == "orca":
			aux = glob.glob("*out")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		elif program == "gaussian":
			aux = glob.glob("*fchk")
			for i in sorted(aux,reverse=True):
				if aux[i][-6:-4] == "cat":
					del aux[i]
				elif aux[i][-7:-4] == "an":
					del aux[i]			
			for i in range(len(aux)):
				text+="2 {0} {1} {2} {3} {4} {5} \n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
	elif option == '3':
		aux = glob.glob("*aux")		
		for i in range(len(aux)):
			text+="3 {0} {1} {2} {3} {4} {5} 0 0 0 0 {6}\n".format(aux[i],lh,gridn,norb,aux[i][:-4]+".pdb",program,bandm)
	
	elif option == '4':
		aux = glob.glob("*aux")
		alist =[]

		for i in range(len(aux)):
			if aux[i][-7:-4] == "cat":
				alist.append(i)
			elif aux[i][-6:-4] == "an":
				alist.append(i)

		for i in sorted(alist,reverse=True):
				del aux[i]
	
		log = glob.glob("*log")
		blist =[]
		for i in range(len(log)):
			if log[i][-6:-4] == "an":
				blist.append(i)
			elif log[i][-7:-4] == "cat":
				blist.append(i)
				
		for i in sorted(blist,reverse=True):			
			del log[i]

		for i in range(len(aux)):
			text+="1 {0} {1} {2} {3} \n".format(aux[i],lh,gridn,"mopac")
		for i in range(len(log)):
			text+="1 {0} {1} {2} {3} \n".format(log[i],lh,gridn,"gamess")
		clist = []
		for i in range(len(aux)):
			if aux[i][-6:-4] == "zy":
				clist.append(i)
		for i in sorted(clist,reverse=True):
				del aux[i]
		for i in range(len(aux)):
			text+="2 {0} {1} {2} {3} {4} {5} {6}\n".format(aux[i][:-4]+".mgf",aux[i][:-4]+"_cat.mgf",aux[i][:-4]+"_an.mgf",lh,gridn,2,"mopac")
		for i in range(len(log)):
			text+="2 {0} {1} {2} {3} {4} {5} {6}\n".format(log[i],log[i][:-4]+"_cat.log",log[i][:-4]+"_an.log",lh,gridn,2,"gamess")


	f.write(text)
	f.close()



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
		self.gstep = []
		self.lstep = []
		self.CT_p1_p2 = []
		self.SPI_p1_p2 = []
		self.HPI_p1_p2 = []
		self.mode = mode 
		self.dist1  = []
		self.dist2  = []
		

		globaln = glob.glob("*global")
		print(globaln[0])
		fgl = open(globaln[0],'r')
	
		#read global descriptors 
		#------------------------------------------------------------------#
		i=0
		for line in fgl:
			if i > 1:
				line2 = line.split()
				self.gstep.append(int(line2[0][6:-7]))
				self.HOF.append(float(line2[2]))
				self.Elec_en.append(float(line2[1]))
				self.hardness.append(float(line2[6]))
				self.ECP.append(float(line2[5]))
				self.electrophilicity.append(float(line2[7]))
			i +=1
		fgl.close()
		
		for j in range(len(self.prds)):
			prd = open(self.prds[j],'r')
			self.lstep.append(int(self.prds[j][6:-15]))
			i = 0
			j = 0
			for line in prd:
				line2 = line.split()
				if i == 0:
					if len(line2) == 7:
						self.pairs = 2
				if i > 0:
					j +=0
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
		
		shof = min(self.HOF)
		selec = min(self.Elec_en)
		for i in range(len(self.gstep)):
			self.HOF[i] = self.HOF[i] - shof
			self.Elec_en[i] = self.Elec_en[i] - selec
		
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
					
		
		fgr_text = "n HOF Energy Hardness ECP Electrophilicity \n"
		fgr = open("global_resume_data",'w')
		for i in range(len(self.gstep)):
			fgr_text += "{} {} {} {} {} {}\n".format(self.gstep[i],self.HOF[i],self.Elec_en[i],self.hardness[i],self.ECP[i],self.electrophilicity[i])
		fgr.write(fgr_text)
		fgr.close()
		
		flr_text = ""
		flr = open("local_resume_data",'w')
		if self.pairs == 1:
			flr_text = "n eas_a1 nas_a1 hardness_a1 chg_a1 chg_a2 eas_a2 nas_a2 hardness_a2 CT SPI HPI\n"
			for i in range(len(self.lstep)):
				flr_text += "{} {} {} {} ".format(self.lstep[i],self.eas_a1[i],self.nas_a1[i],self.hardness_a1[i])
				flr_text += "{} {} ".format(self.chg_a1[i],self.chg_a2[i])
				flr_text += "{} {} {} {} {} {}\n".format(self.eas_a2[i],self.nas_a2[i],self.hardness_a2[i],self.CT_p1[i],self.SPI_p1[i],self.HPI_p1[i])
		elif self.pairs == 2:
			flr_text = "n eas_a1 nas_a1 hardness_a1 eas_a2 nas_a2 hardness_a2 eas_a3 nas_a3 hardness_a3 eas_a4 nas_a4 hardness_a4 chg_a1 chg_a2 chg_a3 chg_a4 CT1 SPI1 HPI1 CT2 SPI2 HPI2\n"
			for i in range(len(self.lstep)):	
				flr_text += "{} {} {} {}".format(self.lstep[i],self.eas_a1[i],self.nas_a1[i],self.hardness_a1[i])
				flr_text += "{} {} {} ".format(self.eas_a2[i],self.nas_a2[i],self.hardness_a2[i])
				flr_text += "{} {} {} ".format(self.eas_a3[i],self.nas_a3[i],self.hardness_a3[i])
				flr_text += "{} {} {} ".format(self.eas_a4[i],self.nas_a4[i],self.hardness_a4[i])
				flr_text += "{} {} {} {} ".format(self.chg_a1[i],self.chg_a2[i],self.chg_a3[i],self.chg_a4[i])
				flr_text += "{} {} {} {} {} {} \n".format(self.CT_p1[i],self.SPI_p1[i],self.HPI_p1[i],self.CT_p2[i],self.SPI_p2[i],self.HPI_p2[i])
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

    
     
