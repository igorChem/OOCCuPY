#!/usr/bin/env python
# -*- coding: utf-8 -*-
# gamessmodule.py

#=======================================================================
#load modules

import os 
import glob 
import numpy as np
	
#=======================================================================

atomnumber = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':36,'BR':36,'CL':17}

#=======================================================================
# A class to storage xyz molecular information
#=======================================================================

class xyz_parser:
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
				self.AtomLabels.append(line2[0][0])
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
			
		return(text_to_return)

		
	def get_origin(self): 
		
		self.origin[0] = min(self.xCoord)-3
		self.origin[1] = min(self.yCoord)-3
		self.origin[2] = min(self.zCoord)-3
		
		self.xvec[0] = max(self.xCoord) +3
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
		
		print (self.ElNumb)
		return(self.ElNumb)
	
		
class log_parser:
	
	def __init__(self,name):
		self.name = name 
		
		
#=======================================================================
# Function to write gamess input 
#=======================================================================
		
def gms_inp_writer(runtyp = 'Optimize',
				   inp_name ='default' ,
			       charge = 0		 ,
			       xyzfile = 'Normal',
			       multiplicity = 1  ,
			       QMmethod = 'HF'   ,
			       name  = ''        ,
			       semiEmpi = False  ,
			       gbasis = 'N31'    ,
			       ngauss = 3        ,
			       ndfunc = 0        ,
			       npfunc = 0        , 
			       diffuseP = False  ,
			       diffuseS = False  , 
			       maxit = 100       ,
			       mw = 20           ,
			       chelpg = False    ,
			       polar = 'popn31'  ,
			       NatOrb = False    ,
			       opttol = 0.005    ,
			       algopt = "QA"     ,
			       pcm = False      ,
			       elsden = False    ):
	
	
	#------------------------------------------------------------------#
	#======================RUNTYP specifications=======================#
	
	if runtyp == 'hess':
		mv = 200
		
	#------------------------------------------------------------------#
	#==========================system group============================#
		
	system_group = '' 
	system_group += ' $system mwords = {0} $end \n'.format(mw)
	
	#------------------------------------------------------------------#
	#========================control group=============================#
	
	control_group  = ''
	control_group += ' $contrl  runtyp = {0} $end \n'.format(runtyp)
	if QMmethod == 'DFT':
			control_group += ' $contrl  dfttyp = b3lyp $end \n'
	elif QMmethod == 'MP2':
			control_group += ' $contrl mplevel = 2  $end \n'
			
	scf = 'rhf'
	if not multiplicity == 1:
		scf = 'rohf'		
				
	control_group += ' $contrl  icharg = {0} mult = {1}  $end \n'.format(charge,multiplicity)
	control_group += ' $contrl  maxit= {0} $end \n'.format(maxit)
	control_group += ' $contrl scftyp={0} $end \n'.format(scf)
	
	#------------------------------------------------------------------#
	#======================optimize group==============================#
	
	opt_group = ''
	if runtyp == 'Optimize':
		opt_group += ' $statpt opttol = 0.0005 nstep = 150 $end \n'
		opt_group += ' $statpt method = QA $end \n'
	#------------------------------------------------------------------#
	#=========================SCF group================================#
	
	scf_group  = ''
	scf_group += ' $scf dirscf =.true.  $end \n'
	scf_group += ' $scf npreo(1)=0,-1,1,9999 $end \n'
	
	if NatOrb == True:
		scf_group += ' $scf uhfnos=.T. $end \n'

	#------------------------------------------------------------------#
	#=========================Basis group==============================#
	
	if  semiEmpi == False:
		basis_group = ''
		basis_group +=  ' $basis gbasis = {0} ngauss = {1} $end \n'.format(gbasis,ngauss)
		basis_group +=  ' $basis ndfunc = {0} npfunc = {1} $end \n'.format(ndfunc,npfunc)
		if not polar == '': 
			basis_group +=  ' $basis polar = {0} $end \n'.format(polar)
		
		if diffuseP == True:
			basis_group +=  ' $basis diffp =.TRUE.  $end \n'
			basis_group +=  ' $basis diffs =.TRUE. $end \n'
	
	else: 
		basis_group = ''
		basis_group = '$basis gbasis = RM1 $end \n' 

	#------------------------------------------------------------------#
	#===================electrostatic potential group==================#
			
	elpot_group = ''		
	if chelpg == True:
		elpot_group += ' $elpot where=PDC $end \n'
		elpot_group += ' $pdc  ptsel= chelpg $end \n'
		
	#------------------------------------------------------------------#
	#===================electric density group=========================#	
		
	elsden_g = ''	
	cube_grid = ''

	if elsden == True:	
		elsden_g =   ' $eldens  ieden=1 morb=0 where=grid output=punch $end \n'
		cube_grid += ' $grid modgrd=1 $end \n'
		cube_grid += ' $grid size = 0.1 $end \n'
		
		obj1 = xyz_parser(name)
		obj1.parse_xyz()
		obj1.get_origin()		
		a = obj1.elsden_text
		cube_grid += a
		
	#------------------------------------------------------------------#
	#===========================data group=============================#

	data_group =  ' $data \n' 	
	data_group += 'Molecule specification \n'
	data_group += 'c1 \n'
		
	obj = xyz_parser(name)
	obj.parse_xyz()
	obj.get_atomnumber()
	b = obj.write_text()		
	
	data_group += b
	data_group += ' $end'
	
	
	
	#------------------------------------------------------------------#
	solv_group =''
	if pcm == True:
		solv_group += ' $pcm solvnt=h2o $end \n'
		
	input_text = '' 
	input_text = system_group + control_group + opt_group + scf_group + basis_group + elpot_group +solv_group+ elsden_g + cube_grid + data_group
	
	if inp_name == 'default':
		inp_name = name[:-4] + '.inp'
			
	inpfile = open(inp_name,'w')
	inpfile.write(input_text)
	inpfile.close()

#=======================================================================
# Class with methods to automatize the gamess run
#=======================================================================	

class run_class:
	def __init__(self,nprocs = 3):
		self.name = ''
		self.QMmethod = ''
		self.nprocs = nprocs
		
	def optimize(self,xyzname='',QMLevel='low',charge = 0,Alg="QA",Optol=0.0005):
				
		self.name = xyzname
			
		Multiplicity = 1
    
		if not charge%2==0:
			Multiplicity=Multiplicity+1
		
		if QMLevel == 'low': 
			gms_inp_writer(name = xyzname,ngauss=6,algopt=Alg, opttol=Optol) 
		elif QMLevel == 'medium':			
			gms_inp_writer(name = xyzname,gbasis='N311',ngauss=6,algopt=Alg, opttol=Optol) 
		elif QMLevel == 'medDFT':
			gms_inp_writer(name = xyzname,QMmethod='DFT',ndfunc=1,gbasis='N311') 
		elif QMLevel == 'high':
			gms_inp_writer(name = xyzname,QMmethod='DFT',ndfunc=1,npfunc=1,polar='popN311',gbasis='N311') 
		elif QMLevel == 'highPol':
			gms_inp_writer(name = xyzname,QMmethod='DFT',ndfunc=2,npfunc=2,polar='popN311',gbasis='N311') 
		elif QMLevel == 'highDiff':
			gms_inp_writer(name = xyzname,QMmethod='DFT',ndfunc=2,npfunc=2,polar='popN311',gbasis='N311',diffuseP=True,diffuseS=True) 
		
			
	def single_point(self,QMMethod='HF',QMLevel='low',xyzname='',Charge = 0,elecDgrid=False,NOB=False,MW=200):
	
		self.name = xyzname
		
		Multiplicity = 1
	    
		if not Charge%2==0:
			Multiplicity=Multiplicity+1
		
		if QMLevel == 'low':
			gms_inp_writer(runtyp='energy',name = xyzname,QMmethod=QMMethod,charge = Charge,multiplicity=Multiplicity,ngauss=6,gbasis='N31',polar='',mw=MW) 	
		elif QMLevel == 'high':
			gms_inp_writer(runtyp='energy',name = xyzname,QMmethod=QMMethod,ngauss=6,polar='popN311',gbasis='N311',elsden = elecDgrid) 
		elif QMLevel == 'highPol':
			gms_inp_writer(runtyp='energy',name = xyzname,QMmethod=QMMethod,ndfunc=2,npfunc=2,polar='popN311',gbasis='N311',elsden = elecDgrid)
		elif QMLevel == 'highDiff':
			gms_inp_writer(runtyp='energy',name = xyzname,QMmethod=QMMethod,ndfunc=2,npfunc=2,polar='popN311',gbasis='N311',diffuseP=True,diffuseS=True,elsden = elecDgrid,NatOrb=NOB)
	
	def protein_sp(self,QMMethod='HF',Charge=0,xyzname='',n=0,Mult = 1):
	
		self.name = xyzname
		inp_list_file = glob.glob('*.inp')
		inp_Name = xyzname[:-4] +str(n)+'.inp'		
					
		Multiplicity = Mult
	    
			
		gms_inp_writer(runtyp='energy',inp_name=inp_Name,name = xyzname,QMmethod=QMMethod,charge = Charge,multiplicity=Multiplicity,ngauss=3,gbasis='N21',polar='',mw=200,elsden=True,pcm=True) 	
	
	
	def run (self,host='igor'): 
		
		jobname = self.name[:-4]
		text_run = '/home/'+host+'/programs/gamess/rungms ' + jobname + ' 00 ' + str(self.nprocs) + ' > ' + jobname + '.log'
	
		os.system(text_run)	
		obj = xyz_parser(jobname +'.log')
		try:
			obj.parse_log()	
			obj.write_text(filename =jobname +'1.xyz', writefile=True)
		except:
			pass
		
#=======================================================================		
	

'''
try :
	text_to_run = 'sudo cp gamessmodule.py /usr/lib/python3.5' 
	text_to_run = 'sudo cp gamessmodule.py /usr/lib/python2.7' 
	os.system(text_to_run)
except :
	pass		
'''


