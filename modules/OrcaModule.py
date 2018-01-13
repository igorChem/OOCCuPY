#!/usr/bin/env python
# -*- coding: utf-8 -*-
# OrcaModule.py

'''
This python script contains functions and classes to run and parse 
ORCA software (Quantum chemistry calculations program) jobs, and 
calculateglobal and local reaction descriptors as well bases on 
the conceptual density functional theory.  
'''


#=======================================================================
#import needed modules

import os 
import glob 
import numpy as np
	
#=======================================================================

#=======================================================================
# class to retrieve information from orca .out files 
#=======================================================================

class out_file: 
	
	'''
	Class to process and parse output files of orca QM runs 
	'''
		
	def __init__(self,filename):
		
		'''
		Initialization method for the out_file class
		'''
		
		self.name = filename
		self.Energy                    = []
		self.MOoccupiedEnergy          = []
		self.MOunnocupiedEnergy        = []
		self.number_elec               = []
		self.total_enthalpy            = []
		self.total_entropy_correction  = []
		self.free_gibbs_enthalpy       = [] 
		self.LUMO                      = []
		self.HOMO                      = []
		self.nLUMO                     = []
		self.nHOMO		               = []
		self.vib_frequences            = []
		self.GRD                       = [] 
		self.Fukui                     = [] 
		self.LRD                       = []
		self.charges				   = []

	
	def read_out(self):
		
		'''
		Method to extract the information required to compute 
		the Reaction Descriptors from orca output file
		'''
		outputFileName = self.name
		outputFile = open(outputFileName, 'r')
		
		for line in outputFile:
			line = line.split()
			if not line:
				continue
			if len(line) >= 3:	
				if line[0] == 'Total' and line[1] == 'Energy':
					self.Energy.append(round(float(line[3]),4))
				elif line[1] == '0.0000':
					self.MOunnocupiedEnergy.append(float(line[3]))
				elif line[1] == '2.0000':
					self.MOoccupiedEnergy.append(float(line[3]))
				elif line[0] == 'Number' and line[2] == 'Electrons':
					self.number_elec.append(int(float(line[-1]))) 
				elif line[0] == 'Total' and line[1] == 'enthalpy':
					self.total_enthalpy.append(line[3])
				elif line[0] == 'Total' and line[1] == 'entropy':
					self.total_entropy_correction.append(line[4])
				elif line[0] == 'Final' and line[1] == 'Gibbs':				
					self.free_gibbs_enthalpy.append(line[5])

		
		print(self.number_elec[0])			
		numOrbOCC = (self.number_elec[0])/2
			
		MOoccupiedEnergy_fUSE = self.MOoccupiedEnergy[-int(numOrbOCC):]
		
		self.LUMO=self.MOunnocupiedEnergy[0]
		self.HOMO=max(MOoccupiedEnergy_fUSE)					
		
		self.nHOMO = numOrbOCC-1
		self.nLUMO = numOrbOCC 
		
		outputFile.close()
		
	def Global_RD(self): 
		
		''' 
		Method for work over the parsed information of the orca out file
		and calculate the Global Reaction Descriptors. 
		 
		The follpwong commands calculates the reactions descriptors 
		bases on the above parsers
		
		Ionization Potential (IP)
		Electron Affinity (EA) ORbasis
		Electronic Chemical Potential (ECP)
		Chemical Hardness (CH)
		Chemical Softness (CS)
		Muliken Electronegativity (ME)
		Electrophicity Total (EPT)
		Maximum Electrons Reciveble (MER)
		Electrodonating (ED)
		Electroaccepting (EAC)
		Nucleophilicity (N)
		Relative Nucleophilicity (RN)
		HOMO-LUMO energy Difference (gap)
		
		'''
	
		IP   = round(-float(self.HOMO),4)
		EA   = round(-float(self.LUMO),4)
		ECP  = round((-(IP+EA)/2),4)
		CH   = round(((IP-EA)/2),4)
		CS   = round((1/CH),4)
		ME   = -ECP
		EPT  = round(((ECP**2)/(2*CH)),4)
		MER  = round((-ECP/(2*CH)),4)
		ED   = round(((IP**2)/(2*CH)),4)
		EAC  = round(((EA**2)/(2*CH)),4)
		N    = round((10/ED),4)
		gap  = round(float(self.LUMO) - float(self.HOMO),4)
	
		results_text = "name Chemical_Potential Hardness Softness IP EA Electronegativity Electrophilicty MER  Energy_(0 K) gap \n"
		results_text2 = self.name +" "+str(ECP) +" "+ str(CH)+" "+str(CS)+" "+str(IP)+" "+str(EA)+" "+str(ME)+" "+str(EPT)+" "+str(MER)+" "+str(ED)+" "+str(EAC)+" "+str(N)+" "+str(self.Energy[-1])+" "+str(gap) + "\n"
		
		self.GRD = [IP,EA,ECP,CH,CS,ME,EPT,MER,ED,EAC,N,gap]
		
		return(results_text+results_text2,results_text2,IP,EA,ECP,CH,CS,ME,EPT,MER,ED,EAC,N,gap)
		
		
	
	def Freq_Analysis(self): 
		
		'''
		Method that analyze parsed information of the FREQ orca runs for 
		ensure the optmized structure in a local minimum in the potential 
		energy surface
		
		'''
		
		outputFileName = self.name + '.out'
		
		phrase='VIBRATIONAL FREQUENCIES'
		phrase2='NORMAL MODES'
	
		ir_init = []
		ir_fin = []
		
		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if phrase in line:
					ir_init.append(i)
				elif phrase2 in line:
					ir_fin.append(i)
	
		if ir_init[0]>0:			
			with open(outputFileName,'r') as text:
				for (i, line) in enumerate(text):
					if i >= ir_init[0] and i <= ir_fin[0]:
						line2 = line.split()
						if len(line2) == 3:
							self.vib_frequences.append(line2[1])
					else:
						continue
		
	
		imaginary = False				
		print(self.vib_frequences)	
		for i in range(0,len(self.vib_frequences)):
			if float(vib_frequences[i]) < 0:
				print(vib_frequences[i])
				imaginary = True
				textfreq = 'Frequency Analysis \n'
				textfreq += str(self.vib_frequences[i]) + '\n'
				
		if	imaginary == True:			
			FA = open(self.name+'freq_analysis','w')
			FA.write(textfreq)
			FA.close()	
			
		outputFile.close()
	
	def thermochemistry(self):
		
		'''
		Method for calculate statistical thermodynamic properties from 
		FREQ analysis from orca QM run 
		'''	
		
		print(self.total_enthalpy,self.total_entropy_correction,self.free_gibbs_enthalpy)
		return(self.total_enthalpy,self.total_entropy_correction,self.free_gibbs_enthalpy)
	
	def xyz_parse(self,argname):
		
		'''
		Method to create xyz file from orca out file 
		'''
		
		outputFileName = argname + '.out'
		
		xyz01 = []
		xyz02 = []
		xyzCO = []
		
		phrase = 'CARTESIAN COORDINATES (ANGSTROEM)'
		phrase2 = 'CARTESIAN COORDINATES (A.U.)'
		
		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if phrase in line:
					xyz01.append(i)
				elif phrase2 in line:
					xyz02.append(i)
		
		
		with open(outputFileName,'r') as text:
			for (i,line) in enumerate(text):
				if i >= xyz01[0] and i <= xyz02[0]:
					line2 = line.split()
					if len(line2) == 4:
						xyzCO.append(line)
		
		
						
		return(xyzCO)		
		
	def Fukui_Extract (self):
	
	
		'''
		Method to extracts the output information for a Fukui 
		indices	run and calculates the condesed forms of local 
		descriptors
		'''
		
		phrase = 'CHELPG Charges'
		phrase2 = 'Total charge:'
		
		outputFileName = self.name
		
		
		chelp_in=[]
		chelp_out=[]
		
		table_index    = []
		table_atom     = []
		neutro_charge  = []
		cation_charge  = []	
		anion_charge   = []

		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if phrase in line:
					chelp_in.append(i)
				elif phrase2 in line:
					chelp_out.append(i)	
					
		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if i >= chelp_in[0] and i <= chelp_out[0]:
					line2 = line.split()
					if len(line2) == 4:
						table_index.append(line2[0])
						table_atom.append(line2[1])
						neutro_charge.append(round(float(line2[3]),4))
				else:
					continue
					
		neutro_charge = np.array(neutro_charge)
		self.charge = neutro_charge 
			
		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if i >= chelp_in[1] and i <= chelp_out[1]:
					line2 = line.split()
					if len(line2) == 4:
						cation_charge.append(round(float(line2[3]),4))
				else:
					continue
					
		caton_charge = np.array(cation_charge)
		
		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if i >= chelp_in[2] and i <= chelp_out[2]:
					line2 = line.split()
					if len(line2) == 4:
						anion_charge.append(round(float(line2[3]),4))
				else:
					continue
			
		anion_charge = np.array(anion_charge)
		
		#--------------------------------------------------------------#
		
		xyz01 = []
		xyz02 = []
		xyzCO = []
				
		phrase = 'CARTESIAN COORDINATES (ANGSTROEM)'
		phrase2 = 'CARTESIAN COORDINATES (A.U.)'
		
		with open(outputFileName,'r') as text:
			for (i, line) in enumerate(text):
				if phrase in line:
					xyz01.append(i)
				elif phrase2 in line:
					xyz02.append(i)
		
		
		with open(outputFileName,'r') as text:
			for (i,line) in enumerate(text):
				if i >= xyz01[0] and i <= xyz02[0]:
					line2 = line.split()
					if len(line2) == 4:
						xyzCO.append(line)
			
		#---------------------------------------------------------------
			
		E_attack = neutro_charge - cation_charge # Electrophile_attack
		N_attack = anion_charge - neutro_charge # Nucleophile_attack
		R_attack = (cation_charge - anion_charge)/2 # Radical_attack
		
		Fukui_Electrophile = round(max(abs(E_attack)),4) # max site
		Fukui_Nucleophile  = round(max(abs(N_attack)),4) # max site
		Fukui_Radical      = round(max(abs(R_attack)),4) # max site
	
		localElec = round(self.GRD[4]*E_attack,4) # electrophilicity sucetibility
		localNuc  = round(self.GRD[4]*N_attack,4) # nucleophilicity  sucetibility
		
		R_localElec = round(localElec/localNuc,4) # relative electrophilicity sucetibility
		R_localNuc  = round(localNuc/localElec,4)  # relative nucleophilicity  sucetibility
			
		
		localAcidityElec  = round(self.GRD[6]*E_attack,4)  # localAcidityElec
		localAcidityNuc   = round(self.GRD[6]*N_attack,4)  # localAcidityNuc
		ElecDonating      = round(self.GRD[8]*E_attack,4)  # local ElecDonating tendency
		ElecAcceping      = round(self.GRD[9]*E_attack,4)  # local Accepting tendency
		Nucleophilicity   = round(self.GRD[10]*N_attack,4) # local nucleophilicity tendency
		

		
		self.Fukui = [Fukui_Electrophile,Fukui_Nucleophile,Fukui_Radical]
		self.LRD   = [E_attack,N_attack,R_attack,localElec,localAcidityNuc,R_localElec,R_localNuc,ElecDonating,ElecAcceping,Nucleophilicity,xyzCO]
		
		text = self.name + '\n'
		text += 'FukuiE FukuiN FukuiR LocalHardnessE LocalHardnessN \n' 
		text += str(Fukui_Electrophile) +' '+ str(Fukui_Nucleophile) +' '+str(Fukui_Radical) + '\n'
		text += 'E_attack N_attack R_attack localElec localAcidityNuc R_localElec R_localNuc ElecDonating Nucleophilicity xyzCO charge\n'
		for i in range(len(self.LRD[0])):
			text += str(self.LRD[0][i]) +' '+ str(self.LRD[1][i]) +' '+str(self.LRD[2][i]) +' '+ str(self.LRD[3][i]) +' '+str(self.LRD[4][i]) +' '+str(self.LRD[5][i]) +' '+str(self.LRD[6][i]) +' '+str(self.LRD[7][i]) +' '+str(self.LRD[8][i])+' '+ str(self.LRD[9][i])+' '+ str(self.LRD[10][i]) +' '+str(self.charge[i])
		
		resultsfukui = open('resultsfukui','w')
		resultsfukui.write(text)
		resultsfukui.close()
			
		
		return(self.Fukui,self.LRD,text,self.charge)
		

		
#=======================================================================
# Cube Extraction and wrtting methods 
#=======================================================================

class Elec_Cube:
	
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
        
	def Read_Elec_Cube (self,filename):
		
		'''
		Method to parse electronic data cube and fill the atributes
		of Elec_cube class 
		'''
		
		if self.Prog == 'ElecDensity':
			self.name = filename[0:-12]
		elif self.Typ == 'MO':
			self.name = filename[0:-11]			
		
		cubefile = open(filename,'r')				
	
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
		
		num = abs(self.natoms) + 5
		for i in range(5,num):
			self.atoms += cubefile.readline()
			
		cubefile2 = open(filename,'r')		
		str = ' '.join(cubefile2.readlines()[NA:])
		self.scalar3d = np.fromstring(str,sep = ' ')
		self.scalar3d.shape =(self.grid, self.grid, self.grid)


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

		for i in range(self.grid):
			for j in range(self.grid):
				for k in range(self.grid):
					text_to_write += str(self.scalar3d[i][j][k]) + '   '
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
		
		z = Elec_Cube()
		
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



#======================================================================= 
#Descriptors input builder Function
#=======================================================================

def RD_inp	   (name	     	  ,
				QMmethod = "DFT"  ,
				Opt = True	      ,
				nprocs = 1	      ,
				multiplicity = 1  ,
				charge = 0        ,
				basis = "polar"   ,
				Freq = False	  ,
				ThermoC = False   ,
				mode = "Normal"  ):
					
		'''
		Building input function for reaction descriptors calculations 
		in EasyHybrid. The main Goal is to set Orca calculation for 
		isolated gasosous species selected from QM region. When 'fukui'
		is setted to 'True' the local reaction descriptors protocol are
		selected
		'''
		
		#=======================================================#
		# selecting options of calculations input based on      #
		# the function args                                     #
		#=======================================================#
		
		#---------------------------------------------------------------
	
		QMorca=""											 
		if QMmethod=="DFT":									
			QMorca="RKS B3LYP";                              
		elif QMmethod=="MP2":
			QMorca="MP2"
		elif QMmethod =="M062X":
			QMorca="M062X"
			
		#---------------------------------------------------------------
		
		ORbasis=""
		if basis=="apolar":
			ORbasis="6-311G*"
		elif basis=="polar":
			ORbasis="6-311++G**"
		elif basis=="ionic":
			ORbasis="6-311++G(2d,2p)"
		elif basis =="Dunning":
			ORbasis= "cc-pVTZ"
				
		#---------------------------------------------------------------
		
		if nprocs==1:
			FL='#'
		else:
			FL=''
			
		#---------------------------------------------------------------
		
		if Freq == True:
			FR=''; Opt=True
		elif Freq == False:
			FR='#'
			
		#---------------------------------------------------------------
		
				
		if ThermoC == True and mode == "SP":
			TC=''
		elif ThermoC == True and QMmethod == "DFT":
			FR=''
		elif ThermoC == True and not QMmethod == "DFT": 
			TC='';Opt=True
		elif ThermoC == False:			
			TC='#'
		
		#---------------------------------------------------------------
		
		#======================================================#
		# Building the text to write in input files			   #
		#======================================================#
		
		#---------------------------------------------------------------
		# Text building for input to calculate ReationDescriptors,
		# geometr Optimizations and others runs
		#---------------------------------------------------------------
					
		if mode=="Normal":						
			TextbodyOpt  = ''			
			TextbodyOpt += '{0}!PAL{1} \n'.format(FL,nprocs) 
			TextbodyOpt += '!HF Opt 6-31G\n'                            
			TextbodyOpt += '%output print[p_mos] 1 end \n'
			TextbodyOpt += '%base "{0}'.format(name) + '_1opt" \n'                            # first opt run name 
			TextbodyOpt += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name) +'.xyz \n' # first opt run xyz file
			TextbodyOpt += '$new_job \n'
			TextbodyOpt += '{0}!PAL{1} \n'.format(FL,nprocs)
			TextbodyOpt += '!{0} Opt PrintBasis {1} \n'.format(QMorca,ORbasis)
			TextbodyOpt += '{0}!NumFreq TIGHTSCF \n'.format(FR)
			TextbodyOpt += '%output print[p_mos] 1 end \n'
			TextbodyOpt += '%base "{0}_{1}_Opt" \n'.format(name,QMmethod)					   # second opt run name 
			TextbodyOpt += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name) +'_1opt.xyz \n' # second opt run xyz file
			
			if not QMmethod=="DFT":
				TextbodyOpt += '$new_job \n'
		
			Textenergy =''             
			Textenergy += '{0}!PAL{1} \n'.format(FL,nprocs)	
			Textenergy += '!{0} PrintBasis {1}  \n'.format(QMorca,ORbasis)
			Textenergy += '{0}!NumFreq \n'.format(TC)	
			Textenergy += '%output 	print[p_mos] 1 end \n'				
			Textenergy += '%base "{0}"'.format(name) + "\n'"                                # single point run name
			Textenergy += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name)+'.xyz'
			
			if not QMmethod=="DFT":
				final_text = TextbodyOpt + Textenergy
			else:
				final_text = TextbodyOpt
				
		#---------------------------------------------------------------
		# Text building for input to calculate the Condensated Fukui
		# indices
		#---------------------------------------------------------------
			
		elif mode=="Fukui":
			Textenergy =''             
			Textenergy += '{0}!PAL{1} \n'.format(FL,nprocs)
			Textenergy += '!{0} PrintBasis {1} CHELPG \n'.format(QMorca,ORbasis)
			Textenergy += '%output 	print[p_mos] 1 end \n'	
			Textenergy += '%base "{0}'.format(name)+'neutro" \n'			
			Textenergy += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name)+'.xyz \n'
			Textenergy += '$new_job \n'
			
			charge=charge+1
			
			if not charge%2==0:
				multiplicity=multiplicity+1
			             
			Textenergy += '{0}!PAL{1} \n'.format(FL,nprocs)
			Textenergy += '!{0} PrintBasis {1}  CHELPG \n '.format(QMorca,ORbasis)
			Textenergy += '%output 	print[p_mos] 1 end \n'
			Textenergy += '%base "{0}'.format(name)+'cation" \n'
			Textenergy += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name)+'.xyz \n'
			Textenergy += '$new_job \n'
			
			charge=charge-2
			multiplicity=1
			
			
			if not charge%2==0:
				multiplicity=multiplicity+1				
			             
			Textenergy += '{0}!PAL{1} \n'.format(FL,nprocs)
			Textenergy += '!{0} PrintBasis {1} CHELPG \n'.format(QMorca,ORbasis)
			Textenergy += '%base "{0}'.format(name)+'_anion" \n'
			Textenergy += '%output 	print[p_mos] 1 end \n'	
			Textenergy += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name)+'.xyz \n'
			final_text = Textenergy
			
		#---------------------------------------------------------------
		# Text building for input to calculate the Reaction Descriptors
		# by Single point calculations
		#---------------------------------------------------------------
	
		elif mode=="SP":
			Textenergy =''             
			Textenergy += '{0}!PAL{1} \n'.format(FL,nprocs)
			Textenergy += '!{0} PrintBasis {1} \n'.format(QMorca,ORbasis)
			Textenergy += '{0}!Freq \n'.format(TC)	
			Textenergy += '%output print[p_mos] 1   end \n'
			Textenergy += '%base "{0}"'.format(name)+'\n'
			Textenergy += '* xyzfile {0} {1} {2}'.format(charge,multiplicity,name)+'.xyz \n'
			final_text =''
			final_text = Textenergy
			
	    #========================================================#
	    # Naming the files and writting them					 #
	    #========================================================#
		
		if mode=="Fukui":
			Opt=False; inpfile = name + '_fukui.inp'
		else:
			inpfile= name + '.inp'
				
				
		name=open(inpfile,'w')	
		name.write(final_text)
		name.close()
		
	

#======================================================================#
# Function to use orca plot                                            #
#======================================================================#

class OR_Plot:
	
	'''
	Class to use the program orca_plot to generate cube electronic 
	density like from QM orca out file 	
	'''
	
	def __init__(self,out_name):
		
		'''
		initialization method for the OR_plot class 
		'''
		
		self.name = out_name
		self.gbw_name = ""
		self.homo =""
		self.lumo =""
		
		obj = out_file(self.name) # creates an object of out_file class
		obj.read_out() # uses method to parser the nHOMO and nLUMO info
		
		self.homo=obj.nHOMO
		self.lumo=obj.nLUMO
		
		self.gbw_name = self.name[:-4] + ".gbw"
		
		
	def ElDens(self,grid=70):
		
		'''
		Method to make orca_plot generate electronic density cube_file
		'''
		
		inp_Eldens =''	
		inp_Eldens += "orca_plot {0} -i <<EOF\n".format(self.gbw_name)
		inp_Eldens +="1\n" #enter type of plot
		inp_Eldens +="2\n"
		inp_Eldens +="y\n"
		inp_Eldens +="4\n"
		inp_Eldens +=str(grid) + '\n'
		inp_Eldens +="5\n"
		inp_Eldens +="7\n" 
		inp_Eldens +="10\n"
		inp_Eldens +="11\n"
		inp_Eldens += "EOF"
		os.system(inp_Eldens)	
		
	def HOMOd(self,grid=70):
		
		'''
		Method to make orca_plot generate HOMO cube_file
		'''
		
		inp_HOMO  =  "orca_plot {0} -i <<EOF\n".format(self.gbw_name)
		inp_HOMO += "1\n" #enter type of plot
		inp_HOMO += "1\n" #molecular orbitals
		inp_HOMO += "2\n" #enter no of orbital to pot
		inp_HOMO += "{0}\n".format(self.homo)
		inp_HOMO += "4\n" # enter number of grid intervals
		inp_HOMO += str(grid) +"\n" # enter number of grid intervals
		inp_HOMO += "5\n" #select output file formt
		inp_HOMO += "7\n" 
		inp_HOMO += "10\n" #generate the plot
		inp_HOMO += "11\n" #exit the program	
		inp_HOMO += "EOF"
		os.system(inp_HOMO)
	
	def LUMOd(self,grid=70):
		
		'''	
		Method to make orca_plot generate LUMO cube_file
		'''

		inp_LUMO =  "orca_plot {0} -i <<EOF\n".format(self.gbw_name)
		inp_LUMO += "1\n" #enter type of plot
		inp_LUMO += "1\n" #molecular orbitals
		inp_LUMO += "2\n" #enter no of orbital to plot
		inp_LUMO += "{0} \n".format(self.lumo)
		inp_LUMO += "4\n" # enter number of grid intervals
		inp_LUMO += str(grid) +"\n" # enter number of grid intervals
		inp_LUMO += "5\n" #select output file formt
		inp_LUMO += "7\n" 
		inp_LUMO += "10\n" #generate the plot
		inp_LUMO += "11\n" #exit the program
		inp_LUMO += "EOF"	
		os.system(inp_LUMO)


class LRDmaps:
	
	'''
	Class to generate cube files with local reaction descriptor 
	informarion 
	'''
	
	def __init__(self,filename,SPname):
		
		'''
		Initialization method for the LRDmaps class
		'''
		
		self.name   = filename
		self.SPname = SPname
		self.Homo = ""
		self.Lumo = ""
		self.descriptors = []
		
		obj = out_file(self.SPname)
		obj.read_out()
		obj.Global_RD()
		self.descriptors = obj.GRD
		
		self.Homo = obj.nHOMO
		self.Lumo = obj.nLUMO
		self.elec_numb = obj.number_elec
		
	def SUC_Electrophile(self,grid2=70):
		
		'''
		Method to generate cube files with Local reaction descriptors 
		for electrophilic succetibility
		'''
		
		name_cube = self.name
		
		HOMO_cube = Elec_Cube(grid=grid2)	
		HOMO_cube.Read_Elec_Cube(filename=name_cube)
		FF_SUC_Elec = HOMO_cube
		FF_SUC_Elec.write_cubeElec(self.name + "_SucE_elecMap")
		
		electrophilicity  = HOMO_cube
		electrophilicity.scalar3d  = (HOMO_cube.scalar3d**2)*self.descriptors[6]
		electrophilicity.write_cubeElec(self.name + "_ElectrophilicityMapSE")
		
		electrodonating = HOMO_cube
		electrodonating.scalar3d  = (HOMO_cube.scalar3d**2)*self.descriptors[8]
		electrodonating.write_cubeElec(self.name +"_Electrodonating_MapSE")
		
		local_softness_SE = HOMO_cube
		local_softness_SE.scalar3d  = (HOMO_cube.scalar3d**2)*self.descriptors[4]
		local_softness_SE.write_cubeElec(self.name + "_local_softness_SE")
		
	def SUC_Nucleophile(self,grid2=70):
		
		'''
		Method to generate cube files with Local reaction descriptors 
		for nucleophilic succetibility
		'''
		
		name_cube = self.name

		LUMO_cube = Elec_Cube(grid=grid2)
		LUMO_cube.Read_Elec_Cube(filename=name_cube)
		FF_SUC_Nuc = LUMO_cube
		FF_SUC_Nuc.write_cubeElec(self.name + "_SucN_elecMap")
		
		electrophilicity  = LUMO_cube
		electrophilicity.scalar3d  = (LUMO_cube.scalar3d**2)*self.descriptors[6]
		electrophilicity.write_cubeElec(self.name + "_ElectrophilicityMapSN")
		
		Electroaccepting = LUMO_cube
		Electroaccepting.scalar3d = (LUMO_cube.scalar3d**2)*self.descriptors[9]
		Electroaccepting.write_cubeElec(self.name + "_Electroacceptingmap")
				
		local_softSNss_SN = LUMO_cube
		local_softSNss_SN.scalar3d  = (LUMO_cube.scalar3d**2)*self.descriptors[4]
		local_softSNss_SN.write_cubeElec(self.name + "_local_softSNss_SN")
		
		
	def Multiphilic(self,grid2=70,lumoname=""):
		
		'''
		Method to generate multifilic moleculr maps cube file
		'''
		
		name_cube = self.name
		
		HOMO_cube = Elec_Cube(grid=grid2)	
		HOMO_cube.Read_Elec_Cube(filename=name_cube)
		FF_SUC_Elec = HOMO_cube
		FF_SUC_Elec.scalar3d =  HOMO_cube.scalar3d**2
		
		name_cube2 = lumoname

		LUMO_cube = Elec_Cube(grid=grid2)
		LUMO_cube.Read_Elec_Cube(filename=name_cube2)
		FF_SUC_Nuc = LUMO_cube
		FF_SUC_Nuc.scalar3d = LUMO_cube.scalar3d**2
		
		FFdual  = FF_SUC_Nuc - FF_SUC_Elec
		FFdual.write_cubeElec(self.name + "_FFdual")
		
		local_softness_SN = LUMO_cube
		local_softness_SN.scalar3d  = (LUMO_cube.scalar3d**2)*self.descriptors[4]
		
		local_softness_SE = HOMO_cube
		local_softness_SE.scalar3d  = (HOMO_cube.scalar3d**2)*self.descriptors[4]
		
		R_Suc_Nuc = local_softness_SN
		R_Suc_Nuc.scalar3d = local_softness_SN.scalar3d/local_softness_SE.scalar3d
		R_Suc_Nuc.write_cubeElec(self.name + "_Relative_Suc_Nucleophilic")
		
		R_Suc_Elec = local_softness_SE
		R_Suc_Elec.scalar3d = local_softness_SE.scalar3d/local_softness_SN.scalar3d
		R_Suc_Elec.write_cubeElec(self.name + "_Relative_Suc_Elecrophilic")

		
	def Local_Harndess(self,eldens_cube="",grid2=70):
		
		'''
		Methodfor generate cube files for mapping local Hardness
		'''
		
		name_cube = self.name
		name_cube2 = eldens_cube
		
		HOMO_cube = Elec_Cube(grid=grid2)	
		HOMO_cube.Read_Elec_Cube(filename=name_cube)
		HOMO_cube.scalar3d = HOMO_cube.scalar3d**2
		
		Eldens_cube = Elec_Cube(grid=grid2,Typ="ELecDensity")
		Eldens_cube.Read_Elec_Cube(filename=name_cube2)
		
		mu = self.descriptors[2]
		eta = self.descriptors[3]
		N_Elec = float(self.elec_numb[0])
		
		local_hardness = HOMO_cube
		local_hardness.scalar3d = (HOMO_cube.scalar3d - (Eldens_cube.scalar3d/N_Elec))*(mu/(2*N_Elec)) + (Eldens_cube.scalar3d/N_Elec)*eta
		local_hardness.write_cubeElec(self.name + "_localhardness")
		
		
		
#======================================================================#
#update in the python3.* libraries                                     #
#======================================================================#

'''
os.system('sudo cp ~/Dropbox/RDscripts/orca module/OrcaModule.py /usr/lib/python3.5')
'''
	

	
	
