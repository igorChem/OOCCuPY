#!/usr/bin/env python
# -*- coding: utf-8 -*-
# OrcaScript.py


#=======================================================================
#load modules

import os 
import glob 
import numpy as np

from OrcaModule import *

#=======================================================================

#=======================================================================
# Function to gererate cube files for the GBW in the dir
#=======================================================================

def Or_plot_run(Grid=40):
	
	'''
	'''
	
	listgbw = glob.glob('*.gbw')
	
	for i in range(len(listgbw)):		
		name_out  = listgbw[i]		
		OR_plot(name_out,grid=Grid,GBW=True)
		
			
	listcube = glob.glob('*eldens.cube')
	
	for i in range(len(listcube)):
		name_cube = listcube[i]
		if 'neutro' in name_cube:            
			name_cube_neutro = Elec_Cube(Typ='ElecDensity',grid = Grid)
			name_cube_neutro.Read_Elec_Cube(name_cube)
			name_cube_GRD = out_file(name_cube[:-18]+'.out')
			name_cube_GRD.read_out()
			name_cube_GRD.Global_RD()
			softness = name_cube_GRD.GRD[4]
			Electrophilicity = name_cube_GRD.GRD[8]
			Nucleofilicity = name_cube_GRD.GRD[9]
			Hardness = name_cube_GRD.GRD[3]
		elif 'cation' in name_cube:			
			name_cube_cation = Elec_Cube(Typ='ElecDensity',grid = Grid)
			name_cube_cation.Read_Elec_Cube(name_cube)
		elif 'anion' in name_cube:
			name_cube_anion = Elec_Cube(Typ='ElecDensity',grid = Grid)
			name_cube_anion.Read_Elec_Cube(name_cube)
			

	fukuiE = name_cube_neutro - name_cube_cation
	fukuiE.write_cubeElec(fukuiE.name + 'fukuiE')
	fukuiN = name_cube_anion - name_cube_neutro
	fukuiN.write_cubeElec(fukuiN.name + 'fukuiN')

	LSE = fukuiE
	LSE.scalar3d = softness*fukuiE.scalar3d
	LSE.write_cubeElec(LSE.name + 'LSE')
	
	LSN = fukuiN
	LSN.scalar3d = softness*fukuiN.scalar3d
	LSN.write_cubeElec(LSN.name + 'LSN')
	
	Elec = fukuiE
	Elec.scalar3d = Electrophilicity*fukuiE.scalar3d
	Elec.write_cubeElec(Elec.name + 'Elec')
	
	Nuc = fukuiN
	Nuc.scalar3d = Electrophilicity*fukuiN.scalar3d
	Nuc.write_cubeElec(Nuc.name + 'Nuc')
	
	LHelec = fukuiE
	LHelec.scalar3d = Hardness*fukuiE.scalar3d
	LHelec.write_cubeElec(LHelec.name + 'LHelec')
	
	LHnuc  = fukuiN
	LHnuc.scalar3d = Hardness*fukuiN.scalar3d
	LHnuc.write_cubeElec(LHnuc.name + 'LHnuc')
		
	return (fukuiE,fukuiN,LSE,LSN,Elec,Nuc)



#======================================================================#
# Function to extract RD results and write these files to a single file#
#======================================================================#

def RD_Extract_run ():
	
	'''
	Extract RD results from all .out files in a directory and writes in 
	a single file
	
	'''
	listout = glob.glob('*.out')

	text = 'RD results \n'
	text += "Name Chemical_Potential Hardness Softness Ionization_Potential Electron_Affinity Electronegativity Electrophilicty Max_electron_acceptable  elecDonating Nucleophilicity Energy_(0 K) gap \n"
	
	for i in range(0,len(listout)):
		obj = out_file(listout[i])
		print(obj.name)
		obj.read_out()
		a = obj.Global_RD()
		text += a[1]
	
	results = open('RDcompiled','w')
	results.write(text)	
	results.close()
			
		
#======================================================================#
# Function to extract fukui indices for all files in a directory       #
#======================================================================#

def Fukui_Extract_run(name):
	
	'''
	Function to use the Fukui_Extract function to all .out files 
	in the directory
	'''
	
	# condensed fukui table 
	
	obj=out_file(name)
	obj.read_out()
	obj.Global_RD()
	obj.Fukui_Extract()
			    
			    	
#======================================================================#
# Function to generate multiple input and run in terminal              #
#======================================================================#
	
def Orca_Run (MODE="Normal"  ,
			  BASIS="polar"  ,
			  QMMETHOD="DFT" , 
			  freq=False     ,
			  TC=False	     ,
			  Charge = 0    ):
				  
	'''
	Function to run all the functions to calculate the reaction 
	descriptors to the molecules in a directory 
	'''

	Multiplicity = 1
    
	if not Charge%2==0:
		Multiplicity=Multiplicity+1

	xyztrue = glob.glob('*xyz')
	
	for i in range(0,len(xyztrue)):
		infilename = xyztrue[i]
		infilename = infilename[0:-4]
		name = infilename		
		RD_inp(name			       ,
		mode  =  MODE              ,
		basis = BASIS              ,
		QMmethod = QMMETHOD        ,
		charge = Charge            ,
		multiplicity = Multiplicity,
		Freq = freq                ,
		ThermoC = TC               )		
		command = ''
		if MODE == "Fukui":
			command = 'orca ' +name+ '_fukui.inp > ' +name+ '_fukui.out'
		else: 
			command = 'orca ' +name+ '.inp > ' +name+ '.out'
		print(command)
		os.system(command)
	
	for file in glob.glob('*prop') + glob.glob('*tmp') + glob.glob('*trj') +glob.glob('*engrad'):
		os.remove(file)	
		
	if MODE == "Fukui":
		Fukui_Extract_run()	
		listfukui = glob.glob('*fukui.out')
		for i in range(len(listfukui)):
			name_fukuiout = listfukui[i]
			newdir = 'mkdir ' + name_fukuiout[0:-4] +'dir'
			os.system(newdir)
			mynewdir = name_fukuiout[0:-4] + 'dir'
			name_molecule = name_fukuiout[0:-10]
			os.system('mv ' + name_molecule +'* '+ mynewdir)
	elif MODE == "Normal":
		RD_Extract_run()
		
#=======================================================================


def Orca_Run (MODE="SP"  ,
			  BASIS="Dunning"  ,
			  QMMETHOD="M062X" , 
			  freq=False     ,
			  TC=False	     ,
			  Charge = 0    ):
				  
	'''
	Function to run all the functions to calculate the reaction 
	descriptors to the molecules in a directory 
	'''

	Multiplicity = 1
    
	if not Charge%2==0:
		Multiplicity=Multiplicity+1

	xyztrue = glob.glob('*xyz')
	
	
	text_to_write = "name Chemical_Potential Hardness Softness IP EA Electronegativity Electrophilicty MER  Energy_(0 K) gap \n"

	
	for i in range(0,len(xyztrue)):
		infilename = xyztrue[i]
		infilename = infilename[0:-4]
		name = infilename		
		RD_inp(name			       ,
		mode  =  MODE              ,
		basis = BASIS              ,
		QMmethod = QMMETHOD        ,
		charge = Charge            ,
		multiplicity = Multiplicity,
		nprocs = 8				   ,
		Freq = freq                ,
		ThermoC = TC               )		
		command = '/home/barden/programs/orca/orca ' + name +'.inp > '+ name + '.out'
		os.system(command)
		a = out_file(name + '.out')
		a.read_out()
		b = a.Global_RD()
		text_to_write += b[1]
		
	d = open("results",'w')
	d.write(text_to_write)
		
	for file in glob.glob('*prop') + glob.glob('*tmp') + glob.glob('*trj') +glob.glob('*engrad'):
		os.remove(file)	
		
		
		
		
		
try :
	text_to_run = 'sudo cp OrcaScripts.py /usr/lib/python3.5' 
	os.system(text_to_run)
except :
	pass		

		
try :
	text_to_run = 'sudo cp OrcaScripts.py /usr/lib/python3.4'
	os.system(text_to_run)
except :
	pass

try :
	text_to_run = 'sudo cp OrcaScripts.py /usr/lib/python2.7'
	os.system(text_to_run)
except :
	pass

