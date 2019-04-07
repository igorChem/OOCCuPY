from Bio.PDB import *
import math, glob
import os, sys, commands, subprocess
import argparse, textwrap
from tqdm import tqdm
import numpy as np
import pandas as pd

def ConvertingXTCFileToPDBFile():
	""" Function doc """

	if args.trajin[-4:] == ".pdb":
    		
		pass 
	
	elif args.trajin[-4:] == ".xtc":

		os.system("gmx trjconv -f "+ args.trajin +" -s "+ args.tprfile + " -pbc nojump -dt "+ args.time +" -o "+ args.trajin[:-4] +".pdb")

		args.trajin = args.trajin[:-4] +".pdb"

	else:
    	
		print """
		
		##########################--Thank you for choose the Gambiarra Script Package--############################
		#                                                                                                         #
		#                      Please check your input file. It should be a PDB or XTC file                       #
		#                                                                                                         # 
		###########################################################################################################
		"""
		sys.exit()

def TrajectorySplitter():
	""" Function doc """									
	
	# Creating a temporary folder
	if not os.path.isdir("temp"):
    	 os.mkdir("temp")

	if os.path.isfile("temp/complex_1.pdb") == True:
		 pass 
	
	else:

		# Get the number of MODELS in the Multi-PDB file
		status, models = commands.getstatusoutput("grep -c MODEL " + args.trajin)
		modelNum = 1
		new_text = ""

		# Setting the progress bar
		pbar = tqdm(total = (int(models) - 1), desc="Extracting PDB files")
		for line in open(args.trajin, "r"):
			line = line.strip()
			if line == "ENDMDL":
				output = open("temp/complex_" + str(modelNum) + ".pdb", "w")
				output.write(new_text.strip())
				output.close()
				modelNum += 1
				new_text = ""

				# Updating progress bar
				pbar.update(1)

			else:
				new_text += line + "\n"
		pbar.close()

def GetAngleBiopython():
	""" Function doc """
	
	dataAll = [] 
	dataArea = []
	dataAngle = []
	dataDihedral = []
	dataRMSD = []
	dataRG = []
	areaLabels = 'Frame\tDist_AB\tDist_AC\tDist_BC\tArea\n'
	angleLabels = 'Frame\tAngle\n'
	dihedralLabels = 'Frame\tDihedral\n'
	rmsLabels = 'Frame\tRMSD\n'
	rgLabels = 'Frame\tRG\n'
	allLabels = 'Frame\tDist_AB\tDist_AC\tDist_BC\tArea\tAngle\tDihedral\tRMSD\tRG\n'
	
	# Calculate the differntial geometry for all residues (arc length, curvature, writhing, phi and psi)
	if args.xgeo == "xgeo":
			
		xgeoCommand = "xgeo --dir="+os.getcwd()+"/temp"
		subprocess.call(xgeoCommand.split())
		files = glob.glob("./temp/*.tbl")
		table = pd.read_table(files[0], delim_whitespace=True)
		for i in range(1,len(files)):
			bufferT = pd.read_table(files[i], delim_whitespace=True)
			table = table.append(bufferT, ignore_index=True)
		table.to_csv(os.getcwd()+"/"+args.xgeoOut+".tbl", sep="\t")

	else:
    	
		parser = PDBParser(PERMISSIVE=1)

		# Setting refence structure to RMSD calculation
		refStruct = parser.get_structure("temp/complex_1.pdb",  os.getcwd()+"/temp/"+"complex_1.pdb")
		refStruct = refStruct[0]
		refModelChain = refStruct[args.chain]
		refCaAtoms = []
		for refRes in refModelChain:
			resnum = refRes.get_id()[1]
			refCaAtoms.append(refRes['CA'])

		# Setting the progress bar2
		status, models = commands.getstatusoutput("grep -c MODEL " + args.trajin)
		pbar2 = tqdm(total = int(models) - 1, desc="Running "+args.option+" calculation")

		for frame in range(1 ,int(models)):
			
			model = parser.get_structure("temp/complex_" + str(frame) + ".pdb",  os.getcwd()+"/temp/"+"complex_" + str(frame) + ".pdb")

			# Updating progress bar2
			pbar2.update(1)

			# Calculate pincer angle
			if args.option == "angle":
				
				angleRes = args.angle.split(',')
				
				atom1 = model[0][args.chain][int(angleRes[0])]['CA']
				atom2 = model[0][args.chain][int(angleRes[1])]['CA']
				atom3 = model[0][args.chain][int(angleRes[2])]['CA']
			
				vector1 = atom1.get_vector()
				vector2 = atom2.get_vector()
				vector3 = atom3.get_vector()
				
				angle = calc_angle(vector1, vector2, vector3)
				angleLabels += "%i\t%.2f\n" %(frame,math.degrees(angle))
			
			# Calculate dihedral angle	
			elif args.option == "dihedral":
				
				dihedralRes = args.dihedral.split(',')
				
				atom1 = model[0][args.chain][int(dihedralRes[0])]['CA']
				atom2 = model[0][args.chain][int(dihedralRes[1])]['CA']
				atom3 = model[0][args.chain][int(dihedralRes[2])]['CA']
				atom4 = model[0][args.chain][int(dihedralRes[3])]['CA']
				
				vector1 = atom1.get_vector()
				vector2 = atom2.get_vector()
				vector3 = atom3.get_vector()
				vector4 = atom4.get_vector()
				
				dihedral = calc_dihedral(vector1, vector2, vector3, vector4)
				dihedralLabels += "%i\t%.2f\n" %(frame,math.degrees(dihedral))
			
			# Calculate the triangle area according to heron equation	
			elif args.option == "triangle_area":
				
				areaRes = args.area.split(',')
				
				atom1 = model[0][args.chain][int(areaRes[0])]['CA']
				atom2 = model[0][args.chain][int(areaRes[1])]['CA']
				atom3 = model[0][args.chain][int(areaRes[2])]['CA']
				
				coordAB = model[0][args.chain][int(areaRes[0])]['CA'].coord - model[0][args.chain][int(areaRes[1])]['CA'].coord
				coordAC = model[0][args.chain][int(areaRes[0])]['CA'].coord - model[0][args.chain][int(areaRes[2])]['CA'].coord
				coordBC = model[0][args.chain][int(areaRes[1])]['CA'].coord - model[0][args.chain][int(areaRes[2])]['CA'].coord
				
				distAB = np.sqrt(np.sum(coordAB * coordAB))
				distAC = np.sqrt(np.sum(coordAC * coordAC))
				distBC = np.sqrt(np.sum(coordBC * coordBC))
				
				sPerimeter = (distAB + distAC + distBC) / 2
				
				area = (sPerimeter*(sPerimeter-distAB)*(sPerimeter-distAC)*(sPerimeter-distBC)) ** 0.5
				areaLabels += "%i\t%.2f\t%.2f\t%.2f\t%.2f\n" %(frame, distAB, distAC, distBC, area)
	
			# Calculate RMSD			
			elif args.option == "rmsd":
    				
				modelStruct = model[0]
				ModelChain = modelStruct[args.chain]
				modelCaAtoms = []
				for modelRes in ModelChain:
					resnum = modelRes.get_id()[1]
					modelCaAtoms.append(modelRes['CA'])
				sup = Superimposer()
				sup.set_atoms(refCaAtoms, modelCaAtoms)
				rmsd = sup.rms
				rmsLabels += "%i\t%.2f\n" %(frame,rmsd)
				
			# Calculate Radius of Gyration (Rg)
			elif args.option == "rg":
    				
				m = model[0]
				
				allCA = []
				for chain in m:
					for residue in chain:
						for atom in residue:
							if atom.get_id() == 'CA':
								caCoord = atom.get_coord()
								allCA.append(caCoord)

				CoM = sum(allCA)/(len(allCA))

				allCA2 = []
				for chain in m:
					for residue in chain:
						for atom in residue:
							if atom.get_id() == 'CA':
								caCoord2 = atom.get_coord() - CoM
								quadCoord = (caCoord2)**2 * 12.04
								allCA2.append(quadCoord)

				rg = math.sqrt((np.sum(allCA2))/((len(allCA2))*12.04))
				rgLabels += "%i\t%.2f\n" %(frame,rg)

			# Perform all analysis
			elif args.option == "all":
				
				# Triangle residues
				areaRes = args.area.split(',')
				
				atom1 = model[0][args.chain][int(areaRes[0])]['CA']
				atom2 = model[0][args.chain][int(areaRes[1])]['CA']
				atom3 = model[0][args.chain][int(areaRes[2])]['CA']
				
				coordAB = model[0][args.chain][int(areaRes[0])]['CA'].coord - model[0][args.chain][int(areaRes[1])]['CA'].coord
				coordAC = model[0][args.chain][int(areaRes[0])]['CA'].coord - model[0][args.chain][int(areaRes[2])]['CA'].coord
				coordBC = model[0][args.chain][int(areaRes[1])]['CA'].coord - model[0][args.chain][int(areaRes[2])]['CA'].coord
				
				distAB = np.sqrt(np.sum(coordAB * coordAB))
				distAC = np.sqrt(np.sum(coordAC * coordAC))
				distBC = np.sqrt(np.sum(coordBC * coordBC))
				
				sPerimeter = (distAB + distAC + distBC) / 2
				
				area = (sPerimeter*(sPerimeter-distAB)*(sPerimeter-distAC)*(sPerimeter-distBC)) ** 0.5
				
				# Dihedral residues
				dihedralRes = args.dihedral.split(',')
				
				atom1 = model[0][args.chain][int(dihedralRes[0])]['CA']
				atom2 = model[0][args.chain][int(dihedralRes[1])]['CA']
				atom3 = model[0][args.chain][int(dihedralRes[2])]['CA']
				atom4 = model[0][args.chain][int(dihedralRes[3])]['CA']
				
				vector1 = atom1.get_vector()
				vector2 = atom2.get_vector()
				vector3 = atom3.get_vector()
				vector4 = atom4.get_vector()
				
				dihedral = calc_dihedral(vector1, vector2, vector3, vector4)

				# Angle residues
				angleRes = args.angle.split(',')
				
				atom1 = model[0][args.chain][int(angleRes[0])]['CA']
				atom2 = model[0][args.chain][int(angleRes[1])]['CA']
				atom3 = model[0][args.chain][int(angleRes[2])]['CA']
			
				vector1 = atom1.get_vector()
				vector2 = atom2.get_vector()
				vector3 = atom3.get_vector()
				
				angle = calc_angle(vector1, vector2, vector3)

				# Calculate RMSD
				modelStruct = model[0]
				ModelChain = modelStruct[args.chain]
				modelCaAtoms = []
				for modelRes in ModelChain:
					resnum = modelRes.get_id()[1]
					modelCaAtoms.append(modelRes['CA'])
				
				sup = Superimposer()
				sup.set_atoms(refCaAtoms, modelCaAtoms)
				rmsd = sup.rms

				# Calculate RG
				m = model[0]

				allCA = []
				for chain in m.get_list():
					for residue in chain.get_list():
						ca = residue["CA"]
						caCoord = ca.get_coord()
						allCA.append(caCoord)

				CoM = sum(allCA)/(len(allCA))

				allCA2 = []
				for chain in m.get_list():
					for residue in chain.get_list():
						ca = residue["CA"]
						caCoord2 = ca.get_coord() - CoM
						quadCoord = (caCoord2)**2 * 12.04
						allCA2.append(quadCoord)

				rg = math.sqrt((np.sum(allCA2))/((len(allCA2))*12.04))
				
				allLabels += "%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" %(frame, distAB, distAC, distBC, area, math.degrees(angle), math.degrees(dihedral), rmsd, rg)

		pbar2.close()

	# Save Triangle area file
	if args.option == "triangle_area":
		
		dataArea.append(areaLabels)
		Triangleoutput = open(args.output+".tbl", "w")
		Triangleoutput.writelines(dataArea)	
		Triangleoutput.close()	

	# Save Angle file
	if args.option == "angle":
		
		dataAngle.append(angleLabels)
		Angleoutput = open(args.output+".tbl", "w")
		Angleoutput.writelines(dataAngle)	
		Angleoutput.close()
		
	# Save Dihedral file
	if args.option == "dihedral":
		
		dataDihedral.append(dihedralLabels)
		Dihedraloutput = open(args.output+".tbl", "w")
		Dihedraloutput.writelines(dataDihedral)	
		Dihedraloutput.close()

	# Save RMSD file
	if args.option == "rmsd":
		
		dataRMSD.append(rmsLabels)
		RMSDoutput = open(args.output+".tbl", "w")
		RMSDoutput.writelines(dataRMSD)	
		RMSDoutput.close()
	
	# Save RG file
	if args.option == "rg":
    		
		dataRG.append(rgLabels)
		RGoutput = open(args.output+".tbl", "w")
		RGoutput.writelines(dataRG)	
		RGoutput.close()

	# Save Triangle area, Angle and Dihedral files
	if args.option == "all":
		
		dataAll.append(allLabels)
		alloutput = open(args.output+".tbl", "w")
		alloutput.writelines(dataAll)	
		alloutput.close()

	# Delete temporary folder
	#os.system("rm -rf temp/")
		
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='''
	
	\t\t\t\t############################################
	\t\t\t\t#                                          #
	\t\t\t\t#         Gambiarra Package Presents       #      
	\t\t\t\t#                                          #
	\t\t\t\t############################################\n\n
	The "Geometric Analysis" script that was developed to carry out geometric analysis on protein structures.
	This script support as input MultiPDB and XTC files, and the avaliable analysis are: 
	
	1 - Pincer angle.
	2 - Dihedral angle.
	3 - Triangle area.
	4 - RMSD
	5 - RG
	6 - Differential geometry - using the program XGEO - XGeo 0.1 [64 bits] - Copyright (c) 2015-2017 Rinaldo Wander Montalvao.
	''', formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('-m', '--TrajFile\t\t\t', action='store', required=False, dest='trajin', help="The trajectory must be a Multi-PDB or XTC file.\n")
	parser.add_argument('-s', '--TPRFile\t\t\t', action='store', required=False, dest='tprfile', help="The '.tpr' file from from GROMACS package.\n")
	parser.add_argument('-dt', '--WritePDB\t\t\t', action='store', required=False, dest='time', help="Frequency to write the MultiPDB file.\n")
	parser.add_argument('-t', '--analysisType\t\t', action='store', required=False, dest='option', type = str, help='''Specify the type of your analysis: "angle", "dihedral", "triangle_area", "rmsd", "rg" or "all". 
	If you select "all" the flags -rA, -rD and -tA must be in your command line. 
	RMSD use the first frame as reference structure.\n''')
	parser.add_argument('-xg', '--differentialType\t\t', action='store', required=False, dest='xgeo', type = str, help='''Differential geometry analysis using the program XGEO. 
	XGeo 0.1 [64 bits] - Copyright (c) 2015-2017 Rinaldo Wander Montalvao. 
	Just specify "xgeo" as argument.\n''')
	parser.add_argument('-c', '--proteinChain\t\t', action='store', dest='chain', default=' ', help="The chain of the PDB file.\n")
	parser.add_argument('-o', '--OutputName\t\t', action='store', dest='output', help="The output name without extension.\n")
	parser.add_argument('-ox', '--xgeoOutput\t\t', action='store', dest='xgeoOut', default='xgeo_analysis',help="The output name without extension, when you are using differential geometry analysis.\n")
	parser.add_argument('-rA', '--angleResidues\t\t', action='store', dest='angle', 
	                    help=textwrap.dedent('''\
	                    
Select the residues that will be used to measure the pincer angle. 
    Make sure to separate the residues by comma (ex. 30,45,60)

                  (1)      (3)
                    \      /
                     \    /
                      \  /
                       \/ 
                      (2)
	
	'''))
	parser.add_argument('-rD', '--dihedralResidues\t', action='store', dest='dihedral', 
	                    help=textwrap.dedent('''\

Select the residues that will be used to measure the dihedral angle. 
    Make sure to separate the residues by comma (ex. 30,45,60,70)

                  (1)         (4)
                   |           |
                   |           |
                   |           |
                   |           |
                  (2)---------(3)
	
	'''))
	parser.add_argument('-tA', '--areaResidues\t\t', action='store', dest='area',
	                    help=textwrap.dedent('''\
	                    
Select the residues that will be used to measure the triangle area. 
    Make sure to separate the residues by comma (ex. 30,45,60)
	                    
                       AC
                 (1)--------(3)
                    \      /
                AB   \    / BC
                      \  /
                       \/
                      (2)	
	
	'''))
	args = parser.parse_args()
	ConvertingXTCFileToPDBFile()
	TrajectorySplitter()
	GetAngleBiopython()
	
