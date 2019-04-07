from Bio.PDB import *
import os, sys, commands, subprocess
from tqdm import tqdm
import argparse, textwrap
import math
from array import *
import numpy as np

def ConvertingXTCFileToPDBFile():
    	""" Function doc """

	if args.trajin[-4:] == ".pdb":
    		
		pass 
	
	elif args.trajin[-4:] == ".xtc":

		os.system("gmx trjconv -f "+args.trajin+" -s "+args.tprfile+\
        " -pbc nojump -dt "+args.time+" -o "+args.trajin[:-4]+".pdb")

		args.trajin = args.trajin[:-4]+".pdb"

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

def RMSD_HetAtm():

    parser = PDBParser(PERMISSIVE=1)

    # Check the total number of models
    status, models = commands.getstatusoutput("grep -c MODEL " + args.trajin)

    # Setting the progress bar
    pbar2 = tqdm(total = int(models) - 1, desc="Running "+args.option+" calculation")

    # Buffer files
    dataRMSD = []
    dataCenter = []
    dataDistance = []
    distLabels = 'Frame\tCenter\n'
    rmsLabels = 'Frame\tRMSD\n'
    distanceLabels = 'Frame\tDistance\n'


    # Reference structure to RMSD calculation
    refStruct = parser.get_structure("temp/complex_1.pdb", os.getcwd()+"/temp/"+"complex_1.pdb")
    #refAtom = refStruct[0][args.chain][int(args.res)]

    for frame in range(1 ,int(models)):
        
        # Setting target structure
        modStruct = parser.get_structure("temp/complex_" + str(frame) + \
        ".pdb", os.getcwd()+"/temp/"+"complex_" + str(frame) + ".pdb")
        
        if args.option == "rmsd":
                
            # Updating progress bar
            pbar2.update(1)   
            
            # Identifing the specific residue
            refAtom = refStruct[0][args.chain][int(args.res)]
            #print refAtom.get_coord()
            modAtom = modStruct[0][args.chain][int(args.res)]
            sup = Superimposer()
            
            # Superimposing reference and target structures
            sup.set_atoms(list(refAtom),list(modAtom))
            rms_hetatm = sup.rms

            # Saving values
            rmsLabels += "%i\t%.2f\n" %(frame,rms_hetatm)
        
        # Distance between two geometric centers
        elif args.option == "dist":
        
            geometricalRes = args.geometrical.split(',')

            # Updating progress bar
            pbar2.update(1) 

            # Calculate the geometric center of three alpha carbons 
            geoAtom1 = modStruct[0][args.chain][int(geometricalRes[0])]['CA']
            geoAtom2 = modStruct[0][args.chain][int(geometricalRes[1])]['CA']
            geoAtom3 = modStruct[0][args.chain][int(geometricalRes[2])]['CA']

            vector1 = geoAtom1.get_coord()
            vector2 = geoAtom2.get_coord()
            vector3 = geoAtom3.get_coord()

            protAtom = (vector1 + vector2 + vector3)/3

            # Calculate the geometric center of the ligand or residue
            allCoord = []
            for residue in modStruct.get_residues():
                if residue.id[1] == int(args.ligand):
                    for atom in residue:
                        ligVec = (atom.get_coord())
                        allCoord.append(ligVec)
                    ligAtom = sum(allCoord)/len(residue)
                    totDist = protAtom - ligAtom
                    distProtLig = math.sqrt(np.sum(totDist * totDist))

                    # Saving values
                    distLabels += "%i\t%.2f\n" %(frame,distProtLig)

                else:
                    pass

        elif args.option == 'simple_dist':
            
            simpleRes = args.simpleRes.split(',')

            # Updating progress bar
            pbar2.update(1) 

            # Calculate the distance between alpha carbons 
            distanceRes = modStruct[0][args.chain][int(simpleRes[0])]['CA'].coord - modStruct[0][args.chain][int(simpleRes[1])]['CA'].coord
            simpleDistance = np.sqrt(np.sum(distanceRes * distanceRes))
            distanceLabels += "%i\t%.2f\n" % (frame, simpleDistance)

    # Save RMSD file
    if args.option == "rmsd":
		
		dataRMSD.append(rmsLabels)
		RMSDoutput = open(args.output+".tbl", "w")
		RMSDoutput.writelines(dataRMSD)	
		RMSDoutput.close()

    # Save distance of geometric centers
    if args.option == "dist":
    		
		dataCenter.append(distLabels)
		Centeroutput = open(args.output+".tbl", "w")
		Centeroutput.writelines(dataCenter)	
		Centeroutput.close()

    # Save distance of geometric centers
    if args.option == "simple_dist":
    		
		dataDistance.append(distanceLabels)
		Distanceoutput = open(args.output+".tbl", "w")
		Distanceoutput.writelines(dataDistance)	
		Distanceoutput.close()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='''
	
	\t\t\t\t############################################
	\t\t\t\t#                                          #
	\t\t\t\t#         Gambiarra Package Presents       #      
	\t\t\t\t#                                          #
	\t\t\t\t############################################\n\n

    RMSD calculation per residue and calculate the distance between two geometric centers

	''', formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('-m', '--TrajFile\t\t\t', action='store', required=False, dest='trajin', help="The trajectory must be a Multi-PDB or XTC file.\n")
    parser.add_argument('-t', '--analysisType\t\t', action='store', required=False, dest='option', type = str, help='''Specify the type of your analysis: "rmsd", "dist" or "simple_dist". 
    If you want to calculate the "rmsd" make sure that you have specified the flag "-res", and remember that RMSD use the first frame as reference structure.
    If you want to measure the distance between two geometric centers make sure that you have specified "-geo_res" and "-lig" flags.
    Otherwise, if you just need to calculate a simple distance between two alpha carbons make sure that you have specified the flag "-sim_res".\n''')
    parser.add_argument('-o', '--OutputName\t\t', action='store', dest='output', help="The output name without extension.\n")
    parser.add_argument('-c', '--ChainName\t\t', action='store', dest='chain', default=' ', help="Specify the chain.\n")
    parser.add_argument('-res', '--ResidueNumber\t\t', action='store', dest='res', help='''Specify the residue number to monitor during the "rmsd".\n''')
    parser.add_argument('-geo_res', '--GeometricResidues\t\t', action='store', dest='geometrical', help="Specify three residues number to calculate the geometric center.\n")
    parser.add_argument('-lig', '--LigandResidues\t\t', action='store', dest='ligand', help="Specify one residue or ligand number to calculate the geometric center.\n")
    parser.add_argument('-sim_res', '--SimpleResidues\t\t', action='store', dest='simpleRes', help="Specify residues numbers to monitor the simple distance.\n")
    args = parser.parse_args()
    ConvertingXTCFileToPDBFile()
    TrajectorySplitter()
    RMSD_HetAtm()