//traj.cpp

/*********************************************************************/
/* This source code file is part of OOCCuPy software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/*********************************************************************/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <sstream>

#include "../include/traj.h"
#include "../include/global.h"

using std::string;
using std::vector;

/***********************************************************************/
traj::traj(){
	
}
/***********************************************************************/
traj::traj(string file_name):
	traj_file(file_name)		{
		
	py_name = remove_extension( file_name.c_str() );
	py_name += "_mdtraj.py";
			
	python_script.open( py_name.c_str() );
	
	python_script	<< " #/usr/bin/env python \n"
					<< " # -*- coding: utf-8 -*-\n\n"
					<< "import mdtraj as md \n"
					<< "import os\n";
}
/***********************************************************************/
traj::traj(string file_name, vector<int> atoms):
	traj_file(file_name)					   ,
	atoms_pairs(atoms)						   {
		
	py_name = remove_extension( file_name.c_str() );
	py_name += "_mdtraj.py";
			
	python_script.open( py_name.c_str() );
	
	python_script	<< " #/usr/bin/env python \n"
					<< " # -*- coding: utf-8 -*- \n\n"
					<< "import mdtraj as md \n"
					<< "import os";

	
}
/***********************************************************************/
traj::~traj(){
}
/***********************************************************************/
void traj::mdtraj_geo(){
	string topname = change_extension( traj_file.c_str(),".pdb" );
	
	python_script << "trj = md.load('" << traj_file << "',top='" << topname << "')\n"
				  << "rg = md.rmsd(trj,trj)\n"
				  << "rmsd = md.compute_rg(trj)\n"
				  << "file_txt = 'Time RMSD RG'\n"
				  << "file_rmsd = open('md_geoD.txt','w')\n"
				  << "for i in range(len(rmsd)):\n"
				  << "	file_txt += str(trj.time[i]) +' '+ str(rmsd[i]) +' '+str(rg[i]) +os.linesep\n"
				  << "file_rmsd.write(file_txt)\n"
				  << "file_rmsd.close()\n";
				  
	
	python_script.close();
	
	string comand = "python3 " + py_name; 
	system( comand.c_str() );	

	vector<double> time;
	vector<double> rmsd;
	vector<double> rg;
	int line = 0;
	char tmp_line[50];
	double tmp_doub;
	int col = 0;
	string dat_name = "md_geoD.txt";
	std::ifstream geo_data( dat_name.c_str() );
	while( !geo_data.eof() ){
		geo_data.getline(tmp_line,50);
		if ( line > 0 ){
			std::stringstream stream(tmp_line);
			while( stream >> tmp_doub ){
				if ( col == 0){
					time.push_back(tmp_doub);
					col++;
				}else if ( col == 1){
					rmsd.push_back(tmp_doub);
					col++;
				}else if ( col == 2 ){
					rg.push_back(tmp_doub);
					col = 0;
				}
			}			
		}
		line++;
	}		
}
/***********************************************************************/
void traj::mdtraj_distances(){
	
}
/***********************************************************************/
void traj::extract_frame(const char* pdb_file, int* frames){
	
}
/***********************************************************************/
std::vector<int> traj::bi_most_probable_point( std::vector<double> v1, std::vector<double> v2 ){
	
}

/***********************************************************************/
