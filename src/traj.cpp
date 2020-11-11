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
#include <vector<
#include <string>
#include <fstream>
#include <include/traj.h>
#include <global.h>

using std::string;

/***********************************************************************/
traj::traj(){
	
}
/***********************************************************************/
traj::traj(std::string file_name):
	traj_file(file_name)		{
		
	std::string py_name =  file_name
		
	python_script.open( file_name.c_str() );
}
/***********************************************************************/
traj::traj(std::string file_name, std<vector>int atoms)	:
	traj_file(file_name)								,
	atoms_pairs(atoms)									{
	python_script.open( file_name.c_str() );
}
/***********************************************************************/
traj::~traj(){
	python_script.close();
}
/***********************************************************************/
void traj::mdtraj_geo(){
	

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
