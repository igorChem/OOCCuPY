//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>

//Including PRIMoRDiA headers
#include "log_class.h"
#include "common.h"
#include "Iaorbital.h"
#include "Iatom.h"
#include "Imolecule.h"
#include "Iline.h"
#include "Ibuffer.h"
#include "QMparser.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>

// Aliases for standard c++ scope functions
using std::move;
using std::vector;
using std::string;
using std::unique_ptr;
using std::cout;
using std::endl;
using std::stoi;
using std::stod;

/************************************************************************************/
QMparser::QMparser()            :
	buffer( new Ibuffer() )     ,
	name_f("nonamed")           ,
	molecule( new Imolecule() ) ,
	parsed(false)               {
}
/************************************************************************************/
QMparser::QMparser(const QMparser& rhs_QMp)     :
	name_f(rhs_QMp.name_f)                      ,
	molecule( new Imolecule(*rhs_QMp.molecule) ),
	buffer( new Ibuffer(*rhs_QMp.buffer) )      ,
	program(rhs_QMp.program)                    {
}
/************************************************************************************/
QMparser& QMparser::operator=(const QMparser& rhs_QMp) {
	if( this!=&rhs_QMp ){
		name_f    = rhs_QMp.name_f;
		*molecule = *rhs_QMp.molecule;
		*buffer   = *rhs_QMp.buffer;
		program   = rhs_QMp.program;
	}
	return *this;
}
/************************************************************************************/
QMparser::QMparser(QMparser&& rhs_QMp) noexcept:
	name_f( move(rhs_QMp.name_f) )             ,
	molecule( move(rhs_QMp.molecule) )         ,
	buffer( move(rhs_QMp.buffer) )             ,
	program( move(rhs_QMp.program) )           { 
}
/************************************************************************************/
QMparser& QMparser::operator=(QMparser&& rhs_QMp) noexcept{
	if( this!=&rhs_QMp ){
		name_f   = move(rhs_QMp.name_f);
		molecule = move(rhs_QMp.molecule);
		buffer   = move(rhs_QMp.buffer);
		program  = move(rhs_QMp.program);
	}
	return *this;
}
/************************************************************************************/
bool QMparser::parse_mopac_aux( const char* file_name){
	
	if ( !check_file_ext(".aux",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	
	m_log->input_message("Starting to parse AUX file from MOPAC.");
	
	unsigned int i,j;
	vector<string> key_wrds_in;
	vector<string> key_wrds_fin;
	vector<double> zetas;
	vector<int> shells;
	vector<int> aoidx;
	
	name_f = file_name;
	molecule->name = remove_extension(file_name);

	buffer.reset( new Ibuffer(file_name,"ATOM_EL[","ATOM_CORE[") );
	for(i=0; i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++){
			Iatom atom;
			atom.set_type( buffer->lines[i].pop_string(0) );
			molecule->add_atom(atom);
		}
	}
	
	if ( molecule->atoms.size() == 0 ){ cout << "Warning! Zero atoms read in aux file. Verify your file!!" << endl; } 

	m_log->input_message("Found number of atoms in the aux file: ");
	m_log->input_message( int(molecule->num_of_atoms) );

	int elecN = 0;
	buffer.reset( new Ibuffer(file_name,"ATOM_CORE[","ATOM_X:ANGSTROMS[") );
	for(i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) elecN += buffer->lines[i].pop_int(0); 
	}
	
	molecule->num_of_electrons = elecN;
	
	m_log->input_message("Number of electron in the system: ");
	m_log->input_message( elecN );
	
	int counter = 0;
	buffer.reset ( new Ibuffer(file_name,"ATOM_X:ANGSTROMS[","AO_ATOMINDEX[") );
	for (i= 0;i<buffer->nLines;i++ ){
			molecule->atoms[counter].xcoord = buffer->lines[i].pop_double(0);
			molecule->atoms[counter].xcoord = buffer->lines[i].pop_double(0);
			molecule->atoms[counter++].xcoord = buffer->lines[i].pop_double(0);	
	}
	
	buffer.reset( new Ibuffer(file_name,"AO_ATOMINDEX[","ATOM_SYMTYPE[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) aoidx.push_back( buffer->lines[i].pop_int(0) );
	}	
	
	m_log->input_message("Number of atomic orbitals: ");
	m_log->input_message( int( aoidx.size() ) );
	
	counter = 0;
	buffer.reset ( new Ibuffer(file_name,"ATOM_SYMTYPE[","AO_ZETA[") );
	for (i=0; i<buffer->nLines;i++ ){
		for (j=0; j<buffer->lines[i].line_len;j++){
			Iaorbital aorb;
			aorb.symmetry = buffer->lines[i].pop_string(0);
			if 		( aorb.symmetry == "PX" ) aorb.powx = 1;
			else if ( aorb.symmetry == "PY" ) aorb.powy = 1;
			else if ( aorb.symmetry == "PZ" ) aorb.powz = 1;
			molecule->atoms[aoidx[counter++]-1].add_orbital(aorb);
		}
	}

	buffer.reset ( new Ibuffer(file_name,"AO_ZETA[","ATOM_PQN[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) zetas.push_back( buffer->lines[i].pop_double(0) );		
	}
	
	buffer.reset ( new Ibuffer(file_name,"ATOM_PQN[","NUM_ELECTRONS=") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) shells.push_back( buffer->lines[i].pop_int(0) );		
	}
	
	counter = 0;
	buffer.reset ( new Ibuffer(file_name,"ATOM_CHARGES[","OVERLAP_MATRIX[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->atoms[counter++].charge = buffer->lines[i].pop_double(0) ;
	}
	
	key_wrds_in.push_back("OVERLAP_MATRIX[");
	key_wrds_fin.push_back("SET_OF_MOS=");
	key_wrds_fin.push_back("SET_OF_ALPHA_MOS=");
	buffer.reset( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=1; i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->m_overlap.push_back( buffer->lines[i].pop_double(0)); 
	}
	
	key_wrds_in[0] = "EIGENVECTORS[";
	key_wrds_in.push_back("ALPHA_EIGENVECTORS[");
	key_wrds_in.push_back("LMO_VECTORS[");
	key_wrds_fin[0] = "TOTAL_DENSITY_MATRIX[";
	key_wrds_fin[1] = "SET_OF_BETA_MOS=";
	key_wrds_fin.push_back("DENSITY_MATRIX");
	
	buffer.reset ( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->coeff_MO.push_back( buffer->lines[i].pop_double(0) );
	}

	buffer.reset( new Ibuffer(file_name,"BETA_EIGENVECTORS[","ALPHA_DENSITY_MATRIX[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->coeff_MO_beta.push_back( buffer->lines[i].pop_double(0) );
	}
	
	if ( molecule->coeff_MO_beta.size() > 0 ) { molecule->betad = true; }
	
	key_wrds_in.resize(3);
	key_wrds_in[0] = "EIGENVALUES[";
	key_wrds_in[1] = "ALPHA_EIGENVALUES[";
	key_wrds_in[2]= "LMO_ENERGY_LEVELS[";
	key_wrds_fin.resize(2);
	key_wrds_fin[0] = "MOLECULAR_ORBITAL_OCCUPANCIES[";
	key_wrds_fin[1] = "BETA_EIGENVALUES[";
	
	buffer.reset( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++){
			molecule->orb_energies.push_back( buffer->lines[i].pop_double(0) );
			molecule->MOnmb++;
		}
	}

	buffer.reset( new Ibuffer(file_name,"BETA_EIGENVALUES[","ALPHA_MOLECULAR_ORBITAL_OCCUPANCIES[") );
	for (i=0;i<buffer->nLines;i++){
		for (unsigned int j=0;j<buffer->lines[i].line_len;j++){
			molecule->orb_energies_beta.push_back( buffer->lines[i].pop_double(0) );
			molecule->MOnmb_beta++;
		}
	}
	
	key_wrds_in[0] = "ALPHA_MOLECULAR_ORBITAL_OCCUPANCIES";
	key_wrds_in[1] = "MOLECULAR_ORBITAL_OCCUPANCIES[";
	key_wrds_in.erase(key_wrds_in.begin()+2,key_wrds_in.end());
	key_wrds_fin[0] = "BETA_MOLECULAR_ORBITAL_OCCUPANCIES[" ;
	key_wrds_fin[1] =  "CPU_TIME:SECONDS[1]=";
	
	buffer.reset( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->occupied.push_back( buffer->lines[i].pop_int(0) );
	}	
		
	buffer.reset( new Ibuffer(file_name,"BETA_MOLECULAR_ORBITAL_OCCUPANCIES[","CPU_TIME:SECONDS[1]=") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->occupied_beta.push_back( buffer->lines[i].pop_int(0) );
	}
	
	buffer.reset(nullptr);
	
	counter = 0;
	for (i=0;i<molecule->atoms.size();i++){
		for (j=0;j<molecule->atoms[i].norb;j++){
			molecule->atoms[i].orbitals[j].alpha = zetas[counter];
			molecule->atoms[i].orbitals[j].shell = shells[counter++];
		}
	}
	
	molecule->norm_orbs();
	molecule->get_homo_energy();
	molecule->get_lumo_energy();
	molecule->num_of_electrons=0;
	
	for(i=0;i<molecule->occupied.size();i++) molecule->num_of_electrons += molecule->occupied[i]; 
	for(i=0;i<molecule->occupied_beta.size();i++){molecule->num_of_electrons += molecule->occupied_beta[i]; }
	
	m_log->input_message("HOMO energy: ");
	m_log->input_message(double_to_string(molecule->homo_energy));
	m_log->input_message("LUMO energy: ");
	m_log->input_message(double_to_string(molecule->lumo_energy));
	
	parsed = true;
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_out_mopac(const char* file_name){
	
	if ( !check_file_ext(".out",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	
	m_log->input_message("Starting to parse out file from MOPAC.");
	
	name_f = file_name;
	molecule->name = remove_extension(file_name);
	
	buffer.reset( new Ibuffer(file_name,true) );
	for (int i=0;i<buffer->nLines;i++){
		if ( buffer->lines[i].IF_line("HEAT",1,"FORMATION",3,10) ){ 
			 molecule->heat_of_formation = buffer->lines[i].pop_double(5);
		}
		else if ( buffer->lines[i].IF_line("ELECTRONIC",0,"ENERGY",1,5) || buffer->lines[i].IF_line("ELECTRONIC",0,"ENERGY",1,8) ){
			molecule->energy_tot = buffer->lines[i].pop_double(3) ;
		}
		else if ( buffer->lines[i].IF_line("HOMO",0,"LUMO",1,7) ){
			molecule->homo_energy = buffer->lines[i].pop_double(5);
			molecule->lumo_energy = buffer->lines[i].pop_double(5);
		}
		else if ( buffer->lines[i].IF_line("SUM",0,5) ){
			for ( int j=0;j<3;j++){ molecule->dipole_moment[j] = buffer->lines[i].pop_double(1); }
			molecule->total_dipmoment = buffer->lines[i].pop_double(1);
		}
	}

	m_log->input_message("Total Energy ");
	m_log->input_message(double_to_string(molecule->energy_tot));
	m_log->input_message("HOMO energy ");
	m_log->input_message(double_to_string(molecule->homo_energy));
	m_log->input_message("LUMO energy ");
	m_log->input_message(double_to_string(molecule->lumo_energy));
	
	parsed = true;
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_mozyme(const char* file_name){
	
	if ( !check_file_ext(".aux",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return parsed;
	}
	
	m_log->input_message("Starting to parse AUX file from MOPAC with MOZYME.");
	unsigned int i,j;
	vector<string> key_wrds_in;
	vector<string> key_wrds_fin;
	vector<double> zetas;
	vector<int> shells;
	vector<int> aoidx;
	int elecN = 0;
	
	name_f = file_name;
	molecule->name = remove_extension(file_name);
	
	buffer.reset( new Ibuffer(file_name,"ATOM_EL[","ATOM_CORE[") );
	for(i=0;i<buffer->nLines;i++){
		for(j=0;j<buffer->lines[i].line_len;j++){
			Iatom atom;
			atom.set_type( buffer->lines[i].pop_string(0) );
			molecule->add_atom(atom);
		}
	}
	
	m_log->input_message("Found number of atoms in the aux file: ");
	m_log->input_message( int(molecule->num_of_atoms) );
	
	buffer.reset( new Ibuffer(file_name,"ATOM_CORE[","ATOM_X:ANGSTROMS[") );
	for(i=0;i<buffer->nLines;i++){
		for (j= 0;j<buffer->lines[i].line_len;j++) elecN += buffer->lines[i].pop_int(0); 
	}
	
	molecule->num_of_electrons = elecN;
	m_log->input_message("Number of electrons: ");
	m_log->input_message( elecN );
	
	int counter = 0;
	buffer.reset ( new Ibuffer(file_name,"ATOM_X:ANGSTROMS[","AO_ATOMINDEX[") );
	for (i=0;i<buffer->nLines;i++ ){
		molecule->atoms[counter].xcoord   = buffer->lines[i].get_double(0);
		molecule->atoms[counter].xcoord   = buffer->lines[i].get_double(1);
		molecule->atoms[counter++].xcoord = buffer->lines[i].get_double(2);
	}
	
	buffer.reset( new Ibuffer(file_name,"AO_ATOMINDEX[","ATOM_SYMTYPE[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) aoidx.push_back( buffer->lines[i].pop_int(0) );
	}	
	
	m_log->input_message("Number of atomic orbitals: ");
	m_log->input_message( int_to_string( aoidx.size() ) );
	
	counter = 0;
	buffer.reset ( new Ibuffer(file_name,"ATOM_SYMTYPE[","AO_ZETA[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++){
			Iaorbital aorb;
			aorb.symmetry = std::move( buffer->lines[i].words[j] );
			if 		( aorb.symmetry == "PX" ) aorb.powx = 1;
			else if ( aorb.symmetry == "PY" ) aorb.powy = 1;
			else if ( aorb.symmetry == "PZ" ) aorb.powz = 1;
			molecule->atoms[aoidx[counter++]-1].add_orbital(aorb);
		}
	}

	buffer.reset ( new Ibuffer(file_name,"AO_ZETA[","ATOM_PQN[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) zetas.push_back( buffer->lines[i].pop_double(0) );
	}
	
	buffer.reset ( new Ibuffer(file_name,"ATOM_PQN[","NUM_ELECTRONS=") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) shells.push_back( buffer->lines[i].pop_int(0) );
	}
	
	counter = 0;
	buffer.reset ( new Ibuffer(file_name,"ATOM_CHARGES[","OVERLAP_MATRIX[") );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->atoms[counter++].charge = buffer->lines[i].pop_double(0);
	}
	
	buffer.reset( new Ibuffer(file_name,"OVERLAP_MATRIX[","SET_OF_MOS=") );
	for (i=1;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->m_overlap.push_back( buffer->lines[i].pop_double(0) );
	}

	key_wrds_in.push_back("EIGENVECTORS[");
	key_wrds_in.push_back("LMO_VECTORS[");
	key_wrds_fin.push_back("DENSITY_MATRIX[");
	
	buffer.reset ( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->coeff_MO.push_back( buffer->lines[i].pop_double(0) );
	}

	key_wrds_in[0]  = "EINGENVALUES[";
	key_wrds_in[1]  = "LMO_ENERGY_LEVELS[";
	key_wrds_fin[0] = "MOLECULAR_ORBITAL_OCCUPANCIES[";
	molecule->MOnmb = 0;
	buffer.reset( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++){
			molecule->orb_energies.push_back( buffer->lines[i].pop_double(0) );
			molecule->MOnmb++;
		}
	}
	
	key_wrds_in[0] = "ALPHA_MOLECULAR_ORBITAL_OCCUPANCIES";
	key_wrds_in[1] = "MOLECULAR_ORBITAL_OCCUPANCIES[";
	key_wrds_in.erase(key_wrds_in.begin()+2,key_wrds_in.end());
	key_wrds_fin[0] = "BETA_MOLECULAR_ORBITAL_OCCUPANCIES[" ;
	key_wrds_fin.emplace_back("CPU_TIME:SECONDS[1]=");
	
	buffer.reset( new Ibuffer(file_name,key_wrds_in,key_wrds_fin) );
	for (i=0;i<buffer->nLines;i++){
		for (j=0;j<buffer->lines[i].line_len;j++){ molecule->occupied.push_back( buffer->lines[i].pop_int(0) );	}
	}	

	buffer.reset(nullptr);
	
	counter = 0;
	for (i=0;i<molecule->atoms.size();i++){
		for (j=0;j<molecule->atoms[j].norb;j++){
			molecule->atoms[i].orbitals[j].alpha = zetas[counter];
			molecule->atoms[i].orbitals[j].shell = shells[counter++];
		}
	}
	
	molecule->norm_orbs();
	molecule->get_homo_energy();
	molecule->get_lumo_energy();
	
	m_log->input_message("HOMO energy: ");
	m_log->input_message(double_to_string(molecule->homo_energy));
	m_log->input_message("LUMO energy: ");
	m_log->input_message(double_to_string(molecule->lumo_energy));
	
	parsed = true;
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_mopac_mgf(const char* file_name){
	
	if ( !check_file_ext(".mgf",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	
	name_f = remove_extension(file_name);
	molecule->name = name_f;
	
	vector<double> zetasS;
	vector<double> zetasP;
	vector<double> zetasD;
	vector<double> inv_mat;
	vector<int>    orbN;
	vector<int>    orbN_beta;
	
	unsigned int i,j = 0;
	int inmat_in      = 0;
	int inmat_fin     = 0;
	int noa           = 0;
	double tot_charge = 0;
 	
	m_log->input_message("Opening mgf file");
	
	buffer.reset( new Ibuffer(file_name,true) );
		
	int counter = 0;
	for(i=0;i<buffer->nLines;i++){
		if( i == 0 ) noa = buffer->lines[i].pop_int(0);
		else if ( i>0 && i<=noa ) {
			string atypr = get_atomic_symbol( buffer->lines[i].pop_int(0) );
			double xcrd  = buffer->lines[i].pop_double(0);
			double ycrd  = buffer->lines[i].pop_double(0);
			double zcrd  = buffer->lines[i].pop_double(0);
			tot_charge  += buffer->lines[i].pop_double(0);
			molecule->add_atom(xcrd,ycrd,zcrd,atypr);
		}
		else if( i>noa && i<=(noa*2)){
			zetasS.push_back( buffer->lines[i].pop_double(0) );
			zetasP.push_back( buffer->lines[i].pop_double(0) );
			zetasD.push_back( buffer->lines[i].pop_double(0) );
		}
		else if( buffer->lines[i].IF_word("ORBITAL",0,7) && inmat_fin ==0 ) orbN.push_back(i);
		else if( buffer->lines[i].IF_word("INVERSE_MATRIX[",0,15) ) inmat_in = i;
		else if( buffer->lines[i].IF_word("Keywords:",0,9)) inmat_fin = i;
		else if( inmat_fin > 0 && buffer->lines[i].IF_word("ORBITAL",0,7) ){
				molecule->betad = true;
				orbN_beta.push_back(i);
			}
		}
		
	string temp = "noname";
	for(int i=0;i<orbN.size();i++){
		int fin_ind = 0;
		if ( i==orbN.size()-1 )	fin_ind = inmat_in;
		else fin_ind = orbN[i+1];
		for(int j=orbN[i];j<fin_ind;j++){
			if ( j == orbN[i] ){
				molecule->occupied.push_back( buffer->lines[j].pop_int(1) );
				molecule->orb_energies.push_back( buffer->lines[j].pop_double(2) );
				molecule->MOnmb++;
			}else{ 
				for(int k=0;k<buffer->lines[j].line_len;k++){
					if (k == 0 &&  buffer->lines[j].words[k][0] == '-' ){
						if ( buffer->lines[j].words[k].size() > 15 ){
							temp = buffer->lines[j].words[k].substr(0,15);
							molecule->coeff_MO.push_back( D_E_conv(temp) );
							temp = buffer->lines[j].words[k].substr( 15,buffer->lines[j].words[k].size() );
							if ( temp.size() > 15) {
								temp = temp.substr(0,15);
								molecule->coeff_MO.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 30,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 45,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15 ){
										temp = temp.substr(0,15);
										molecule->coeff_MO.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr(60,buffer->lines[j].words[k].size() );
										molecule->coeff_MO.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO.push_back( D_E_conv(temp) );
						}else molecule->coeff_MO.push_back( D_E_conv(buffer->lines[j].words[k]) );
					}else{
						if ( buffer->lines[j].words[k].size() > 14 ){
							temp = buffer->lines[j].words[k].substr(0,14);
							molecule->coeff_MO.push_back( D_E_conv(temp) );
							temp = buffer->lines[j].words[k].substr( 14,buffer->lines[j].words[k].size() );
							if ( temp.size() > 15) {
								temp = temp.substr(0,15);
								molecule->coeff_MO.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 29,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 44,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15 ){
										temp = temp.substr(0,15);
										molecule->coeff_MO.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr(59,buffer->lines[j].words[k].size() );
										molecule->coeff_MO.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO.push_back( D_E_conv(temp) );
						}else molecule->coeff_MO.push_back( D_E_conv(buffer->lines[j].words[k]) );						
					}	
				}
			}
		}
	}
	for(int i=0;i<noa;i++){
		int sh = 0;
		Iaorbital orbS;
		orbS.alpha = zetasS[i];
		if 		( molecule->atoms[i].atomicN > 1  && molecule->atoms[i].atomicN < 11 ) sh = orbS.shell = 2;
		else if ( molecule->atoms[i].atomicN > 10 && molecule->atoms[i].atomicN < 19 ) sh = orbS.shell = 3;
		else if ( molecule->atoms[i].atomicN > 18 && molecule->atoms[i].atomicN < 37 ) sh = orbS.shell = 4;
		molecule->atoms[i].add_orbital(orbS);
		if ( zetasP[i] > 1e-05 ){
			Iaorbital orbPx, orbPy, orbPz;
			orbPx.alpha = orbPy.alpha = orbPz.alpha = zetasP[i];
			orbPx.shell = orbPy.shell = orbPz.shell = sh;
			orbPx.symmetry = "PX";
			orbPy.symmetry = "PY";
			orbPz.symmetry = "PZ";
			orbPx.powx  = 1;  
			orbPy.powy  = 1;  
			orbPz.powz  = 1;
			molecule->atoms[i].add_orbital(orbPx);
			molecule->atoms[i].add_orbital(orbPy);
			molecule->atoms[i].add_orbital(orbPz);
		}
		if ( zetasD[i] > 1e-05 ){
			Iaorbital orbDxx, orbDyy, orbDzz, orbDxy, orbDyz, orbDxz;
			orbDxx.alpha = orbDyy.alpha = orbDzz.alpha = orbDxy.alpha = orbDyz.alpha = orbDxz.alpha = zetasD[i];
			orbDxx.shell = orbDyy.shell = orbDzz.shell = orbDxy.shell = orbDyz.shell = orbDxz.shell = sh;
			orbDxx.powx  = 2;
			orbDxx.symmetry  = "XX";
			orbDyy.powy  = 2;
			orbDyy.symmetry  = "YY";
			orbDzz.powz  = 2;
			orbDzz.symmetry  = "zz";
			orbDxy.powx  = orbDxy.powy  = 1;
			orbDxy.symmetry  = "XY";
			orbDyz.powy  = orbDyz.powz  = 1;
			orbDyz.symmetry  = "YZ";
			orbDxz.powx  = orbDxz.powz  = 1;
			orbDxz.symmetry  = "XZ";
			molecule->atoms[i].add_orbital(orbDxx);
			molecule->atoms[i].add_orbital(orbDyy);
			molecule->atoms[i].add_orbital(orbDzz);
			molecule->atoms[i].add_orbital(orbDxy);
			molecule->atoms[i].add_orbital(orbDyz);
			molecule->atoms[i].add_orbital(orbDxz);
		}
	}
	
	for(int i=inmat_in+1;i<inmat_fin;i++){
		for(int k=0;k<buffer->lines[i].line_len;k++){
			temp = buffer->lines[i].words[k];
			if ( k == 0 && temp[0] == '-' ){
				if ( temp.size() > 15 ){ 
					temp = temp.substr(0,15);
					inv_mat.push_back( D_E_conv(temp) );
					temp = buffer->lines[i].words[k].substr( 15,buffer->lines[i].words[k].size() );
					if ( temp.size() > 15 ){
						temp = temp.substr(0,15);
						inv_mat.push_back( D_E_conv(temp) );
						temp = buffer->lines[i].words[k].substr( 30,buffer->lines[i].words[k].size() );
						if ( temp.size() > 15) {
							temp = temp.substr(0,15);
							inv_mat.push_back( D_E_conv(temp) );
							temp = buffer->lines[i].words[k].substr( 45,buffer->lines[i].words[k].size() );
							if ( temp.size() > 15 ){
								temp = temp.substr(0,15);
								inv_mat.push_back( D_E_conv(temp) );
								temp = buffer->lines[i].words[k].substr( 60,buffer->lines[i].words[k].size() );
								inv_mat.push_back( D_E_conv(temp) );
							}else inv_mat.push_back( D_E_conv(temp) );
						}else inv_mat.push_back( D_E_conv(temp) );
					}else inv_mat.push_back( D_E_conv(temp) );
				}else inv_mat.push_back( D_E_conv(temp) );
			}else{
				if( temp.size() > 14 ){
					temp = temp.substr(0,14);
					inv_mat.push_back( D_E_conv(temp) );
					temp = buffer->lines[i].words[k].substr( 14,buffer->lines[i].words[k].size() );
					if ( temp.size() > 15 ){
						temp = temp.substr(0,15);
						inv_mat.push_back( D_E_conv(temp) );
						temp = buffer->lines[i].words[k].substr( 29,buffer->lines[i].words[k].size() );
						if ( temp.size() > 15 ) {
							temp = temp.substr(0,15);
							inv_mat.push_back( D_E_conv(temp) );
							temp = buffer->lines[i].words[k].substr( 44,buffer->lines[i].words[k].size() );
							if ( temp.size() > 15 ){
								temp = temp.substr(0,15);
								inv_mat.push_back( D_E_conv(temp) );
								temp = buffer->lines[i].words[k].substr( 59,buffer->lines[i].words[k].size() );
								inv_mat.push_back( D_E_conv(temp) );
							}else inv_mat.push_back( D_E_conv(temp) );
						}else inv_mat.push_back( D_E_conv(temp) );
					}else inv_mat.push_back( D_E_conv(temp) );
				}else inv_mat.push_back( D_E_conv(temp) ); 
			}
		}
	}	
	int nmo = molecule->MOnmb;
	int k   = 0;
	
	
	for(int i=0;i<orbN_beta.size();i++){		
		int fin_ind = 0;
		if ( i==orbN_beta.size()-1 ) fin_ind = buffer->nLines;
		else fin_ind = orbN_beta[i+1];		
		for(int j=orbN_beta[i];j<fin_ind;j++){
			if ( j == orbN_beta[i] ){
				molecule->occupied_beta.push_back( buffer->lines[j].pop_int(1) );
				molecule->orb_energies_beta.push_back( buffer->lines[j].pop_double(2) );
				molecule->MOnmb_beta++;
			}else{
				for(int k=0;k<buffer->lines[j].line_len;k++){
					if (k == 0 &&  buffer->lines[j].words[k][0] == '-' ){
						if ( buffer->lines[j].words[k].size() > 15 ){
							temp = buffer->lines[j].words[k].substr(0,15);
							molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
							temp = buffer->lines[j].words[k].substr( 15,buffer->lines[j].words[k].size() );
							if ( temp.size() > 15) {
								temp = temp.substr(0,15);
								molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 30,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 45,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15 ){
										temp = temp.substr(0,15);
										molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr(60,buffer->lines[j].words[k].size() );
										molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
						}else molecule->coeff_MO_beta.push_back( D_E_conv(buffer->lines[j].words[k]) );
					}else{
						if ( buffer->lines[j].words[k].size() > 14 ){
							temp = buffer->lines[j].words[k].substr(0,14);
							molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
							temp = buffer->lines[j].words[k].substr( 14,buffer->lines[j].words[k].size() );
							if ( temp.size() > 15) {
								temp = temp.substr(0,15);
								molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 29,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 44,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15 ){
										temp = temp.substr(0,15);
										molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr(59,buffer->lines[j].words[k].size() );
										molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
						}else molecule->coeff_MO_beta.push_back( D_E_conv(buffer->lines[j].words[k]) );						
					}	
				}
			}
		}
	}
	buffer.reset();
	
	Eigen::MatrixXd inv(nmo,nmo);
	Eigen::MatrixXd coeff_corrected(nmo,nmo);

	int bnmo = molecule->MOnmb_beta;
 
	Eigen::MatrixXd coeff_corrected_beta(bnmo,bnmo);
	
	for (i=0;i<nmo;i++){
		for (j=0;j<nmo;j++){  coeff_corrected(i,j) = molecule->coeff_MO[i*nmo+j]; }
	}
	
	for (i=0;i<bnmo;i++){
		for (j=0;j<bnmo;j++) coeff_corrected_beta(i,j) = molecule->coeff_MO_beta[i*nmo+j];
	}
	
	for(i=0;i<nmo;i++){
		for(j=0;j<=i;j++){
			if (i==j)  inv(i,j) = inv_mat[k++];
			else       inv(j,i) = inv(i,j) = inv_mat[k++];
		}
	}
	
	coeff_corrected = coeff_corrected*inv;
	for (int i=0;i<nmo;i++){
		for (int j=0;j<nmo;j++)  molecule->coeff_MO[i*nmo+j] = coeff_corrected(i,j);
	}
	
	if ( molecule->betad ){
		coeff_corrected_beta = coeff_corrected_beta*inv; 
		for (i=0;i<nmo;i++){
			for (j=0;j<nmo;j++)	molecule->coeff_MO_beta[i*nmo+j] = coeff_corrected_beta(i,j);
		}
	}
	
	molecule->get_homo_energy();
	molecule->get_lumo_energy();
	
	m_log->input_message("HOMO energy: ");
	m_log->input_message(double_to_string(molecule->homo_energy));
	m_log->input_message("LUMO energy: ");
	m_log->input_message(double_to_string(molecule->lumo_energy));
	
	molecule->norm_orbs();
	parsed = true;
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_log_gamess(const char* file_name){
	
	if ( !check_file_ext(".log",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	
	m_log->input_message("Starting to parse log file from GAMESS.");

	unsigned int i,j,k,l;

	int atom_in   = 0;
	int atom_fin  = 0;
	int fmo_in    = 0;
	int fmo_fin   = 0;
	int fmob_in   = 0;
	int fmob_fin  = 0;
	int basis_in  = 0;
	int basis_out = 0;
	int chg_in    = 0;
	int chg_fin   = 0;
	int noe       = 0;
	
	name_f = file_name;
	molecule->name = remove_extension(file_name);

	buffer.reset( new Ibuffer(file_name,true) );
	for(i=0;i<buffer->nLines;i++){
		if      ( buffer->lines[i].IF_line("COORDINATES",2,"(BOHR)",3,4)) atom_in = i;  
		else if ( buffer->lines[i].IF_line("INTERNUCLEAR",0,"(ANGS.)",2,3)) atom_fin = i;
		else if ( buffer->lines[i].IF_line("ATOMIC",0,"BASIS",1,3)) {
			if ( atom_fin == 0 ) atom_fin = i;
		}  
		else if ( buffer->lines[i].IF_line("SHELL",0,"COEFFICIENT(S)",5,6) ) basis_in = i;
		else if ( buffer->lines[i].IF_line("BASIS",3,"SHELLS",5,8) ) basis_out = i;
		else if ( buffer->lines[i].IF_line("EIGENVECTORS",0,1) )  {
			if ( fmo_in == 0 ) fmo_in  = i;
		}
		else if ( buffer->lines[i].IF_line("BETA",1,"SET",2,4) ) fmob_in = fmo_fin = i;
		else if ( buffer->lines[i].IF_line("END",1,"CALCULATION",4,6) ) { 
			if ( fmo_fin > 0 ) fmob_fin = i; 
			else fmo_fin = i;
		}
		else if ( buffer->lines[i].IF_line("NUMBER",0,"ELECTRONS",2,5)) noe = buffer->lines[i].pop_int(4);
		else if ( buffer->lines[i].IF_line("MULLIKEN",1,"LOWDIN",3,6) ) chg_in = i;
		else if ( buffer->lines[i].IF_line("BOND",0,"ANALYSIS",4,8) ) chg_fin = i;
		else if ( buffer->lines[i].IF_line("SOLVENT",4,"A.U.",7,8) ) {
			molecule->energy_tot = stod(buffer->lines[i].words[6]);
		}
		else if ( buffer->lines[i].IF_line("TOTAL",0,"ENERGY",1,4) ) {
			molecule->energy_tot = buffer->lines[i].pop_double(3);
		}
	}
	
	buffer.reset( new Ibuffer(file_name,atom_in,atom_fin));
	for(j=1;j<buffer->nLines;j++){
		if ( buffer->lines[j].line_len == 5 ) {
			string symb = buffer->lines[j].words[0];
			double xcrd = buffer->lines[j].pop_double(2); 
			double ycrd = buffer->lines[j].pop_double(2); 
			double zcrd = buffer->lines[j].pop_double(2); 
			molecule->add_atom(xcrd,ycrd,zcrd,symb);
		}
	}
	m_log->input_message("Found number of atoms in the aux file: ");
	m_log->input_message( int(molecule->num_of_atoms) );
	
	vector<int> atom_n_basis;
	vector<int> shell_n;
	vector<int> shell_n_size;
	vector<string> shell_t;
	vector<double> exponents;
	vector<double> c_coefficients;
	vector<double> cP_coefficients;
	
	int jj = -1;
	buffer.reset( new Ibuffer(file_name, basis_in,basis_out)  );
	for(i=1;i<buffer->nLines;i++){
		if ( buffer->lines[i].line_len == 5 || buffer->lines[i].line_len == 6 ){
			atom_n_basis.push_back(jj);
			shell_n.push_back( buffer->lines[i].pop_int(0) );
			shell_t.push_back( buffer->lines[i].words[0] );
			exponents.push_back( buffer->lines[i].pop_double(2) );
			if ( buffer->lines[i].line_len == 5 ) 	c_coefficients.push_back( buffer->lines[i].pop_double(2) );			
			else if ( buffer->lines[i].line_len == 6 ){
				c_coefficients.push_back( buffer->lines[i].pop_double(2) );
				cP_coefficients.push_back( buffer->lines[i].pop_double(2) );
			}
		}else if ( buffer->lines[i].line_len == 1 ) jj++;
	}
	
	int kk = 1;
	int ii = 0;
	for( i=0;i<shell_n.size();i++){		
		if 		( i == shell_n.size()-1  ){
			shell_n_size.push_back(ii);
			if ( shell_n[i] == ++kk )  shell_n_size.push_back(1);
		}
		else if ( shell_n[i] == kk ) ii++; 		
		else{
			kk++;
			shell_n_size.push_back(ii);
			ii=1;
		}		
	}	
	
	buffer.reset(nullptr);
	
	int cont = 1;
	int rr   = 0;
	int ss   = 0;
	int pp   = 0;
	for( i=0;i<shell_n.size();i++){
		if ( shell_n[i] == cont ){
			if ( shell_t[i] == "S" ){
				Iaorbital orbS;
				for ( j=0;j<shell_n_size[shell_n[i]-1];j++){ orbS.add_primitive( exponents[ss++],c_coefficients[rr++] ); }
				molecule->atoms[atom_n_basis[i]].add_orbital(orbS);			
				cont++;
			}else if ( shell_t[i] == "L" ){
				Iaorbital orbS;
				Iaorbital orbpx, orbpy, orbpz;
				for ( j=0;j<shell_n_size[shell_n[i]-1];j++){
					orbS.add_primitive(exponents[ss],c_coefficients[rr++]);
					orbpz.add_primitive(exponents[ss++],cP_coefficients[pp++]);
				}
				orbpx = orbpy = orbpz;
				orbpx.powx = 1;
				orbpy.powy = 1;
				orbpz.powz = 1;
				orbpx.symmetry = "PX";
				orbpy.symmetry = "PY";
				orbpz.symmetry = "PZ";
				molecule->atoms[atom_n_basis[i]].add_orbital(orbS);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbpx);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbpy);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbpz);
				cont++;
			}else if ( shell_t[i] == "D" ){
				Iaorbital orbdxx, orbdyy, orbdzz, orbdxy, orbdxz, orbdyz;
				for (j=0;j<shell_n_size[shell_n[i]-1];j++){ orbdyz.add_primitive(exponents[ss++],c_coefficients[rr++]); }
				orbdxx = orbdyy =  orbdzz =  orbdxy =  orbdxz =  orbdyz;
				orbdxx.powx = 2;
				orbdyy.powy = 2;
				orbdzz.powz = 2;
				orbdxy.powy = orbdxy.powx = 1;
				orbdxz.powz = orbdxz.powx = 1;
				orbdyz.powz = orbdyz.powy = 1;
				orbdxx.symmetry = "XX";
				orbdyy.symmetry = "YY";
				orbdzz.symmetry = "ZZ";
				orbdxy.symmetry = "XY";
				orbdxz.symmetry = "XZ";
				orbdyz.symmetry = "YZ";
				molecule->atoms[atom_n_basis[i]].add_orbital(orbdxx);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbdyy);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbdzz);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbdxy);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbdxz);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbdyz);
				cont++;
			}else if ( shell_t[i] == "F" ){
				Iaorbital orbfxxx, orbfyyy, orbfzzz, orbfxxy, orbfyyx, orbfyyz;
				Iaorbital orbfxxz, orbfzzy, orbfzzx, orbfxyz;
				for (j=0;j<shell_n_size[shell_n[i]-1];j++){ orbfxxx.add_primitive(exponents[ss++],c_coefficients[rr++]); }
				orbfyyy=orbfzzz=orbfxxy=orbfyyx=orbfyyz=orbfxxz=orbfzzy=orbfzzx=orbfxyz=orbfxxx;
				orbfxxx.powx = 3;
				orbfyyy.powy = 3;
				orbfzzz.powz = 3;
				orbfxyz.powx = orbfxxz.powy = orbfxxz.powz = 1;
				orbfxxy.powx = 2;
				orbfxxy.powy = 1;
				orbfyyz.powy = 2;
				orbfyyz.powz = 1;
				orbfyyx.powy = 2;
				orbfyyx.powx = 1;
				orbfxxz.powx = 2;
				orbfxxz.powz = 1;
				orbfzzy.powy = 1;
				orbfzzy.powz = 2;
				orbfzzx.powx = 1;
				orbfzzx.powz = 2;
				orbfxxx.symmetry = "XXX";
				orbfyyy.symmetry = "YYY";
				orbfzzz.symmetry = "ZZZ";
				orbfxyz.symmetry = "XYZ";
				orbfxxy.symmetry = "XXY";
				orbfxxz.symmetry = "XXZ";
				orbfyyx.symmetry = "YYX";
				orbfyyz.symmetry = "YYX";
				orbfzzx.symmetry = "ZZX";
				orbfzzy.symmetry = "ZZY";
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfxxx);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfyyy);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfzzz);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfxxy);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfxxz);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfyyx);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfyyz);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfzzx);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfzzy);
				molecule->atoms[atom_n_basis[i]].add_orbital(orbfxyz);
				cont++;
			}
		}
	}

	vector<int>().swap(atom_n_basis);
	vector<int>().swap(shell_n);
	vector<int>().swap(shell_n_size);
	vector<string>().swap(shell_t);
	vector<double>().swap(exponents);
	vector<double>().swap(cP_coefficients);
	vector<double>().swap(c_coefficients);	
	
	int col_n = 0;
	int row_n = 0;
	int col_c = 0;
	int line_indicator = 0;	
	
	int aonum = molecule->get_ao_number();
	
	//-----------------------------------------------------
	// filling the occupied orbitals vector of the molecule
	if ( noe%2 == 0 ){
		molecule->occupied.resize(aonum);
		for(i=0;i<noe/2;i++) molecule->occupied[i] = 2;
	}else{
		molecule->occupied.resize(aonum);
		molecule->occupied_beta.resize(aonum);
		molecule->betad = true;
		for(i=0;i<(noe+1)/2;i++) molecule->occupied[i] = 1;
		for(i=0;i<(noe-1)/2;i++) molecule->occupied_beta[i] = 1;
			
	}
	//-----------------------------------------------------
	
	buffer.reset( new Ibuffer(file_name,fmo_in,fmo_fin));
	for( j=1;j<buffer->nLines;j++){		
		if ( buffer->lines[j].line_len > 0 && line_indicator == 0 ) {
			col_n = buffer->lines[j].pop_int(0);
			line_indicator++;
		}
		else if ( buffer->lines[j].line_len > 0 && line_indicator == 1){
			for( k=0;k<buffer->lines[j].line_len;k++){
				molecule->orb_energies.push_back(buffer->lines[j].pop_double(0) );
				molecule->MOnmb++;
			}
			line_indicator++;
		}
		else if ( buffer->lines[j].line_len > 0 && line_indicator == 2 )
			line_indicator = 3;
		else if ( buffer->lines[j].line_len >= 5 && line_indicator == 3 ){
			row_n = buffer->lines[j].pop_int(0);
			if ( molecule->coeff_MO.size() < molecule->MOnmb*aonum )
				molecule->coeff_MO.resize(molecule->MOnmb*aonum);
			for(l=0;l<buffer->lines[j].line_len-4;l++){ 
				molecule->coeff_MO[(col_n+l-1)*aonum+row_n-1] = buffer->lines[j].pop_double(3);
			}
			if ( row_n == aonum ) line_indicator = 0;
		}
	}
	
	if ( fmob_in > 0 ){
		buffer.reset( new Ibuffer(file_name,fmob_in,fmob_fin));
		for( j=4;j<buffer->nLines;j++){	
			if ( buffer->lines[j].line_len > 0 && line_indicator == 0 ) {
				col_n = buffer->lines[j].pop_int(0);
				line_indicator++;
			}
			else if ( buffer->lines[j].line_len > 0 && line_indicator == 1 ){
				for( k=0;k<buffer->lines[j].line_len;k++){
					molecule->orb_energies_beta.push_back( buffer->lines[j].pop_double(0) );
					molecule->MOnmb_beta++;
				}
				line_indicator++;
			}
			else if ( buffer->lines[j].line_len >  0 && line_indicator == 2 ) line_indicator = 3;
			else if ( buffer->lines[j].line_len >= 5 && line_indicator == 3 ){
				row_n = buffer->lines[j].pop_int(0);
				if ( molecule->coeff_MO_beta.size() < molecule->MOnmb*aonum )
					molecule->coeff_MO_beta.resize(molecule->MOnmb*aonum);
				for(l=0;l<buffer->lines[j].line_len-4;l++){
					molecule->coeff_MO_beta[(col_n+l-1)*aonum+row_n-1] = buffer->lines[j].pop_double(3);
				}
				if ( row_n == aonum ) line_indicator = 0;
			}
		}
	}
	
	buffer.reset(nullptr);
	
	molecule->get_homo_energy();
	molecule->get_lumo_energy();
	
	molecule->energy_tot  = molecule->energy_tot*27.2114;
	molecule->homo_energy = molecule->homo_energy*27.2114; // conversion to electronvolt
	molecule->lumo_energy = molecule->lumo_energy*27.2114;
	
	m_log->input_message("Total Energy ");
	m_log->input_message(double_to_string(molecule->energy_tot));
	m_log->input_message("HOMO energy: ");
	m_log->input_message(double_to_string(molecule->homo_energy));
	m_log->input_message("LUMO energy: ");
	m_log->input_message(double_to_string(molecule->lumo_energy));
	
	//molecule->print();
	molecule->norm_orbs();
	//molecule->print_basis();
	molecule->bohr_to_ang();
	
	parsed = true;
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_orca_out(const char* file_name){	
	if ( !check_file_ext(".out",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	
	m_log->input_message("Starting to parse out file from ORCA.");
	
	unsigned i,j   = 0;
	int in_coords  = 0;
	int fin_coords = 0;
	int orbs_in    = 0;
	int orbs_fin   = 0;
	int chg_in     = 0;
	int chg_fin    = 0;
	
	name_f = file_name;
	molecule->name = remove_extension(file_name);
	
	buffer.reset( new Ibuffer (file_name,true) );
	for (i=0;i<buffer->nLines;i++){
		if ( buffer->lines[i].IF_line("CARTESIAN",0,"(ANGSTROEM)",2,3) ){ in_coords = i; }
		else if ( buffer->lines[i].IF_line("CARTESIAN",0,"(A.U.)",2,3) ){ fin_coords = i; }
		else if ( buffer->lines[i].IF_line("Number",0,"Electrons",2,6) ){
			std::stringstream nofel(buffer->lines[i].words[5]);
			nofel >> molecule->num_of_electrons;
		}
		else if ( buffer->lines[i].IF_line("Total",0,"Energy",1,7) ){
			std::stringstream energ(buffer->lines[i].words[5]);
			energ >> molecule->energy_tot;
		}
		else if ( buffer->lines[i].IF_line("NO",0,"OCC",1,4) ) { orbs_in = i; }
		else if ( buffer->lines[i].IF_line("MULLIKEN",1,"ANALYSIS",3,5) ) { orbs_fin = i; }
		else if ( buffer->lines[i].IF_line("MULLIKEN",0,"CHARGES",2,3) ) { chg_in = i; }
		else if ( buffer->lines[i].IF_line("Sum",0,"charges:",2,5) ) { chg_fin = i; }
	}
	for(j = in_coords;j<fin_coords;j++){
		if ( buffer->lines[j].line_len == 4 ){
			double xx, yy, zz;
			xx = buffer->lines[j].get_double(1);
			yy = buffer->lines[j].get_double(2);
			zz = buffer->lines[j].get_double(3);
			string type_= buffer->lines[j].get_string(0);
			molecule->add_atom(xx,yy,zz,type_);
		}
	}
	
	for(i=(orbs_in +1);i<orbs_fin;i++){
		if ( buffer->lines[i].line_len == 4 ) 
			molecule->orb_energies.push_back( buffer->lines[i].get_double(3) ); 
	}
	
	int counter = 0;
	for(i=chg_in;i<chg_fin;i++){
		if ( buffer->lines[i].line_len == 4 )
			molecule->atoms[counter++].charge = buffer->lines[i].get_double(3);
	}

	buffer.reset();
	molecule->homo_energy = molecule->orb_energies[molecule->num_of_electrons/2 - 1];
	molecule->lumo_energy = molecule->orb_energies[molecule->num_of_electrons/2];
	parsed = true;
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_gaussian_fchk(const char* file_name){
	if ( !check_file_ext(".fchk",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	
	name_f         = file_name;	
	molecule->name = remove_extension(file_name);
	
	unsigned int i,j,k,l;
	int atomic_n_i = 0;
	int atomic_n_f = 0;
	int coord_n_i  = 0;
	int coord_n_f  = 0;
	int shell_t_i  = 0;
	int shell_t_f  = 0;
	int prim_b_i   = 0;
	int prim_b_f   = 0;
	int shell_m_i  = 0;
	int shell_m_f  = 0;
	int prim_e_i   = 0;
	int prim_e_f   = 0;
	int cont_c_i   = 0;
	int cont_c_f   = 0;
	int conspt_c_i = 0;
	int conspt_c_f = 0;
	int alpha_e_i  = 0;
	int alpha_e_f  = 0;
	int beta_e_i   = 0;
	int beta_e_f   = 0;
	int alpha_c_i  = 0;
	int alpha_c_f  = 0;
	int beta_c_i   = 0;
	int beta_c_f   = 0;
	int chgs_i     = 0;
	int chgs_f     = 0;
	int dens_i     = 0;
	int dens_f     = 0;
	vector<double> coords;
	vector<int> shell_t;
	vector<int> ngtos;
	vector<int> shell_map;
	vector<double> cont_c;
	vector<double> cont_c_p;
	vector<double> expos;
	
	int mult = 1;

	m_log->input_message("Starting to parse fchk file from gaussian.");
	
	unique_ptr<Ibuffer> b_file( new Ibuffer(file_name,true) );
	for(i=0;i<b_file->lines.size();i++){
		if ( b_file->lines[i].IF_line("Multiplicity",0,"I",1,3) ) { mult = b_file->lines[i].pop_int(2); }
		if ( b_file->lines[i].IF_line("Number",0,"electrons",2,5) ) { molecule->num_of_electrons = b_file->lines[i].pop_int(4); } 
		else if ( b_file->lines[i].IF_line("Number",0,"alpha",2,6) ) {
			int nalfa = b_file->lines[i].pop_int(5);
			if ( mult%2 == 0) for(j=0;j<nalfa;j++) { molecule->occupied.push_back(1); }
			else for(j=0;j<molecule->num_of_electrons;j++) { molecule->occupied.push_back(1); }
		}
		else if ( b_file->lines[i].IF_line("Number",0,"beta",2,6) ) {
			if ( mult%2 == 0 ) {
				int nbeta = b_file->lines[i].pop_int(5);
				for(j=0;j<nbeta;j++) { molecule->occupied_beta.push_back(1); }
				molecule->betad = true;
			}
		}
		else if ( b_file->lines[i].IF_line("Number",0,"beta",2,5) ) { molecule->num_of_electrons = b_file->lines[i].pop_int(4); } 
		else if	( b_file->lines[i].IF_line("Atomic",0,"numbers",1,5) ){ atomic_n_i = i;	}
		else if ( b_file->lines[i].IF_line("Nuclear",0,"charges",1,5) ){ atomic_n_f = i; }
		else if ( b_file->lines[i].IF_line("cartesian",1,"coordinates",2,6) ){ coord_n_i = i; }
		else if ( b_file->lines[i].IF_line("Force",0,"Field",1,4) ){ if ( coord_n_f == 0 ) coord_n_f =i; }
		else if ( b_file->lines[i].IF_line("Shell",0,"types",1,5) ){ shell_t_i = i; }
		else if ( b_file->lines[i].IF_line("primitives",2,"shell",4,8) ){ shell_t_f = prim_b_i = i; }
		else if ( b_file->lines[i].IF_line("Shell",0,"map",3,7) ){ shell_m_i = prim_b_f = i; }
		else if ( b_file->lines[i].IF_line("Primitive",0,"exponents",1,5) ){ shell_m_f = prim_e_i = i; }
		else if ( b_file->lines[i].IF_line("Contraction",0,"coefficients",1,5) ){ prim_e_f = cont_c_i = i; }
		else if ( b_file->lines[i].IF_line("Contraction",1,"coefficients",2,6) ){ conspt_c_i = cont_c_f = i; }
		else if ( b_file->lines[i].IF_line("Coordinates",0,"shell",3,7) ){ conspt_c_f = i; }		
		else if ( b_file->lines[i].IF_line("Total",0,"Energy",1,4) ){ molecule->energy_tot = b_file->lines[i].pop_double(3); }
		else if ( b_file->lines[i].IF_line("Alpha",0,"Energies",2,6) ){ alpha_e_i = i; }
		else if ( b_file->lines[i].IF_line("Beta",0,"Energies",2,6) ){ beta_e_i  = alpha_e_f = i; }
		else if ( b_file->lines[i].IF_line("Alpha",0,"coefficients",2,6) ){
			if ( beta_e_i > 0 ) beta_e_f  = alpha_c_i = i;
			else alpha_e_f = alpha_c_i = i;
		}
		else if ( b_file->lines[i].IF_line("Beta",0,"coefficients",2,6) ){  beta_c_i  = alpha_c_f = i; }
		else if ( b_file->lines[i].IF_line("Total",0,"Density",2,6) ){
			dens_i = i;
			if ( beta_c_i > 0 )	beta_c_f = i;
			else alpha_c_f = i;
		}
		else if ( b_file->lines[i].IF_line("Mulliken",0,"Charges",1,5) ){ dens_f = chgs_i = i;}
		else if ( b_file->lines[i].IF_line("Optimization",0,"MaxStp",1,4) ){ chgs_f = i;}
		else if ( b_file->lines[i].IF_line("ONIOM",0,"Charges",1,5) ){ if ( chgs_f == 0 ) chgs_f = i;}
	}
	k=0;
	int noe = molecule->num_of_electrons;
	
	for(i=0;i<b_file->lines.size();i++){
		if ( i>atomic_n_i && i<atomic_n_f){
			for(j=0;j<b_file->lines[i].line_len;j++ ){
				Iatom atom;
				atom.set_type( get_atomic_symbol( b_file->lines[i].pop_int(0) ) );
				molecule->add_atom(atom);
			}		
		}
		else if ( i>coord_n_i && i<coord_n_f ){
			for(j=0;j<b_file->lines[i].line_len;j++ ){
				coords.push_back(b_file->lines[i].pop_double(0));
			}
		}
		else if( i>shell_t_i && i< shell_t_f){
			for(j=0;j<b_file->lines[i].line_len;j++){
				shell_t.push_back(b_file->lines[i].pop_int(0));
			}
		}
		else if( i>prim_b_i && i<prim_b_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				ngtos.push_back(b_file->lines[i].pop_int(0));
			}
		}	
		else if( i>shell_m_i && i<shell_m_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				shell_map.push_back(b_file->lines[i].pop_int(0));
			}
		}
		else if( i>prim_e_i && i<prim_e_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				expos.push_back(b_file->lines[i].pop_double(0));
			}
		}
		else if( i>cont_c_i && i<cont_c_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				cont_c.push_back(b_file->lines[i].pop_double(0));
			}
		}
		else if( i>conspt_c_i && i<conspt_c_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				cont_c_p.push_back(b_file->lines[i].pop_double(0));
			}
		}
		else if( i>alpha_e_i && i<alpha_e_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				molecule->orb_energies.push_back(b_file->lines[i].pop_double(0));
				molecule->MOnmb++;
			}
		}
		else if( i>beta_e_i && i<beta_e_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				molecule->orb_energies_beta.push_back(b_file->lines[i].pop_double(0));
				molecule->MOnmb_beta++;
			}
		}
		else if( i>alpha_c_i && i<alpha_c_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				molecule->coeff_MO.push_back(b_file->lines[i].pop_double(0));
			}
		}
		else if( i>beta_c_i && i<beta_c_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				molecule->coeff_MO_beta.push_back(b_file->lines[i].pop_double(0));
			}
		}
		else if( i>dens_i && i<dens_f){
			for(j=0;j<b_file->lines[i].line_len;j++){
				molecule->m_dens.push_back( b_file->lines[i].pop_double(0) );
			}
		}
		else if( i>chgs_i && i<chgs_f ){
			for(j=0;j<b_file->lines[i].line_len;j++){
				molecule->atoms[j].charge = b_file->lines[i].pop_double(0);
			}
		}
	}

	j=0;
	for(i=0;i<molecule->atoms.size();i++){
		molecule->atoms[i].xcoord = coords[j++];
		molecule->atoms[i].ycoord = coords[j++];
		molecule->atoms[i].zcoord = coords[j++];
	}
	k = 0;
	int r = 0;
	for(i=0;i<shell_map.size();i++){
		if 		( shell_t[i] == 0 ){
			Iaorbital orbS;
			orbS.symmetry = "S";
			for (j=0;j<ngtos[i];j++){ orbS.add_primitive(expos[k++],cont_c[r++]); }
			molecule->atoms[shell_map[i]-1].add_orbital(orbS);
		}
		else if ( shell_t[i] == -1 ){
			
			vector<double> ex_p;
			vector<double> conc_p; 
			
			Iaorbital orbS;
			orbS.symmetry = "S";
			for (j=0;j<ngtos[i];j++){
				double tmp_ex = expos[k++];
				double tmp_cc = expos[r++];
				orbS.add_primitive(tmp_ex,tmp_cc);
				ex_p.push_back(tmp_ex);
				conc_p.push_back(cont_c_p[k-1]);
			}
			molecule->atoms[shell_map[i]-1].add_orbital(orbS);
			
			Iaorbital orbPX,orbPY,orbPZ;
			for (j=0;j<ngtos[i];j++){ orbPX.add_primitive(ex_p[j],conc_p[j]); }
			orbPY = orbPZ = orbPX;
			orbPX.symmetry = "PX";
			orbPX.powx     = 1;
			orbPY.symmetry = "PY";
			orbPY.powy     = 1;
			orbPZ.symmetry = "PZ";
			orbPZ.powz     = 1;
			molecule->atoms[shell_map[i]-1].add_orbital(orbPX);
			molecule->atoms[shell_map[i]-1].add_orbital(orbPY);
			molecule->atoms[shell_map[i]-1].add_orbital(orbPZ);
		}
		else if ( shell_t[i] == 1 ){
			Iaorbital orbPX,orbPY,orbPZ;
			for (j=0;j<ngtos[i];j++){ orbPX.add_primitive(expos[k++],cont_c[r++]); }
			orbPY = orbPZ = orbPX;
			orbPX.symmetry = "PX";
			orbPX.powx     = 1;
			orbPY.symmetry = "PY";
			orbPY.powy     = 1;
			orbPZ.symmetry = "PZ";
			orbPZ.powz     = 1;
			molecule->atoms[shell_map[i]-1].add_orbital(orbPX);
			molecule->atoms[shell_map[i]-1].add_orbital(orbPY);
			molecule->atoms[shell_map[i]-1].add_orbital(orbPZ);
		}
		else if ( shell_t[i] == 2 ){						
			Iaorbital orbDx2,orbDy2,orbDz2,orbDxy,orbDyz,orbDxz;
			for (j=0;j<ngtos[i];j++){ 
				orbDx2.add_primitive(expos[k++],cont_c[r++]);
			}
			orbDy2 = orbDz2 = orbDxy = orbDyz = orbDxz = orbDx2;
			orbDx2.symmetry = "XX";
			orbDx2.powx     = 2;
			orbDy2.symmetry = "YY";
			orbDy2.powy     = 2;
			orbDz2.symmetry = "ZZ";
			orbDz2.powz     = 2;
			orbDxy.symmetry = "XY";
			orbDxy.powx     = 1;
			orbDxy.powy     = 1;
			orbDyz.symmetry = "YZ";
			orbDyz.powy     = 1;
			orbDyz.powz     = 1;
			orbDxz.symmetry = "XZ";
			orbDxz.powx     = 1;
			orbDxz.powz     = 1;			
			molecule->atoms[shell_map[i]-1].add_orbital(orbDx2);
			molecule->atoms[shell_map[i]-1].add_orbital(orbDy2);
			molecule->atoms[shell_map[i]-1].add_orbital(orbDz2);
			molecule->atoms[shell_map[i]-1].add_orbital(orbDxy);
			molecule->atoms[shell_map[i]-1].add_orbital(orbDxz);
			molecule->atoms[shell_map[i]-1].add_orbital(orbDyz);
		}
		//else if ( shell_t[i] == 3 ){
			
		//}
	}
	
	int aonum = molecule->get_ao_number();
	if ( noe%2 == 0 ){
		molecule->occupied.resize(aonum);
		for(i=0;i<aonum;i++) molecule->occupied[i] = 0;
		for(i=0;i<noe/2;i++) molecule->occupied[i] = 2;
	}else{
		molecule->occupied.resize(aonum);
		molecule->occupied_beta.resize(aonum);
		molecule->betad = true;
		for(i=0;i<aonum;i++) molecule->occupied[i] = 0;
		for(i=0;i<aonum;i++) molecule->occupied_beta[i] = 0;
		for(i=0;i<(noe+1)/2;i++) molecule->occupied[i] = 1;
		for(i=0;i<(noe-1)/2;i++) molecule->occupied_beta[i] = 1;
			
	}
	
	molecule->norm_orbs();
	molecule->get_homo_energy();
	molecule->get_lumo_energy();
	
	m_log->input_message("HOMO energy: ");
	m_log->input_message(molecule->homo_energy);
	m_log->input_message("LUMO energy: ");
	m_log->input_message(molecule->lumo_energy);
	
	molecule->print_basis();
	molecule->bohr_to_ang();
	
	parsed = true;	
	return parsed;
}
/************************************************************************************/
bool QMparser::parse_gaussian_overlap(const char* file_name){
	if ( !check_file_ext(".log",file_name) ) {
		cout << "Warning! The file has wrong etension name!" << endl;
		m_log->input_message("Warning! The file has wrong etension name!");
		return false;
	}
	int over_in  = 0;
	int over_fin = 0;
	
	molecule->m_overlap.resize((molecule->MOnmb*(molecule->MOnmb+1))/2);
		
	buffer.reset(new Ibuffer(file_name,false));
	for(int i=0;i<buffer->lines.size();i++){
		if ( buffer->lines[i].IF_line("***",0,"Overlap",1,3) ) over_in  = i;
		if ( buffer->lines[i].IF_line("***",0,"Kinetic",1,4) ) over_fin = i;
	}
	
	int col_n = 0;
	int row_n = 0;
	int col_c = 0;
	
	for(int i=over_in+1;i<over_fin;i++){
		if ( buffer->lines[i].words.size() == 1 ) col_n = stoi(buffer->lines[i].words[0]);
		else if( buffer->lines[i].words.size() > 1 && buffer->lines[i].words[1].size() < 6) col_n = stoi(buffer->lines[i].words[0]);
		else{
			
			row_n = stoi(buffer->lines[i].words[0]) -1;
			col_c = col_n -1;
			for(int j=1;j<buffer->lines[i].line_len;j++){
				molecule->m_overlap[col_c + (row_n*(row_n+1))/2] = buffer->lines[i].pop_double_f(1);
				cout << row_n << " " << col_c << " " <<molecule->m_overlap[col_c + (row_n*(row_n+1))/2] << endl;
				col_c++;
			}
		}
	}
	
	return true;
}
/************************************************************************************/
Imolecule& QMparser::get_molecule() { return *molecule; }
/************************************************************************************/
QMparser::~QMparser(){}
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////END OF CLASS//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
