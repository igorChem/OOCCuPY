// source file for the Imolecule class to represent and deal with moleculae information
// imolecule.cpp 

#include <iostream> 
#include <vector>
#include <string>  
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "log_class.h"
#include "common.h"
#include "Iaorbital.h"
#include "Iatom.h"
#include "Imolecule.h"

using std::cout; 
using std::endl;
using std::move;
using std::string;
using std::vector;

//=======================================================================================
/***************************************************************************************/
Imolecule::Imolecule()    :
	name("nonamed")       ,
	num_of_atoms(0)       ,
	num_of_electrons(0)   ,
	molar_mass(0.0)       , 
	energy_tot(0.0)       ,
	homo_energy(0.0)      ,
	lumo_energy(0.0)      ,
	total_dipmoment(0.0)  ,
	heat_of_formation(0.0),
	MOnmb(0)              ,
	MOnmb_beta(0)         ,
	homoN(0.0)            ,
	lumoN(0.0)            ,
	normalized(false)     ,
	bohr(false)     ,
	betad(false)          {
	
	for (unsigned int i=0; i<3;i++){
		ver_inf[i]       = 0.0;
		ver_sup[i]       = 0.0;
		dipole_moment[i] = 0.0;
	}
} 
/***************************************************************************************/
Imolecule::Imolecule(const Imolecule& rhs_molecule)  :
	name(rhs_molecule.name)                          ,
	num_of_atoms(rhs_molecule.num_of_atoms)          ,
	num_of_electrons(rhs_molecule.num_of_electrons)  ,
	molar_mass(rhs_molecule.molar_mass)              ,
	energy_tot(rhs_molecule.energy_tot)              ,
	homo_energy(rhs_molecule.homo_energy)            ,
	lumo_energy(rhs_molecule.lumo_energy)            , 
	heat_of_formation(rhs_molecule.heat_of_formation),
	MOnmb(rhs_molecule.MOnmb)                        ,
	MOnmb_beta(rhs_molecule.MOnmb_beta)              ,
	homoN(rhs_molecule.homoN)                        ,
	lumoN(rhs_molecule.lumoN)                        ,
	normalized(rhs_molecule.normalized)              ,
	bohr(rhs_molecule.bohr)                          ,
	betad(rhs_molecule.betad)                        ,
	orb_energies(rhs_molecule.orb_energies)          ,
	orb_energies_beta(rhs_molecule.orb_energies_beta),
	coeff_MO(rhs_molecule.coeff_MO)                  ,
	coeff_MO_beta(rhs_molecule.coeff_MO_beta)        ,
	m_dens(rhs_molecule.m_dens)                      ,
	beta_dens(rhs_molecule.beta_dens)                ,
	m_overlap(rhs_molecule.m_overlap)                ,
	occupied(rhs_molecule.occupied)                  ,
	occupied_beta(rhs_molecule.occupied_beta)        ,
	atoms(rhs_molecule.atoms)                        {
	
	for(int i = 0 ; i < 3; i++) {
		dipole_moment[i] = rhs_molecule.dipole_moment[i];
		ver_inf[i]       = rhs_molecule.ver_inf[i];
		ver_sup[i]       = rhs_molecule.ver_sup[i];
	}
}
/***************************************************************************************/
Imolecule& Imolecule::operator=(const Imolecule& rhs_molecule){
	if(this != &rhs_molecule) {	
		name              = rhs_molecule.name;
		num_of_atoms      = rhs_molecule.num_of_atoms;
		num_of_electrons  = rhs_molecule.num_of_electrons;
		molar_mass        = rhs_molecule.molar_mass;
		energy_tot        = rhs_molecule.energy_tot;
		homo_energy       = rhs_molecule.homo_energy;
		lumo_energy       = rhs_molecule.lumo_energy;
		heat_of_formation = rhs_molecule.heat_of_formation;
		MOnmb             = rhs_molecule.MOnmb;
		MOnmb_beta        = rhs_molecule.MOnmb_beta;
		homoN             = rhs_molecule.homoN;
		lumoN             = rhs_molecule.lumoN;
		normalized        = rhs_molecule.normalized;
		bohr              = rhs_molecule.bohr;
		betad             = rhs_molecule.betad;
		orb_energies      = rhs_molecule.orb_energies;
		orb_energies_beta = rhs_molecule.orb_energies_beta;
		coeff_MO          = rhs_molecule.coeff_MO;
		coeff_MO_beta     = rhs_molecule.coeff_MO_beta;
		m_dens            = rhs_molecule.m_dens;
		beta_dens         = rhs_molecule.beta_dens;
		m_overlap         = rhs_molecule.m_overlap;
		occupied          = rhs_molecule.occupied;
		occupied_beta     = rhs_molecule.occupied_beta;
		atoms             = rhs_molecule.atoms;	
		
		for(int i = 0 ; i < 3; i++) {
			dipole_moment[i] = rhs_molecule.dipole_moment[i];
			ver_inf[i]       = rhs_molecule.ver_inf[i];
			ver_sup[i]       = rhs_molecule.ver_sup[i];
		}
	}
	return *this;
}
/***************************************************************************************/
Imolecule::Imolecule(Imolecule&& rhs_molecule) noexcept      :
	name(move (rhs_molecule.name) )                          ,
	num_of_atoms(rhs_molecule.num_of_atoms)                  ,
	num_of_electrons(rhs_molecule.num_of_electrons)          ,
	molar_mass(rhs_molecule.molar_mass)                      ,
	energy_tot(rhs_molecule.energy_tot)                      ,
	homo_energy(rhs_molecule.homo_energy)                    ,
	lumo_energy(rhs_molecule.lumo_energy)                    , 
	heat_of_formation(rhs_molecule.heat_of_formation)        ,
	MOnmb(rhs_molecule.MOnmb)                                ,
	MOnmb_beta(rhs_molecule.MOnmb_beta)                      ,
	homoN(rhs_molecule.homoN)                                , 
	lumoN(rhs_molecule.lumoN)                                ,
	normalized(rhs_molecule.normalized)                      ,
	bohr(rhs_molecule.bohr)                                  ,
	betad(rhs_molecule.betad)                                ,
	orb_energies( move(rhs_molecule.orb_energies) )          ,
	orb_energies_beta( move(rhs_molecule.orb_energies_beta) ),
	coeff_MO( move(rhs_molecule.coeff_MO) )                  ,
	coeff_MO_beta( move(rhs_molecule.coeff_MO_beta) )        ,
	m_dens( move(rhs_molecule.m_dens) )                      ,
	beta_dens( move(rhs_molecule.beta_dens) )                ,
	m_overlap( move(rhs_molecule.m_overlap) )                ,
	occupied( move(rhs_molecule.occupied))                   ,
	occupied_beta( move(rhs_molecule.occupied_beta) )        ,
	atoms( move(rhs_molecule.atoms) )                        {
		
	for(int i = 0 ; i < 3; i++) {
		dipole_moment[i] = rhs_molecule.dipole_moment[i];
		ver_inf[i]       = rhs_molecule.ver_inf[i];
		ver_sup[i]       = rhs_molecule.ver_sup[i];
	}
}
/***************************************************************************************/
Imolecule& Imolecule::operator=(Imolecule&& rhs_molecule) noexcept {
	if(this != &rhs_molecule) {	
		name              = move(rhs_molecule.name);
		num_of_atoms      = rhs_molecule.num_of_atoms;
		num_of_electrons  = rhs_molecule.num_of_electrons;
		molar_mass        = rhs_molecule.molar_mass;
		energy_tot        = rhs_molecule.energy_tot;
		homo_energy       = rhs_molecule.homo_energy;
		lumo_energy       = rhs_molecule.lumo_energy;
		heat_of_formation = rhs_molecule.heat_of_formation;
		MOnmb             = rhs_molecule.MOnmb;
		MOnmb_beta        = rhs_molecule.MOnmb_beta;
		homoN             = rhs_molecule.homoN;
		lumoN             = rhs_molecule.lumoN;
		normalized        = rhs_molecule.normalized;
		bohr              = rhs_molecule.bohr;
		betad             = rhs_molecule.betad;
		orb_energies      = move(rhs_molecule.orb_energies);
		orb_energies_beta = move(rhs_molecule.orb_energies_beta);
		coeff_MO          = move(rhs_molecule.coeff_MO);
		coeff_MO_beta     = move(rhs_molecule.coeff_MO_beta);
		m_dens            = move(rhs_molecule.m_dens);
		beta_dens         = move(rhs_molecule.beta_dens);
		m_overlap         = move(rhs_molecule.m_overlap);
		occupied          = move(rhs_molecule.occupied);
		occupied_beta     = move(rhs_molecule.occupied_beta);
		atoms             = move(rhs_molecule.atoms);
		
		for(int i = 0 ; i < 3; i++) {
			dipole_moment[i] = rhs_molecule.dipole_moment[i];
			dipole_moment[i] = rhs_molecule.dipole_moment[i];
			ver_inf[i]       = rhs_molecule.ver_inf[i];
			ver_sup[i]       = rhs_molecule.ver_sup[i];
		}
	}
	return *this;
}
/***************************************************************************************/
void Imolecule::add_atom(double x         ,
						 double y         , 
						 double z         , 
						 string typ)      {
							 
	atoms.emplace_back( x,y,z,move(typ) );
	num_of_atoms = atoms.size();
	this->update();
}
/***************************************************************************************/
void Imolecule::add_atom(Iatom atom) {
	atoms.emplace_back( move(atom) ) ; 
	num_of_atoms = atoms.size();

	mol_charge       = 0;
	molar_mass       = 0;
	num_of_electrons = 0;
	this->update();
}
/***************************************************************************************/
void Imolecule::print_coordinates() {
	for(int i=0;i<num_of_atoms;i++){
			 cout << "x # " << i << ": " << atoms[i].xcoord << " "
				  << "y # " << i << ": " << atoms[i].ycoord << " "
				  << "z # " << i << ": " << atoms[i].zcoord << endl;
	}
}
/***************************************************************************************/
void Imolecule::write_xyz() {
	string xyz_name;
	xyz_name = name + ".xyz";
	std::ofstream xyz_file(xyz_name.c_str());
	xyz_file << num_of_atoms << "\n";
	xyz_file << "xyz file generated by Imolecule class by Barden/2018" << "\n";
	m_log->input_message("Writing xyz!");
	for(int i=0;i<num_of_atoms;i++){ 
		xyz_file.precision(7);
		xyz_file << std::fixed;
		xyz_file << std::setw(3) << std::left << atoms[i].element
		<< "  " 
	    << std::setw(11) << std::right << atoms[i].xcoord 
	    << "  "
	    << std::setw(11) << std::right << atoms[i].ycoord
	    << "  "
	    << std::setw(11) << std::right << atoms[i].zcoord
	    << "\n";
	}
}
/***************************************************************************************/
void Imolecule::mol_vert_up(){
	ver_sup[0] = atoms[0].xcoord;
	ver_sup[1] = atoms[0].ycoord;
	ver_sup[2] = atoms[0].zcoord;
	
	ver_inf[0] = atoms[0].xcoord;
	ver_inf[1] = atoms[0].ycoord;
	ver_inf[2] = atoms[0].zcoord;
	
	m_log->input_message("Definind the grid vertices.");
	
	for (unsigned int i=0;i<atoms.size();i++){
		if ( atoms[i].xcoord > ver_sup[0] ) ver_sup[0] = atoms[i].xcoord;
		if ( atoms[i].ycoord > ver_sup[1] ) ver_sup[1] = atoms[i].ycoord;
		if ( atoms[i].zcoord > ver_sup[2] ) ver_sup[2] = atoms[i].zcoord;
		if ( atoms[i].xcoord < ver_inf[0] ) ver_inf[0] = atoms[i].xcoord;
		if ( atoms[i].ycoord < ver_inf[1] )	ver_inf[1] = atoms[i].ycoord;
		if ( atoms[i].zcoord < ver_inf[2] )	ver_inf[2] = atoms[i].zcoord;
	}
}
/***************************************************************************************/
vector<double> Imolecule::extract_MO(int MO,bool beta){
	int nMO = 0;
	if ( !beta ) nMO = MOnmb; 
	else         nMO = MOnmb_beta;
	vector<double> res_mo(nMO);	 
	for (int i=0;i<nMO;i++){
		if ( !beta ) res_mo[i] = coeff_MO[MO*MOnmb+i];
		else         res_mo[i] = coeff_MO_beta[MO*MOnmb_beta+i];
	}
	return res_mo;
}
/***************************************************************************************/
void Imolecule::output_DOS(){
	string Name = name +".DOS";
	string Name2 = name +"_DOS.gnus";
	std::ofstream dos_file( Name.c_str() );
	std::ofstream dos_file_gnu( (name +"_DOS.gnus").c_str() );
	if (!betad){ for(int i=0;i<MOnmb;i++) dos_file << orb_energies[i] << "\t " << i << endl; }
	
	m_log->input_message("Outputing Density of States information to file.");
	
	dos_file_gnu << "set grid \n"
				 << "set ylabel 'Density of States' \n"
				 << "set xlabel 'Orbital energies' \n"
				 << "set terminal pngcairo enhanced font \"arial,10\" fontscale 1.0 size 1024, 900 \n"
				 << "set output '" << Name << ".png' \n"
				 << "plot '" << Name << "'" << " using 1:2 smooth kdensity  with filledcurves above y lt 10 title ''";
	dos_file.close();
	dos_file_gnu.close();
}
/***************************************************************************************/
void Imolecule::print(){
	std::cout << "Molecule's name: " << name << std::endl;
	this->print_coordinates();
	std::cout << "The coeff are normalized :"   << normalized  << std::endl;
	std::cout << "LUMO energy  : "              << lumo_energy << std::endl;
	std::cout << "HOMO energy  : "              << homo_energy << std::endl;
	std::cout << "Total energy : "              << energy_tot  << std::endl;
	std::cout << "Printing the first ten values of the vector in this object" << std::endl;
	for (unsigned int o=0;o<num_of_atoms;o++) atoms[o].print();
	for (unsigned int i=0;i<MOnmb;i++) {
		std::cout << orb_energies[i] << " " 
				  << coeff_MO[i]  	 << " " 
				  << occupied[i]     << std::endl;
	}
}
/***************************************************************************************/
void Imolecule::ang_to_bohr(){
	const double angtobohr = 1.0/0.52917726;
	for(unsigned int i=0;i<num_of_atoms;i++){
		atoms[i].xcoord = atoms[i].xcoord*angtobohr;
		atoms[i].ycoord = atoms[i].ycoord*angtobohr;
		atoms[i].zcoord = atoms[i].zcoord*angtobohr;
	}
	bohr = true;
	m_log->input_message("Cartesian molecular coordinates converted to bohr.");
}
/***************************************************************************************/
void Imolecule::bohr_to_ang(){
	const double bohrtoang = 0.52917726;
	for(unsigned int i=0;i<num_of_atoms;i++){
		atoms[i].xcoord = atoms[i].xcoord*bohrtoang;
		atoms[i].ycoord = atoms[i].ycoord*bohrtoang;
		atoms[i].zcoord = atoms[i].zcoord*bohrtoang;
	}
	bohr = false;
	m_log->input_message("Cartesian molecular coordinates converted to angstron.");
}
/***************************************************************************************/
void Imolecule::norm_orbs(){
	for (unsigned int i=0;i<atoms.size();i++){
		for (unsigned int j=0;j<atoms[i].orbitals.size();j++){
			atoms[i].orbitals[j].normalize();
		}
	}
}
/***************************************************************************************/
int Imolecule::get_ao_number(){
	int aon = 0;
	for(unsigned int i=0;i<atoms.size();i++){
		for (unsigned int j=0;j<atoms[i].orbitals.size();j++) aon++;
	}
	return aon;
}
/***************************************************************************************/
int Imolecule::get_homo_number(){
	int homo_nalfa = homoN;
	int homo_nbeta = homoN;
	
	for (unsigned int i=0;i<occupied.size();i++){ 
		if ( occupied[i]>=1 )
			homo_nalfa = i;
	}
	for (unsigned int j=0;j<occupied_beta.size();j++) { 
		if ( occupied_beta[j]>=1 )
			homo_nbeta = j;
	}
	
	double homo_alfa  = orb_energies[homo_nalfa];
	double homo_beta  = -1000.0;
	
	if ( homo_nbeta != 0 ) homo_beta  = orb_energies_beta[homo_nbeta]; 
	
	if ( homo_alfa >= homo_beta ) {
		homoN = homo_alfa;
		return homo_nalfa;
	}
	else if ( homo_beta >  homo_alfa ){
		homoN = -homo_beta;
		return homo_nbeta;
	}
	return 0;
}
/***************************************************************************************/
double Imolecule::get_homo_energy(){
	homoN = this->get_homo_number();
	if ( homoN > 0)      return homo_energy = orb_energies[homoN];
	else if ( homoN < 0 ) return homo_energy = orb_energies_beta[-1*homoN]; 
}
/***************************************************************************************/
int Imolecule::get_lumo_number(){
	int lumo_nalfa = lumoN;
	int lumo_nbeta = lumoN;
	
	for (int i=0;i<occupied.size();i++){ 
		if ( occupied[i] < 1) {
			lumo_nalfa = i;
			break;
		}
	}
	for (int j=0;j<occupied_beta.size();j++) {
		if ( occupied_beta[j] < 1 ) {
			lumo_nbeta = j;
			break;
		}
	}
	
	double lumo_alfa = orb_energies[lumo_nalfa];
	double lumo_beta = 10000.0;
	if ( lumo_nbeta > 0) { lumo_beta = orb_energies_beta[lumo_nbeta]; }
	
	if ( lumo_alfa <= lumo_beta ) {
		lumoN = lumo_alfa;
		return lumo_nalfa;
	}
	else if ( lumo_beta <  lumo_alfa ){
		lumoN = -lumo_beta;
		return lumo_nbeta;
	}
}
/***************************************************************************************/
void Imolecule::print_basis(){
	for(int i=0;i<atoms.size();i++){
		for(int j=0;j<atoms[i].orbitals.size();j++)	atoms[i].orbitals[j].print();
	}
}
/***************************************************************************************/
void Imolecule::update(){
	for(int i=0;i<num_of_atoms;i++) {
		mol_charge       += atoms[i].charge;
		molar_mass       += atoms[i].atomic_mass;
		num_of_electrons += atoms[i].atomicN;
	}
}
/***************************************************************************************/
void Imolecule::clear(){
	vector<double>().swap(orb_energies);      
	vector<double>().swap(orb_energies_beta);
	vector<double>().swap(coeff_MO);
	vector<double>().swap(coeff_MO_beta);
	vector<double>().swap(m_dens);
	vector<double>().swap(beta_dens);
	vector<double>().swap(m_overlap);
	vector<int>().swap(occupied);
	vector<int>().swap(occupied_beta);
	vector<Iatom>().swap(atoms);
}
/***************************************************************************************/
double Imolecule::get_lumo_energy(){
	lumoN = this->get_lumo_number();
	if ( lumoN > 0)      return lumo_energy = orb_energies[lumoN];
	else if ( lumoN < 0 ) return lumo_energy = orb_energies_beta[-1*lumoN]; 
}
/***************************************************************************************/
Imolecule::~Imolecule(){}
/***************************************************************************************/
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////end of Imolecule.cpp/////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////