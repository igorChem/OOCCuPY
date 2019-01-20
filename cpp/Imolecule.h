// header file for Imolecule.cpp
// Imolecule.H

/* header file with the class declaration for chemical and geometrical information representation */

#ifndef IMOLECULE
#define IMOLECULE

//including c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

//including PRIMoRDiA headers
#include "common.h"

class Iaorbital;
class Iatom; // foward declaration
//========================================================================
/**
 * @class Imolecule
 * @author Igor Barden Grillo
 * @date 20/03/18
 * @file Imolecule.h 
 * @brief Imolecule class to hold abstract representation for molecules and its quantum chemical informations
 * get from computational methods.
 * 
 * Main class of the program to hold quantum chemical molecular information, store atomic representations
 * objects.
 */
class Imolecule {
	public:
		//--------------------------------------------------------------------
		/**
		 * @brief String with molecule's name.
		 */ 
		std::string name;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Integer with number of atoms in the molecule.
		 */
		unsigned int num_of_atoms;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Integer value that hold the total number of electrons in the molecule from
		 * the sum of atomic numbers of Iatom objects.
		 */
		unsigned int num_of_electrons;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Float value with molecular mass calculated from the atomic masses
		 * from Iatom objects.
		 */
		float molar_mass;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Float value with net molecular charge calculatedc by the partial charges 
		 * from Iatom objects.
		 */
		float mol_charge;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Double value of total energy of the molecule calculated with some quantum mechanical method.
		 */
		double energy_tot;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Double value of energy of the Highest energy Occupied Molecular Orbital.
		 */
		double homo_energy;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Double value of energy of the Loawest energy Unoccupied Molecular Orbital.
		 */
		double lumo_energy;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Double value to hold the total dipole moment.
		 */
		double total_dipmoment;
		
		//----------------------------------------------------------------------
		/**
		 * @brief Double value to hold the molecule's heat of formation 
		 */
		double heat_of_formation;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Integer number indicating the number of molecular orbitals in the object. 
		 */
		unsigned int MOnmb;	
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Integer number indicating the number of molecular orbitals in the object. 
		 */
		unsigned int MOnmb_beta;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Integer number indicating the number of molecular orbitals in the object. 
		 */
		int homoN;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Integer number indicating the number of molecular orbitals in the object. 
		 */
		int lumoN;
		
		//-------------------------------------------------------------------
		/**
		 * @brief Bool type with information if the MO are normalized
		 */
		bool normalized;
		
		//-------------------------------------------------------------------
		/**
		 * @brief Bool type with information if the MO are normalized
		 */
		bool bohr;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Bool indicating if beta density matrix was stored.
		 */
		bool betad;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Double  with coordinate of inferior vertice of cube that 
		 * enclose the molecule to calculate quantum proprieties in a grid.
		 */
		double ver_inf[3];
		
		//--------------------------------------------------------------------
		/**
		 * @brief Double  with coordinate of superior vertice of cube that
		 * enclose the molecule to calculate quantum proprieties in a grid.
		 */
		double ver_sup[3];
		
		//--------------------------------------------------------------------
		/**
		 * @brief  Double array with the three dimensional components of dipole moment calculated with quantum
		 * mechanical calculations.
		 */
		double dipole_moment[3];
				
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector to doubles values of orbital energies.
		 */
		std::vector<double> orb_energies;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief Integer number indicating the number of molecular orbitals in the object. 
		 */
		std::vector<double> orb_energies_beta;
		
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector to doubles values of the molecular orbitals coefficients array.
		 */
		std::vector<double> coeff_MO;
				
		//--------------------------------------------------------------------
		/**
		 *  @brief Integer number indicating the number of molecular orbitals in the object. 
		 */
		std::vector<double> coeff_MO_beta;
		
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector to doubles with the values of density matrix.
		 */
		std::vector<double> m_dens;
		
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector to doubles with values of beta density matrix.
		 */
		std::vector<double> beta_dens;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Double pointer to array of overlap matrix
		 */ 
		std::vector <double> m_overlap;		
		
		//--------------------------------------------------------------------
		/**
		 *  @brief  STL vector of integers with the occupation of the molecular orbitals
		 */
		std::vector <int> occupied;
		
		//--------------------------------------------------------------------
		/**
		 *  @brief STL vector of integers with the occupation of the molecular orbitals for
		 * the beta set.
		 */
		std::vector <int> occupied_beta;
		
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector of Iatom objects with the instantiated objects that 
		 * represent the atoms of the molecule.
		 */
		std::vector<Iatom> atoms;
		
		//--------------------------------------------------------------------
		/**
		 * Constructor of the Imolecule class to instantiate a empty object to be filled later.
		 * @brief Default constructor of Imolecule.
		 */
		Imolecule(); 
		
		//--------------------------------------------------------------------
		/**
		 * Instantiates a class object from a copy of another one.
		 * @brief Imolecule copy constructor. 
		 * @param Constant reference to Imolecule object to be copied.
		 */
		Imolecule(const Imolecule& rhs_molecule);
		
		//--------------------------------------------------------------------
		/**
		 * Overwrite current object from member variables from another Imolecule object in the
		 * RHS assigment operators.
		 * @brief Assigmet operator overloading. 
		 * @param Constant reference to an Imolecule object to assign to the current object.
		 * @return Reference of the current object.
		 */
		Imolecule& operator=(const Imolecule& rhs_molecule);
		
		//--------------------------------------------------------------------
		/**
		 * @brief Move constructor
		 * @param Imolecule rvalue reference to be moved.
		 */
		Imolecule(Imolecule&& rhs_molecule) noexcept;
		
		/**
		 * @brief Move assigment operator overloading.
		 * @param Imolecule rvalue to be moved.
		 * @return Imolecule reference to this.
		 */
		Imolecule& operator=(Imolecule&& rhs_molecule) noexcept;
		
		//--------------------------------------------------------------------
		/**
		 * Instantiates an Iatom object with the basic information provided and stores in the atoms STL 
		 * vector of the current object followed of a molecule information up-date.  
		 * @brief Member function with four-args of the class to add an Iatom object to the current object.
		 * @param Double with x coordinate to the atom object.
		 * @param Double with y coordinate to the atom object.
		 * @param Double with z coordinate to the atom object.
		 * @param String reference with element symbol of the atom.
		 * @return None.
		 */
		void add_atom(double x,double y,double z,std::string typ);
		
		//--------------------------------------------------------------------
		/**
		 * Member function with one-arg of the class to add an Iatom object to the current object. 
		 * @brief Copies an Iatom object to the atoms STL vector of the current object 
		 * @param Reference to Iatom to be stored in the Imolecule object.
		 * @return None. 
		 */
		void add_atom(Iatom atom);
		
		//--------------------------------------------------------------------
		/**
		 * Contant Member function to print all cartesian coordinates of the molecule object.
		 * @brief Print to the console all molecule cartesian coordinates. 
		 * @return None.
		 */
		void print_coordinates();
		
		//--------------------------------------------------------------------
		/**  
		 * Creates a file in xyz coordinate styles to be readed from 
		 * molecular visualizarion programs.
		 * @brief Memebr function to write xyz file style from molecular info.
		 * @return none.
		 */
		void write_xyz();
		
		//--------------------------------------------------------------------
		/**
		 * From the molecular coordinates values this function calculates the
		 * inferior ans superior vertices coordinates to build a cube that enclose the molecule
		 * to further quantum mechanical calculations.		 * 
		 * @brief Calculates the max and minimum values of vertices cooridantes of the molecule. 
		 * @return none.
		 */
		void mol_vert_up();
		
		//--------------------------------------------------------------------
		/**
		 * @brief Return a STL vector of doubles with the coefficients of a specific molecular orbital.
		 * @param Integer with the number of the MO,
		 * @param Bool indicating if the molecular orbital is from beta set.
		 * @return STL vector of doubles.
		 */
		std::vector<double> extract_MO(int MO,bool beta);
		
		//--------------------------------------------------------------------
		/**
		 * @brief Class method to compute normalized  molecular orbitals coefficients 
		 * @return None.
		 */
		void output_DOS();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Print molecular information to the console of the current object.
		 * @return None.
		 */
		void print();
		  
		//-------------------------------------------------------------------
		/**
		 * @brief Convert the atoms coordinates from agnstrom to bohr.
		 * @return None.
		 */
		void ang_to_bohr();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Convert the atoms coordinates from bohr to angstrom.
		 * @return None.
		 */		
		void bohr_to_ang();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Class method force the calculation of normalization factors
		 * of the atomic orbitals.
		 * @return None.
		 */
		void norm_orbs();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Count and return the total number of atomic orbitals. 
		 * @return Integer with the number of atomic orbitals in the molecule.
		 */
		int get_ao_number();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Search in the orb_energies and orb_energies_beta vectors the index of 
		 * the highest occupied energy molecular orbital.
		 * @return Integer with the index number of the HOMO orbital.
		 */
		int get_homo_number();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Get the double value with the energy of HOMO orbital.
		 * @return Double with the energy(eV) of the HOMO.
		 */
		double get_homo_energy();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Search in the orb_energies and orb_energies_beta vectors the index of 
		 * the lowest energy unnoccupied molecular orbital.
		 * @return Integer with the index number of the LUMO orbital.
		 */
		int get_lumo_number();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Get the double value with the energy of LUMO orbital.
		 * @return Double with the energy(eV) of the LUMO.
 		 */
		double get_lumo_energy(); 
		
		//-------------------------------------------------------------------
		/**
		 * @brief Print to the console the atomic orbitals information and its 
		 * gaussian primitives if has any.
		 * @return None.
 		 */
		void print_basis();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Update the information dependent of the atoms in the molecule.
		 * @return None.
 		 */
		void update();
		
		//-------------------------------------------------------------------
		/**
		 * @brief Erases the elements of all STL vectors in the curret Imolecule
		 * object and releases memory.
		 * @return None.
 		 */
		void clear();
		
		//--------------------------------------------------------------------
		/**
		 * Deallocate memory allocated in the constructors
		 *  @brief Destructor of the class.   
		 */
		~Imolecule();
};

#endif