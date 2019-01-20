//Iatom.h

#ifndef IATOM
#define IATOM

//including c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <cstring>

//including PRIMoRDiA headers
#include "common.h"

class Iaorbital; // foward declaration

//=================================================================================================
/**
 * This class is meant to hold values for atomic properties, also manipulate and modify, from quantum
 * chemistry programs outputs
 * 
 * @class Iatom
 * @author Igor Barden Grillo - barden.igor@gmail.com
 * @date 20/03/18
 * @file Iatom.h
 * @brief A Iatom Class to instatiate atom representation object to Quantum Mechanical properties calculations. 
 */
class Iatom {
	public:
		//--------------------------------------------------------------------
		/**
		 * @brief Double with x cartesian coordinate for the atom. 
		 */
		double xcoord; 
		
		//--------------------------------------------------------------------
		/**
		 * @brief Double with y cartesian coordinate for the atom.
		 */
		double ycoord;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Double with z cartesian coordinate for the atom.
		 */
		double zcoord;
		
		//--------------------------------------------------------------------
		/**
		 * @brief String with element symbol for the atom.
		 */
		std::string element;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Float with value of partial charge of the atom.
		 */
		float charge;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Float with value of atomic mass of the atom.
		 */
		float atomic_mass;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Integer value of atomic number.
		 */
		unsigned int atomicN;
		
		//--------------------------------------------------------------------
		/**
		 * @brief Integer with the number of atomic orbitals for the atom.
		 */
		unsigned int norb;
		
		//--------------------------------------------------------------------
		/**
		 * @brief STL vector container with Iaorbital objects of atomic orbitals for the atom.
		 */	 
		std::vector<Iaorbital> orbitals;
		
		//--------------------------------------------------------------------
		/**
		 * Instantiate a empty Iatom object to be filled later and initializes the member data to 
		 * default values
		 * @brief Default constructor for Iatom object.
		 */
		Iatom();
		
		//--------------------------------------------------------------------
		/**
		 * Four arg constructor to Iatom object from the basic information required
		 * for label atoms.
		 * @brief Instatiate a Iatom object with minimal information of 
		 * coordinates and element symbol.
		 * @param double value with x coordinate of the atom.
		 * @param double value with y coordinate of the atom.
		 * @param double value with z coordinate of the atom.
		 * @param string with symbol with of the atom. 
		 */
		Iatom(double x,double y,double z, std::string typ); 
	
		//--------------------------------------------------------------------
		/** 
		 * Instantiate a Iatom object as a copy of another Iatom object.
		 * @brief Class copy constructor. 
		 * @param Iatom constant reference to a object to be copied.
		 */
		Iatom(const Iatom& rhs_atom); 
		
		//--------------------------------------------------------------------
		/** 
		 * Overwrite the current iatom object information with another object.
		 * @brief Assigment operator overloading.
		 * @param Iatom constant reference to a object to be copied.
		 * @return Iatom reference to this object.
		 */
		Iatom& operator=(const Iatom& rhs_atom);
		
		//--------------------------------------------------------------------
		/**
		 * @brief Move constructor.
		 * @param Iatom rvalue reference of object to be moved to this.
		 */
		Iatom(Iatom&& rhs_atom) noexcept;
 
		//--------------------------------------------------------------------
		/**
		 * @brief  Move assigment operator overloading.
		 * @param  Iatom rvalue reference of object to be moved to this.
		 * @return Iatom reference to this object.
		 */
		Iatom& operator=(Iatom&& rhs_atom) noexcept;
		
		//--------------------------------------------------------------------
		 /**
		 * @brief Print basic atom information to the console.
		 * @return None. 
		 */
		void print();
		
		//--------------------------------------------------------------------
		/**
		 * Memeber function to set change the symbol of the Iatom element.
		 * @brief Set element symbol of the Iatom object. 
		 * @param String reference to atom symbol to be setted.
		 * @return None.
		 */			 		
		void set_type(const std::string typ);		
		
		//--------------------------------------------------------------------
		/**
		 * Memeber function to set the three dimensional coordinates of the Iatom objects.
		 * @brief Set cartesian coordinates to Iatom object. 
		 * @param double  to x coordinate value.
		 * @param double  to y coordinate value.
		 * @param double  to z coordinate value.
		 * @return None.
		 */
		void set_coord(double x, double y, double z);	
	
		//--------------------------------------------------------------------	
		/**
		 * Member function to create and store an Iaorbital object in the orbitals memeber function of 
		 * the current Iatom object.
		 * @brief Create an Iaorbital Object and store it in the current object.
		 * @param Ineger with the level of orbital shell type.
		 * @param String reference with symmetry of the atomic orbital to be created.
		 * @param Double to zeta coefficient of atomic orbital. 
		 * @return None.
		 */		 
		void add_orbital(int level,const std::string sym, double coef);
		
		//--------------------------------------------------------------------
		/**
		 * Member function to store a Iaorbital object in the orbitals member variable 
		 * of the current object.. 
		 * @brief Store Iaorbital object passed.
		 * @param Const reference to a Iaorbital oject. 
		 * @return None. 
		 */
		void add_orbital(Iaorbital orb);
		
		//--------------------------------------------------------------------
		/**
		 * Destructor to the Iatom class.
		 */
		~Iatom();
};

#endif