//Iaorbital.h

#ifndef IAORBITAL
#define IAORBITAL

//inclding c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
//including PRIMoRDiA headers
#include "common.h"

//-------------------------------------------------------------------------
/**
 * This class is meant to hold objects with gaussian primitive atomic orbitals function
 * information, and manipulate and change this information, like calculating the normalization constants.
 * @class Iprimitive
 * @author Igor Barden Grillo
 * @date 11/09/18
 * @file Iaorbital.h
 * @brief Gaussian type orbital primitive for the contracted GTO atomic orbital class.
 */
class Iprimitive {
	public:
		//----------------------------------------------------------------
		/**
		 * @brief Double for the normalization factor 
		 * as a product of calculated normalization factor with the contraction coefficient.
		 */
		double n_fact;
		
		//----------------------------------------------------------------
		/**
		 * @brief Double for the exponent value readed from QM output for the current gaussian primitive.
		 */
		double exponent;
		
		//----------------------------------------------------------------
		/**
		 * @brief Double holding the contraction coefficient readed from QM output file for the current 
		 * gaussian primitive.
		 */
		double c_coef;
		
		//----------------------------------------------------------------
		/**
		 * Instantiates an gaussian primitve object without any parameter and 
		 * initializes the member variables to default values.
		 * @brief Default constructor.
		 */
		Iprimitive();
		
		//----------------------------------------------------------------
		/**
		 * Instantiates the object given the values of primitive gaussian exponent
		 * and its contraction coefficient. Initializes the other member variables
		 * @brief Constructor with the exponent and contraction coefficient parameters.
		 * @param Double with the exponent value.
		 * @param Double with contraction coefficient.
		 */
		Iprimitive(double exp,double contrac_coeff);
		
		//----------------------------------------------------------------
		/**
		 * Instantiates the object from a copy of another.
		 * @brief Copy constructor.
		 * @param Iprimitive constant reference of the object to be copied.
		 */
		Iprimitive(const Iprimitive& prim_rhs);
		
		//----------------------------------------------------------------
		/**
		 * Instantiate the object from a moved information from another 
		 * object.
		 * @brief Move constructor.
		 * @param Iprimitive constant rvalue reference of the object to be moved.
		 */
		Iprimitive(Iprimitive&& prim_rhs) noexcept;
		
		//----------------------------------------------------------------
		/**
		 * Overwrites the object in LHS (the current object) of the assigment operator 
		 * with the member varibale information of the object in the RHS.
		 * @brief Assigment operator overloading for Iprimitive class.
		 * @param Iprimitive constant reference of the object to be copied.
		 * @return Iprimitive reference object.
		 */
		Iprimitive& operator=(const Iprimitive& prim_rhs);
		
		//----------------------------------------------------------------
		/**
		 * Overwrites the object in LHS (the current object) of the assigment operator 
		 * with moved member varibale information of the object in the RHS.
		 * @brief Move assigment operat0or overloading. 
		 * @param Iprimitive constant rvalue reference of the object to be moved.
		 * @return Iprimitive reference object.
		 */
		Iprimitive& operator=(Iprimitive&& prim_rhs) noexcept;
		
		//----------------------------------------------------------------
		/**
		 * @brief Class method to print to the console object Iaorbital information.
		 * @return None.
		 */
		void print();
		
		//----------------------------------------------------------------
		/**
		 * @brief Destructor.
		 */
		~Iprimitive();
	
};
//========================================================================
/** 
 * This class is meant to be hold atom object abstraction to hold, manipulate and modify
 * information extracted from quantum chemical output programs.
 * @class Iaorbital
 * @author Igor Barden Grillo.
 * @date 20/03/18
 * @file Iatom.h
 * @brief  A Iaorbital Class declaration to instatiate atomic orbital object to Quantum Mechanical 
 * properties calculations.  
 *   
 */ 
class Iaorbital {
	public:
		/**
		 * @brief Integer unsigned with the level of the atomic shell.
		 */
		unsigned int shell;
		
		//----------------------------------------------------------------
		/**
		 * @brief Bool with info if the orbital is gaussian type. 
		 * If is not the oarbital is slater.
		 */
		bool gto;
		
		//----------------------------------------------------------------
		/**
		* @brief String with atomic orbital symmetry symbol.
		*/
		std::string symmetry; 
		
		//----------------------------------------------------------------
		/**
		 * @brief Double with normalization constant of the current atomic orbital object.
		 */
		double n_factor;
		
		//----------------------------------------------------------------
		/**
		 * @brief Double with zeta coefficient information for atomic orbital.
		 */
		double alpha; 
		
		//----------------------------------------------------------------
		/**
		 * @brief Integer with exponential value for symmetry accounting in cartersian atomic orbital for x dimension.
		 */
		unsigned int powx;
		
		//----------------------------------------------------------------
		/**
		 * @brief Integer with exponential value for symmetry accounting in cartersian atomic orbital for y dimension.
		 */
		unsigned int powy;
		
		//----------------------------------------------------------------
		/**
		 * @brief Integer with exponential value for symmetry accounting in cartersian atomic orbital for z dimension.
		 */
		unsigned int powz;
		
		//----------------------------------------------------------------
		/**
		 * @brief STL vector of Iprimitive objects of the gaussian type orbitals of the atomic orbital.
		 */
		std::vector<Iprimitive> gtos;
  
		//----------------------------------------------------------------
		/**
		 * Instantiate a empty atomic orbital object to be filled later and 
		 * initializes the member varibale to default values.
		 * @brief Default constructor for Iaorbital class. 
		 */		
		Iaorbital();
		
		//----------------------------------------------------------------
		/**
		 * Constructor with three arguments with basic information to represent an atomic orbital.
		 * @brief Instantiate the object with information provided by Quantum Mechanical program output.
		 * @param Integer with shell level for the atomic orbital object. 
		 * @param string reference with symmetry symbol for the atomic orbital object.
		 * @param double  with zeta value for the atomic orbital object.
		 */
		Iaorbital(unsigned int level, std::string sym, double coef); 
		
		//----------------------------------------------------------------
		/**  
		 * Instantiates a object from a copy of another Iatom object.
		 * @brief Copy constructor for the class.  
		 * @param Iaorbital const reference with to be copied to this object. 
		 */
		Iaorbital(const Iaorbital& rhs_orb);
		
		//----------------------------------------------------------------
		/**
		 * Overloads the Assigment operator to the class.
		 * @brief Assigment operator overloading
		 * @param Iaorbital const reference to be assigned to this object. 
		 * @return a reference of Iaorbital object.
		 */
		Iaorbital& operator=(const Iaorbital& rhs_orb);
		
		//----------------------------------------------------------------
		/**
		 * Instantiates an object from the moved ownership member data from another class
		 * object.
		 * @brief move constructor of the class.
		 * @param Iaorbital rvalue reference to be moved.
		 */
		Iaorbital(Iaorbital&& rhs_orb) noexcept;
		
		//----------------------------------------------------------------
		/**
		 * Overwrites the member data with stolen ownership of memeber data from a object
		 * to be moved.
		 * @brief Move assigment operator overloading.
		 * @param Iaorbital rvalue reference to be moved.
		 * @return Iaorbital reference of this.
		 */
		Iaorbital& operator=(Iaorbital&& rhs_orb) noexcept;
		
		//----------------------------------------------------------------
		/**
		 * Create a gaussian type orbital primitive object from the exponent and contraction 
		 * coefficient values and stores as a element of STL gtos member variable.
		 * @brief Member function to create a gto primitive and add to the gtos vector.
		 * @param Double with the exponent value.
		 * @param Double with the contraction coefficient value.
		 * @return None.
		 */
		void add_primitive(double expo, double contrac);
		
		//-----------------------------------------------------------------
		/**
		 * Stores a gaussian type orbital primitive object in the as a element STL
		 * gtos member variable.
		 * @brief Member function to create a gto primitive and add to the gtos vector.
		 * @param Iprimitive oject.
		 * @return None.
		 */
		void add_primitive(Iprimitive gorb);
		
		//----------------------------------------------------------------
		/**
		 * Member function that calculates the normalization constant of the current Iaorbital object,
		 * using the other member varibles : shell, symmetry and exponent coefficient. 
		 * The calculation is implemented for cartesian atomic orbitals of Slater type.
		 * @brief Class method to calculate the normalization factor of the current aomtic
		 * orbital object.
		 * @return None.
		 */
		bool normalize();
		
		//----------------------------------------------------------------
		/**
		 * @brief Class method to print to the console object Iaorbital information.
		 * @return None.
		 */
		void print();
		
		//----------------------------------------------------------------
		/**
		 * @brief Destructor for the Iaorbital class.  
		 */
		~Iaorbital();
};

#endif