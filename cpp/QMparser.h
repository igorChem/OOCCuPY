//header file for QM output file parser
// currently done for mopac and GAMESS
// QMparser.h

#ifndef  QMPARSER
#define QMPARSER

//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <memory>

#include "common.h"

class Iline; //fowards declarations 
class Ibuffer;
class Imolecule;

/**
 * @class QMparser
 * @author igor Barden Grillo
 * @date 24/03/18
 * @file QMparser.h
 * @brief Class to read and stores files from output of quantum chemical calculations programs
 * to generate objects with info to calculate reactivity descriptors.
 * 
 * This class instatiates a Ibuffer object to load all the lines of a given file from a quantum
 * chemistry output calculation and parse. For now the member functions to parse file of gamess and mopac
 * are implemented, but the class make available the needed tools to parse formated file and store the chemical
 * relevant information. 
 */
class QMparser{
	public:
	/**
	 * @brief String with file name.
	 */
	std::string name_f;
	//---------------------------------------------------------------------------------------------------
	/**
	 * @brief std::string specifying the program type.
	 */
	std::string program;
	 
	//---------------------------------------------------------------------------------------------------
	/**
	 * @brief Ibuffer object class with the buffered fille.
	 */
	std::unique_ptr<Ibuffer> buffer;

	//---------------------------------------------------------------------------------------------------
	/**
	 * @brief Imolecule pointer to the object to be filled with molecular information.
	 */
	std::unique_ptr<Imolecule> molecule;
	 
	//---------------------------------------------------------------------------------------------------
	/**
	 * @brief Indicating if the file was parsed. 
	 */
	bool parsed;
	
	//-------------------------------------------------------------------------------------------------
	/**
	 * @brief Default constructor.
	 */
	QMparser();

	//--------------------------------------------------------------------------------------------------
	/**
	 * @brief Copy constructor. 
	 * @param QMparser reference object to be copied to the current object. 
	 */
	QMparser(const QMparser& rhs_QMp);
	 
	//-------------------------------------------------------------------------------------------------
	/**
	 * @brief Assigment operator overloading.
	 * @param QMparser object to be copied.
	 * @return QMparser reference to this object.
	 */
	QMparser& operator=(const QMparser& rhs_QMp);
	
	//--------------------------------------------------------------------------------------------------
	/**
	 * @brief Move constructor. 
	 * @param QMparser reference object to be moved to the current object. 
	 */
	QMparser(QMparser&& rhs_QMp) noexcept;
	 
	//-------------------------------------------------------------------------------------------------
	/**
	 * @brief Move assigment operator overloading.
	 * @param QMparser object to be moved.
	 * @return QMparser reference to this object.
	 */
	QMparser& operator=(QMparser&& rhs_QMp) noexcept;
	 
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse mopac aux file for molecular and quantum information.
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_mopac_aux(const char* file_name);
	
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse mopac out file for molecular and quantum information.
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_out_mopac(const char* file_name);
	
	//-----------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse MOZYME mopac aux file
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_mozyme(const char* file_name);
	
	//-----------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse MOZYME mopac mgf file
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_mopac_mgf(const char* file_name);
	
	//-----------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse orca output file.
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_orca_out(const char* file_name);
	
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Method to parse log gamess program.
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_log_gamess(const char* file_name);
	
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse gaussian fchk file.
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_gaussian_fchk(const char* file_name);
	
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to parse gaussian log file to get overlap matrix.
	 * @param Char pointer constant indicating the name of the file to be parsed.
	 * @return Bool indicating if the file was parsed.
	 */
	bool parse_gaussian_overlap(const char* file_name);
	
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Class method to return a Imolecule object.
	 * @return Imolecule reference of the molecule member variable.
	 */
	Imolecule& get_molecule();
	 
	//------------------------------------------------------------------------------------------------
	/**
	 * @brief Destructor of the class.
	 */
	~QMparser();
};
 
#endif
