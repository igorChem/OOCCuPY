//commom.h
//Declarations for functions and source file with constant values. 

#ifndef COMMON
#define COMMON
//--------------------------------------------------------------------------------
//Including header from the c++
#include <iostream>
#include <string>
#include <vector> 
#include <map>
#include <memory>
//--------------------------------------------------------------------------------
//Including header from PRIMoRDiA library. 
#include "log_class.h"
#include "Itimer.h"

//===========================================================
// GLOBAL VARIABLES for internal usage: DECLARION
//===========================================================
//---------------------------------------------------------
/**
 * @brief Global variable to indicate if gnuplot scripts for radial distribution 
 * of scalar field data for gnuplot.
 */
extern bool gnuplot;

//---------------------------------------------------------
/**
 * @brief Global variable
 */ 
extern bool extra_RD;

//---------------------------------------------------------
/**
 * @brief Global variable
 */
extern bool pymol_script;

//---------------------------------------------------------
/**
 * @brief Global variable
 */
extern bool timer;

//---------------------------------------------------------
/**
 * @brief Global variable
 */
extern unsigned int NP; 

//---------------------------------------------------------
/**
 * @brief Global variable
 */
extern Itimer chronometer;

//---------------------------------------------------------
/**
 * @brief Global variable
 */
extern std::unique_ptr<Ilog> m_log;

//===========================================================
// GLOBAL VARIABLES for main function usage DECLARION
//===========================================================
//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_file_in;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_ED;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_MO;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_verbose;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern bool M_logfile;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern int  M_gridsize;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern int  M_orbNum;

//---------------------------------------------------------
/**
 * @brief Global variable for main function.
 */
extern std::string M_program;


//------------------------------------------------------------------------------------------
/**
 * @fn get_atom_mass
 * @brief Function to get the atomic mass from array initialized in source file.
 * @param String indicating the element symbol.
 * @return Float with atomic mass value.
 */
float get_atom_mass(std::string sym);

//-----------------------------------------------------------------------------------------
/**
 * @fn get_atomic_number
 * @brief Function to get the atomic number from array initialized in source file.
 * @param String indicating the element symbol.
 * @return Integer with the atomic number.
 */
int get_atomic_number(std::string sym);

//-----------------------------------------------------------------------------------------
/**
 * @fn get_atomic_symbol
 * @brief Function to get the atomic symbol from array initialized in source file.
 * @param Integer indicating the atomic nuber.
 * @return String with the atomic number.
 */
std::string get_atomic_symbol(int i);

//---------------------------------------------------------------------------------------
/**
 * @fn factorial
 * @brief Calculate the factorial of an integer value.
 * @param Integer with the value to get the factorial from.
 * @return Integer with of the factorial of the integer passed. 
 */
int factorial(int n);

//---------------------------------------------------------------------------------------
/**
 * @fn double factorial
 * @brief Calculate the semi-factorial/double factorial of an integer value.
 * @param Integer with the value to get the factorial from.
 * @return Integer with of the factorial of the integer passed. 
 */
int doublefactorial(int n);

//---------------------------------------------------------------------------------------
/**
 * @fn int_to_string
 * @brief Convert integer value to an std string.
 * @param Integer to be converted to string.
 * @return String from integer passed.
 */
std::string int_to_string(int val);

//---------------------------------------------------------------------------------------
/**
 * @fn int_to_string
 * @brief Convert integer value to an std string.
 * @param Integer to be converted to string.
 * @return String from integer passed.
 */
std::string double_to_string(double val);

//---------------------------------------------------------------------------------------
/**
 * @fn string_to_int
 * @brief Convert std string to integer type.
 * @param String to be converted.
 * @return Integer from string passed.
 */
int string_to_int(std::string str);

//---------------------------------------------------------------------------------------
/**
 * @fn IF_file
 * @brief Test if file can be open.
 * @param Const char with the name of the file.
 * @return Bool if the file can be open.
 */
bool IF_file(const char* name); 

//---------------------------------------------------------------------------------------
/**
 * @fn check_file_ext
 * @brief Test if file the file extension is the one that is required.
 * @param Const char with the name of the file.
 * @return Bool if the file can be open.
 */
bool check_file_ext(std::string extension,const char* file_name);

//---------------------------------------------------------------------------------------
/**
 * @fn remove_extension
 * @brief return a string with the filename without the extension part 
 * @param constant char* containing the file name.
 * @return string with the filename without the extension.
 */
std::string remove_extension(const char* file_name);

//---------------------------------------------------------------------------------------
/**
 * @fn change_extension
 * @brief Return a string with a file name with a new extension
 * @param char pointer constant with the file name
 * @param string with the new extension
 * @return string of the file name without the extension
 */
std::string change_extension(const char* file_name,std::string new_ext);

//---------------------------------------------------------------------------------------
/**
 * @brief Convert double in string format from file with "D" to "E"
 * @param string with the double to be converted
 * @return double with the converted values
 */
double D_E_conv(std::string sc_not);


#endif