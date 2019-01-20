//commmon.cpp
/*
 * Source file for common utilities for the primordia libs program, 
 */
//--------------------------------------------------------------------------------
//Including header from the c++
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
//--------------------------------------------------------------------------------
//Including header from PRIMoRDiA library. 
#include "common.h"
//--------------------------------------------------------------------------------
//Including standard statements.
using std::string;
/*********************************************************************************/
bool gnuplot      = false;
bool extra_RD     = false;
bool pymol_script = false;
bool timer        = false;
unsigned int NP   = omp_get_max_threads();
/*********************************************************************************/
Itimer chronometer;
std::unique_ptr<Ilog> m_log ( new Ilog() );
/*********************************************************************************/
bool M_file_in    = false;
bool M_ED         = false;
bool M_MO         = false;
bool M_verbose    = false;
bool M_logfile    = false;
int  M_gridsize   = 0;
int  M_orbNum     = 0;
string M_program  = "noprogram";
/************************************************************************************************************/
string atomType[] = {
					 "H",		       																    "He",
					 "Li","Be",											  			"B","C","N","O","F","Ne",
					 "Na","Mg",										  			"Al","Si","P","S","Cl", "Ar",
					 "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Cu","Zn",			"Ga","Ge","As","Se","Br","Kr",
					 "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
					 "Cs","Ba",
							  "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
							 "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Ti","Pb","Bi","Po","At","Rn",
					 "Fr","Ra",
							  "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bh","Bk","Cf","Es","Fm","Md","No","Lr"
					 };
/************************************************************************************************************/
double atomMass[] = {
					1.00794,																				   		 				 			4.0026,
					6.941,9.012182,																	10.811,	12.0107,14.0067,15.9994,18.9984,	20.1797,
					22.9897,24.305,																	26.981,28.0855,30.973762,32.065,35.453,39.948,39.948,
					40.078,44.9559,47.867,50.9415,51.9961,54.938,55.845,58.6934,58.9332,			63.546,	65.39,69.723,72.64,74.9216,78.96,79.904,83.8,
					85.467,87.62l,88.9059,91.224,92.9064,95.94,98,101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,121.76,126.9045,127.6,131.293,
					132.9055,137.327,
									138.9055,140.116,140.9077,144.24,145,150.36,151.964,157.25,158.9253,162.5,164.9303,167.259,168.9342,173.04,174.967,
									178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,204.3833,207.2,208.9804,209,210,222,
					223,226,
							227,231.0359,232.0381,237,238.0289,243,244,247,247,251,252,257,258,259,261,262
					};
/************************************************************************************************************/
int factorial(int n){  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
/***********************************************************************************/
int doublefactorial(int n){  return (n == 1 || n == 0) ? 1 : factorial(n - 2) * n;}
/***********************************************************************************/
float get_atom_mass(string sym){
	int indx = 0;
	for (unsigned int i=0;i<103;i++)	if ( sym == atomType[i] ) indx = i;
	return atomMass[indx];
}
/*********************************************************************************/
int get_atomic_number(string sym){
	int atomic_num;
	for (unsigned int i=0;i<103;i++) if ( sym == atomType[i] ) atomic_num = i+1;
	return atomic_num;
}
/*********************************************************************************/
string get_atomic_symbol(int i){ return atomType[i-1]; }
/*********************************************************************************/
string int_to_string(int val){
	string result;
	std::ostringstream convert;
	convert << val;
	result = convert.str();
	return result;
}
/*********************************************************************************/
string double_to_string(double val){
	string result;
	std::ostringstream convert;
	convert << val;
	result = convert.str();
	return result;
}
/*********************************************************************************/
int string_to_int(string str){
	int result;
	std::stringstream convert(str);
	convert >> result;
	return result;
}
/*********************************************************************************/
bool IF_file(const char* name){ return ( access( name, F_OK ) != -1 );}
/*********************************************************************************/
bool check_file_ext(string extension,const char* file_name){
	string fname = file_name;
	int len_ext  = extension.size();
	int len_file = fname.size();
	fname = fname.substr(len_file-len_ext,len_file);
	if ( extension == fname ) return true;
	else return false;	
}
/*********************************************************************************/
string remove_extension(const char* file_name){
	string file_name_string  = file_name;
	int point = 0;
	char dot = '.';
	for(unsigned int i = 0;i<file_name_string.size();i++){
		char character = file_name_string[i];
		if ( character == dot){	point = i; }
	}
	int pos = file_name_string.size()-point;
	return file_name_string.substr(0,file_name_string.size()-pos);
}
/*********************************************************************************/
string change_extension(const char* file_name,string new_ext){
	string name_wth_ext = remove_extension(file_name);
	name_wth_ext += new_ext;
	return name_wth_ext;
}
/********************************************************************************/
double D_E_conv(string sc_not){
	std::replace(sc_not.begin(),sc_not.end(),'D','E');
	double result = std::stod(sc_not);
	return result;
}
/********************************************************************************/
//END OF FILE
//================================================================================