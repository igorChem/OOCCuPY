//global.h
//Declarations for functions and source file with constant values. 

/*********************************************************************/
/* This source code file is part of OOCCuPy++ software project created 
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

#ifndef GLOBAL
#define GLOBAL

#include <iostream>
#include <string>
#include <vector> 
#include <memory>
#include <experimental/filesystem>

bool IF_file(const char* name); 
bool IF_file(std::experimental::filesystem::path& name); 
bool check_file_ext(std::string ext,const char* file_name);
std::string remove_extension(const char* file_name);
std::string change_extension(const char* file_name,std::string new_ext);
void rename_file(const char* file_name,std::string new_file_name);
std::string get_file_name(const char* path_name);
std::string str_array(std::string& line, int in, int fin);

#endif
