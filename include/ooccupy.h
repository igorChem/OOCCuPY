//ooccupy.h

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

#ifndef OOCCUPY
#DEFINE OOCCUPY

#include <string>
#include <vector>

class ooccupy{
	public:
		int m_argc;
		std::vector<std::string> m_argv;
		ooccupy();
		~ooccupy();
		ooccupy(int argc, char* argv[]);
		ooccupy(const ooccupy& rhs) = delete;
		ooccupy& operator=(const ooccupy& rhs) = delete;		
		void run();
		void print_options();
		void help();
};


#endif