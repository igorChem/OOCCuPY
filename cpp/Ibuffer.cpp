#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "common.h"
#include "log_class.h"
#include "Iline.h"
#include "Ibuffer.h"

using std::move;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::unique_ptr;

/******************************************************************/
Ibuffer::Ibuffer():
	nLines(0)     ,
	name("buffer"){   
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name ,
				bool 			parse ):
	nLines(0)                          ,
	name(file_name)                    {
	
	if ( parse ){
		if ( IF_file(file_name) ){
			m_log->input_message("Storing all lines in Ibuffer object to be parsed.\n");
			ifstream buf(file_name);
			name  = file_name;
			string tmp_line;
			while( !buf.eof() ){
				getline(buf,tmp_line);
				lines.emplace_back( tmp_line );
				nLines++;
			}
			buf.close();
			parsed = true;
		}else{
			string message = "Not possible to open the file: ";
			message+= name;
			cout << message << endl;
			m_log->input_message(message);
			parsed = false;
		}
	}else{
		if ( IF_file(file_name) ){
			m_log->input_message("Checking if the file can be paser and how many lines it has.\n");
			ifstream buf(file_name);
			name  = file_name;
			string tmp_line;
			while( !buf.eof() ){
				getline(buf,tmp_line);
				nLines++;
			}
			buf.close();
		}else{
			string message = "Not possible to open the file: ";
			message+= name;
			cout << message << endl;
			m_log->input_message(message);
		}
		parsed = false;
	}
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name,
				int in                ,
				int fin)              :
				nLines(0)             ,
				name(file_name)       {
	
	int in_indx  = in;
	int fin_indx = fin;
	
	m_log->input_message("Searching for block in the file starting in the line: ");
	m_log->input_message(in_indx);
	m_log->input_message("And ending in the line: ");
	m_log->input_message(fin_indx);
	
	if ( IF_file(file_name) ){
		if ( in_indx !=0 ){
			m_log->input_message("Starting to store the block in a Ibuffer object to be parsed.");
			ifstream buf(file_name);
			nLines = 0;
			string tmp_line;
			int real_nlines = 0;
			while( !buf.eof() ){
				getline(buf,tmp_line);
				if ( nLines > in_indx && nLines < fin_indx  ){
					lines.emplace_back( move (tmp_line) );
					real_nlines++;
				}
				if ( nLines == fin_indx ) {	break; }
				nLines++;
			}
			nLines = real_nlines;
			buf.close();
			parsed = true;
		}else{
			m_log->input_message("Nothing to read!\n");
			nLines = 0;
			parsed = false;
		}
	}else{
		string message = "Not possible to open the file: ";
		message+= name;
		message+= "\n";
		cout << message << endl;	
		m_log->input_message(message);
		parsed = false;
	}
	
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name,
				 string wrdin         ,
				 string wrdfin       ):
	nLines(0)                         ,
	name(file_name)                   {
		
	int in_indx  = -1;
	int fin_indx = 0;
	
	m_log->input_message("Searching for block in the file starting with: ");
	m_log->input_message(wrdin);
	m_log->input_message("And ending with: ");
	m_log->input_message(wrdfin);
	
	if ( IF_file(file_name) ){
		std::ifstream buf(file_name);
		string tmp_line;
		while( !buf.eof() ){
			getline(buf,tmp_line);
			unique_ptr<Iline> Line ( new Iline(tmp_line) );
			if ( Line->IF_word( wrdin,0,wrdin.size() ) ) {
				in_indx = nLines;
			}
			if ( Line->IF_word( wrdfin,0,wrdfin.size() ) ) {
				fin_indx = nLines;
				break;
			}
			nLines++;
		}
		buf.close();
		if ( in_indx >= 0 ){
			m_log->input_message("Starting to store the block in a Ibuffer object to be parsed.");
			ifstream buf2(file_name);
			nLines = 0;
			string tmp_line2;
			int real_nlines = 0;
			while( !buf2.eof() ){
				getline(buf2,tmp_line2);
				if ( nLines > in_indx && nLines < fin_indx  ){
					lines.emplace_back( move (tmp_line2) );
					real_nlines++;
				}
				if ( nLines == fin_indx ) {	break; }
				nLines++;
			}
			nLines = real_nlines;
			buf2.close();
			parsed = true;
		}else{
			m_log->input_message("Nothing to read!\n");
			nLines = 0;
			parsed = false;
		}
	}else{
		string message = "Not possible to open the file: ";
		message+= name;
		message+= "\n";
		cout << message << endl;	
		m_log->input_message(message);
		parsed = false;
	}	
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name   , 
				vector<string>& wrds_in  ,
				vector<string>& wrds_fin ):
	nLines(0)                            ,
	name(file_name)                      {
	
	int in_indx  = -1;
	int fin_indx = 0;
	
	m_log->input_message("Searching for block in the file starting with: ");
	m_log->input_message(wrds_in[0]);
	m_log->input_message("And ending with: ");
	m_log->input_message(wrds_fin[0]);
	
	if ( IF_file(file_name) ){
		ifstream buf(file_name);
		string tmp_line;
		while( !buf.eof() ){
			getline(buf,tmp_line);
			unique_ptr<Iline> Line ( new Iline(tmp_line) );
			for(unsigned int i=0;i<wrds_in.size(); i++){
				if ( Line->IF_word( wrds_in[i],0,wrds_in[i].size() ) ) {
					if ( in_indx == 0 ) in_indx = nLines;
				}
			}
			for(unsigned int j=0;j<wrds_fin.size();j++){
				if ( Line->IF_word( wrds_fin[j],0,wrds_fin[j].size() ) ) {
					if ( fin_indx == 0 ) fin_indx = nLines;
					break;
				}	
			}
			nLines++;
		}
		if (in_indx >=0 ){
			buf.close();
			ifstream buf2(file_name);
			nLines = 0;
			string tmp_line2;
			int real_nlines = 0;
			while( !buf2.eof() ){
				getline(buf2,tmp_line2);
				if ( nLines > in_indx && nLines < fin_indx  ){
					lines.emplace_back( move(tmp_line2) );
					real_nlines++;
				}
				if ( nLines == fin_indx ) {	break; }
				nLines++;
			}
			nLines = real_nlines;
			buf2.close();
		}else{
			m_log->input_message("Nothing to read!\n");
			nLines = 0;
		}
	}else{
		string message = "Not possible to open the file: ";
		message += name;
		message +="\n";
		cout << message << endl;
		m_log->input_message(message);
	}
}
/******************************************************************/
Ibuffer::Ibuffer(const Ibuffer& rhs_ibuf):
	name(rhs_ibuf.name)                  ,
	lines(rhs_ibuf.lines)                ,
	nLines(rhs_ibuf.nLines)              {
}
/******************************************************************/
Ibuffer& Ibuffer::operator=(const Ibuffer& rhs_ibuf){	
	if( this!=&rhs_ibuf ){
		name   = rhs_ibuf.name;
		lines  = rhs_ibuf.lines;
		nLines = rhs_ibuf.nLines;
	}
	return *this;
}  
/******************************************************************/
Ibuffer::Ibuffer(Ibuffer&& rhs_ibuf) noexcept:
	name( move(rhs_ibuf.name) )              ,
	lines( move(rhs_ibuf.lines) )            ,
	nLines( rhs_ibuf.nLines )                {
}
/******************************************************************/
Ibuffer& Ibuffer::operator=(Ibuffer&& rhs_ibuf) noexcept {
	if( this!=&rhs_ibuf ){
		name   = move(rhs_ibuf.name);
		lines  = move(rhs_ibuf.lines);
		nLines = rhs_ibuf.nLines;
	}
	return *this;
}
/******************************************************************/
Ibuffer& Ibuffer::get_block(int in, int fin){
	if ( parsed ){
		vector<Iline> temp;
		for(int i=in;i<fin;i++) { temp.emplace_back( move(lines[i] ) ); }
		lines = move(temp);
		return *this;
	}else{
		if ( IF_file(name) ){
			if ( in !=0 ){
				ifstream buf(name);
				nLines = 0;
				string tmp_line;
				int real_nlines = 0;
				while( !buf.eof() ){
					getline(buf,tmp_line);
					if ( nLines > in && nLines < fin  ){
						lines.emplace_back( move (tmp_line) );
						real_nlines++;
					}
					if ( nLines == fin ) {	break; }
					nLines++;
				}
				nLines = real_nlines;
				buf.close();
				parsed = true;
			}else{
				nLines = 0;
				parsed = false;
			}
		}else{
			string message = "Not possible to open the file: ";
			message+= name;
			message+= "\n";
			cout << message << endl;	
			parsed = false;
		}
		return *this;
	}
}
/******************************************************************/
Iline* Ibuffer::pop_front_line(){
	Iline* Line = new Iline(lines[0]);
	lines.erase(lines.begin());
	return Line;
}
/******************************************************************/
void Ibuffer::shrink_lines(){ lines.shrink_to_fit(); }
/******************************************************************/
void Ibuffer::print(){
	cout << name   << endl;
	cout << nLines << endl;
}
/******************************************************************/
Ibuffer::~Ibuffer(){}
/******************************************************************/