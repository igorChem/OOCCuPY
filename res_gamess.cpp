#include <iostream>
#include <string>
#include <vector>
#include <fstream>


using std::move;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;


int main(int argc,char** argv){
	
	char tmp_line[1024];

	string input_lines;
	string data_group;
	string vec_group;

	string oname = argv[1];
	oname  = oname.substr(0,oname.size()-4);
	oname += "_vec.inp";
	std::ofstream final_input(oname.c_str());
	ifstream buf(argv[1]);
	bool dat_gr = false;
	bool inp_gr = true;

	int nLines = 0;
	int ini,fin = 0;
	while( !buf.eof() ){
		buf.getline(tmp_line,1024);
		string word(tmp_line,0,6);
		if ( word == " $DATA" || word == " $data"   ){
			dat_gr = true;
			inp_gr = false;
		}
		if ( inp_gr ){
			final_input << tmp_line << endl;
		}else if ( dat_gr ){
			data_group += tmp_line;
			data_group += "\n";
		}
	}
	buf.close();

 	final_input << " $guess guess=moread $end" << endl;
	ifstream buf2(argv[2]);
	while( !buf2.eof() ){
		buf2.getline(tmp_line,1024);
			string word(tmp_line,0,5);
			if ( word == " $VEC" ){
				ini = nLines;
			}
			if ( ini > 0 ){
				if ( word == " $END" ){
					fin = nLines;	
				}
			}
		nLines++;
	}
	buf2.close();

	nLines = 0;
	ifstream buf3(argv[2]);
	while( !buf3.eof() ){
		buf3.getline(tmp_line,1024);
		if ( nLines >= ini && nLines <= fin ){
			final_input << tmp_line << endl;
		}			
		nLines++;
	}
	string out_name = argv[1];
	final_input << data_group << endl;
	final_input.close();

	return 0;
}