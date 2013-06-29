/* File name: Fasta.h
 * Author:	Weilong GUO
 * Created Date:		2010-07-24
 * Last Modification:	2010-07-24
 */

#ifndef _FASTA_H
#define _FASTA_H

#include <string>
#include <iostream>

using namespace std;

struct seqline{
	string name;
	string seq;
};

/** Read in a file which is in FASTA format, and store in
 *  a 'vector<seqline>' structure fasta.
 */
int ReadFasta( string filename, vector<seqline> & fasta )
{
	fasta.clear();
	ifstream fa(filename.c_str());
	if(!fa) {
		cerr << "Cannot open fasta file " << filename << "\n";
		return -1;
	}
	seqline sl;
	sl.name = "";
	sl.seq = "";
	while(!fa.eof()){
		char buffer[10000 + 1];
		fa.getline(buffer, 10000);
		// For UNIX system, the end for line might be '\n\r'
		if(buffer[strlen(buffer) - 1] == '\r'){
			buffer[strlen(buffer) - 1] = '\0';
		}
		string tmp = buffer;
		if(buffer[0] == '>'){
			if( sl.name.compare("") ) {
				fasta.push_back(sl);
			}
			sl.name = tmp.substr(1);
			sl.seq = "";
		}else{
			sl.seq += tmp;
		}
	}
	fasta.push_back(sl);
	fa.close();
	return fasta.size();
}

int ReadFasta( char* infile, vector<seqline> & fasta  )
{
	return ReadFasta( string(infile), fasta );
}

/** Output the structure in a FASTA-format file. By default, the length of each
 *  line is 50 bp.
 */
int OutputFasta( string filename, vector<seqline> & fasta, int len = 50 )
{
	ofstream fa(filename.c_str());
	if(!fa) {
		cerr << "Cannot create fasta file " << filename << "\n";
		return -1;
	}
	vector<seqline>::iterator iter;
	for( iter = fasta.begin(); iter != fasta.end(); iter++ ) {
		fa << ">" << iter->name << endl;
		string seq = iter->seq;
		while( seq.length() > len ) {
			fa << seq.substr(0,len) << endl;
			seq = seq.substr(len);
		}
		fa << seq << endl;
	}
	fa.close();
	return fasta.size();
}

int OutputFasta( char* infile, vector<seqline> & fasta  )
{
	return OutputFasta( string(infile), fasta );
}

#endif
