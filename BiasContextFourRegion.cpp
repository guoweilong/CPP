/* BiasContextFourRegion.cpp
 * 2011-12-01 From BiasContext.cpp 2011-10-27
 * Modification:
 * 
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <iomanip>

using namespace std;
//#include "math.h"

struct parameter
{
	string fastafile;	// -f
	string refFlatfile;	// -r
};

parameter param;

void exit_with_help( void )
{
	printf(
		"Usage:	BiasContextFourRegion -f <Fasta> -r <refFlat> \n"
		"Author: Guo, Weilong  guoweilong@gmail.com  2011-12-1\n"
		"Options:\n"
		"-f	Input, FASTA format, store sequence information\n"
		"-r	Input, refFlat file \n"
		"Output file, format \n"
		"	transID	exon[A,C,G,T] 5SS[A,C,G,T] MI[A,C,G,T] 3SS[A,C,G,T]\n"
		);
	exit(1);
}

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	int i;

	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
		case 'f':
			if(i == argc)	exit_with_help();
			param.fastafile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'r':
			if(i == argc)	exit_with_help();	
			param.refFlatfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.fastafile.length() || !param.refFlatfile.length() ){
		exit_with_help();
	}
}


vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
inline vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty) {
	// Skip delimiters at beginning.
	string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	vector<string> result;
	result.clear();

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		//__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");

		//if (pos == lastPos) result.push_back("");
		result.push_back(str.substr(lastPos, pos - lastPos));

		if (pos == string::npos) break;
		if (pos == str.length() - 1) {
			if (!skip_empty) result.push_back("");
			break;
		}
		// Skip delimiters.  Note the "not_of"
		lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
	return result;
}

void init(){
	// default values
	param.fastafile = "";
	param.refFlatfile = "";
}

#define MAX 10000

class RefFlat {
public:
	string 	gene;
	string 	trans;
	string 	chr;
	char	strand;
	unsigned int tStart;	// Start of transcript
	unsigned int tEnd;	// End of transcript
	unsigned int cStart;	// Start of coding region
	unsigned int cEnd;	// End of coding region
	short	nExon;		// Number of exons
	vector<unsigned int> eStart;	// Starts of exons
	vector<unsigned int> eEnd;	// Ends of exons
};

/*
 format of refFlat file
 DLX1	NM_178120	chr2	+	172658453	172662647	172658651	172661231	3	172658453,172659627,172660976,	172658964,172659827,172662647,
*/

int ReadRefFlatFile ( list<RefFlat> & Annotation ) {
	ifstream Rfile(param.refFlatfile.c_str());
	if(!Rfile) {
		cout << "cannot open input file" << param.refFlatfile.c_str() << endl;
		return -1;
	}
	while ( !Rfile.eof() ) {
		char buffer[MAX+1];
		Rfile.getline(buffer, MAX);
		if(!strlen(buffer)) continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);

		RefFlat rf;
		rf.gene = tokens[0];
		rf.trans = tokens[1];
		rf.chr = tokens[2];
		rf.strand = tokens[3][0];
		rf.tStart = atoi ( tokens[4].c_str() );
		rf.tEnd = atoi ( tokens[5].c_str() );
		rf.cStart = atoi ( tokens[6].c_str() );
		rf.cEnd = atoi ( tokens[7].c_str() );
		rf.nExon = atoi ( tokens[8].c_str() );
		vector<string> st = string_tokenize( tokens[9], "," );
		int n = st.size();
		for (int i = 0; i < n; i++ ) {
			rf.eStart.push_back ( atoi( st[i].c_str() ) );
		}
		vector<string> et = string_tokenize( tokens[10], "," );
		n = et.size();
		for (int i = 0; i < n; i++ ) {
			rf.eEnd.push_back ( atoi( et[i].c_str() ) );
		}
		Annotation.push_back(rf);
	}
	Rfile.close();
	
	return 1;
}


// To test function { ReadRefFlatFile }
int WriteRefFlatFile ( list<RefFlat> & Annotation ) {
	list<RefFlat>::iterator liter;
	for ( liter = Annotation.begin(); liter != Annotation.end(); liter++ ) {
		cout << liter->gene << "\t" << liter->trans << "\t" << liter->chr << "\t"
		     << liter->strand << "\t" << liter->tStart << "\t" << liter->tEnd << "\t"
		     << liter->cStart << "\t" << liter->cEnd << "\t"  << liter->nExon << "\t";
		int n = liter->eStart.size();
		for ( int i = 0; i < n; i++ ) {
			cout << liter->eStart[i] << ",";
		}
		cout << "\t"; n = liter->eEnd.size();
		for (int i = 0; i < n; i++ ) {
			cout << liter->eEnd[i] << ",";
		}
		cout << endl;
	}
	return 1;
}



/** Read in a file which is in FASTA format, and store in
 *  a 'map<string,string>' structure fasta.
 */
int ReadFasta( string filename, map<string,string> & fasta )
{
	fasta.clear();
	ifstream fa(filename.c_str());
	if(!fa) {
		cerr << "Cannot open fasta file " << filename << "\n";
		return -1;
	}
	string chr = "";
	string seq = "";
	while(!fa.eof()){
		char buffer[MAX + 1];
		fa.getline(buffer, MAX);
		// For UNIX system, the end for line might be '\n\r'
		if(buffer[strlen(buffer) - 1] == '\r'){
			buffer[strlen(buffer) - 1] = '\0';
		}
		string tmp = buffer;
		if(buffer[0] == '>'){
			if( chr.compare("") ) {
				fasta[chr]=seq;
			}
			chr = tmp.substr(1);
			seq = "";
		}else{
			seq += tmp;
		}
	}
	fasta[chr]=seq;
	fa.close();
	return fasta.size();
}

#define STEP 50

int WriteFasta( map<string,string> & fasta )
{
	map<string,string>::iterator fiter;
	for ( fiter = fasta.begin(); fiter != fasta.end(); fiter++ ) {
		string chr = fiter->first;
		cout << ">" << chr << endl;
		unsigned long len = fasta[chr].length();
		for (unsigned long i = 0; i < len; i+= STEP ) {
			cout << fasta[chr].substr(i,STEP) << endl ;
		}
	}
}

string GetSequence ( string chr, unsigned long start, int length, map<string, string> & Genome) {
	if ( length <= 0 || Genome[chr].length() < start + length)
		return string("");
	return Genome[chr].substr(start, length);
}

class Content {
public:
	float A;
	float C;
	float G;
	float T;
};

Content GetContent (string seq) {
	int len = seq.length();
	int A, C, G, T;
	A = C= G = T = 0;
	ToUpperString(seq);
	for (int i = 0; i < len; i++ ) {
		switch (seq[i]) {
		case 'A': A++; break;
		case 'C': C++; break;
		case 'G': G++; break;
		case 'T': T++; break;
		}
	}
	float all = A + C + G + T;
	Content ct;
	ct.A = A / all;
	ct.C = C / all;
	ct.G = G / all;
	ct.T = T / all;
	return ct;
}

int PrintIn50 (string seq) {
	int n = seq.length();
	for (int i = 0; i < n; i+=50 ) {
		cout << seq.substr(i,50) << endl;
	}
	return 1;
}

int BiasContext (map<string,string> & Genome, list<RefFlat> & Annotation) {
	list<RefFlat>::iterator liter;
	for ( liter = Annotation.begin(); liter != Annotation.end(); liter++ ) {
		int N = liter->nExon;
		string es, ls, ms, rs;
		es = ls = ms = rs = "";

		for ( int i = 1; i < N; i++ ) {
			int LE = liter->eStart[i-1]; // left of exon
			int LI = liter->eEnd[i-1]; // left of intron
			int RI = liter->eStart[i]; // right of intron
			int MI = ( LI + RI ) /2;
			
			es += GetSequence(liter->chr, LE, LI-LE, Genome);
			ls += GetSequence(liter->chr, LI+5, 120, Genome);
			ms += GetSequence(liter->chr, MI-60, 120, Genome);
			rs += GetSequence(liter->chr, RI-125, 120, Genome);

		}
		es += GetSequence(liter->chr, liter->eStart[N-1], liter->eEnd[N-1] - liter->eStart[N-1], Genome);
		Content ect = GetContent( es );
		Content lct = GetContent( ls );
		Content mct = GetContent( ms );
		Content rct = GetContent( rs );

		cout << setprecision(4);
		if (liter->strand == '+') {
			cout << liter->gene << "\t" << liter->trans; 
			cout  << "\t"<< ect.A << "\t" << ect.C << "\t" << ect.G << "\t" << ect.T;
			cout  << "\t"<< lct.A << "\t" << lct.C << "\t" << lct.G << "\t" << lct.T;
			cout  << "\t"<< mct.A << "\t" << mct.C << "\t" << mct.G << "\t" << mct.T;
			cout  << "\t"<< rct.A << "\t" << rct.C << "\t" << rct.G << "\t" << rct.T;
			cout << endl;
		} else {
			cout << liter->gene << "\t" << liter->trans; 
			cout  << "\t"<< ect.T << "\t" << ect.G << "\t" << ect.C << "\t" << ect.A;
			cout  << "\t"<< rct.T << "\t" << rct.G << "\t" << rct.C << "\t" << rct.A;
			cout  << "\t"<< mct.T << "\t" << mct.G << "\t" << mct.C << "\t" << mct.A;
			cout  << "\t"<< lct.T << "\t" << lct.G << "\t" << lct.C << "\t" << lct.A;
			cout << endl;
		}
	}
	return 1;
}

int main(int argc, char* argv[])
{
	init ();
	parse_command_line ( argc, argv );

	list<RefFlat> Annotation;
	ReadRefFlatFile ( Annotation );
//	WriteRefFlatFile ( Annotation );
	
	map<string,string> Genome;
	ReadFasta( param.fastafile, Genome );
//	WriteFasta ( Genome );

	BiasContext(Genome, Annotation);

	return 1;
}



