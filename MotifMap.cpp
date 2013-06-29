/* Usage: MotifMap -s <filename.fa> -M <motif.txt> -m 1 -c 50 -A -o <output.fa> -O 3 -H 1 -p <pos.txt>
 * Weilong GUO 2010-08-01
 * *********************
 * Modification:
 * 2010-08-10: Change the model to calculate the weight of motif binding site
 *     as one count for each site of every hit
 * 2010-08-26: Combine the above two model together to a parameter let user to choose
 * 2010-09-05: Modify the behavour of the output, output the mapped result as FASTA format
 * 2010-09-06: Fix a bug caused by 0905 modification, which cofused size_type(unsigned) with int
 * 2010-09-13: Modify the mechanism of indentifying the mapped reads considering the overlaps
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
#include "math.h"

struct parameter
{
	int m;			// -m
	int c;			// -c
	int H;			// -H
	int Overlap;	// -O
	string motif;	// -M
	string infile;	// -s
	string outfile;	// -o
	string pos;		// -p
	bool Anti;	// -A Antisense Check
};

struct Seq{
	string name;
	string content;
	string status;
};

struct Tread{
	int start;
	int length;
};

parameter param;
vector<string> MotifVec;	// Vector to store the motifs

void exit_with_help( void )
{
	printf(
		"Usage: MotifMap -s <filename.fa> -M <motif.txt> -m 1 -c 50 -A -o <output.fa> -O 3 -H 1 -p <pos.txt>\n"
		"options:\n"
		"-s Source file : Fasta file to be read in\n"
		"-M Motif file : Read in as motif, each line is a possible word\n"
		"-m <int> : to specify the mismatch number\n"
		"-c <int> : to specify the center of the fasta sequence\n"
		"-O <int> : to specify the overlap that allowed for different sequence\n"
		"-A : The Antisense of the motif will also be checked\n"
		"-o outfile : Output all the mapped sequences in FASTA format\n"
		"-H <int> : sepecify the model to accumulate the hit site; \n"
		"           0 for only center[default setting], 1 for every site\n"
		"-p : File to store the positions\n"
	);
	exit(1);
}

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	// MotifMap -s <filename.fa> -M <motif.txt> -m 1 -c 50 -A -o <output.fa> -O 3 -H 1 -p <pos.txt>
	int i;

	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 'M':
				param.motif = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 's':
				param.infile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'p':
				param.pos = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'A':
				param.Anti = true;
				break;
			case 'm':
				param.m = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'c':
				param.c = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'H':
				param.H = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'O':
				param.Overlap = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.infile.length() || !param.motif.length() ){
		exit_with_help();
	}
}

int ReadSequence ( string filename, vector<Seq> & SeqVec )
{
	ifstream infile(filename.c_str());
	if(!infile)
		return 0;
	string str;
	vector<Seq>::iterator iter = SeqVec.begin();
	Seq seq;
	seq.name = "";
	while ( getline( infile, str ) ) {
		if( str[0] == '>' ) {
			if( !seq.name.empty() ){
				SeqVec.push_back(seq);
			}
			seq.name = str.substr(1,str.length());
			seq.content = "";
			seq.status = "";
		} else {
			ToUpperString(str); //only toupper the sequence line
			seq.content += str;
			seq.status += str;
		}
	}
	
	if( !seq.name.empty() ){
		SeqVec.push_back( seq );
	} // push back the last element
	infile.close();
	return 1;
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

int ReadMotif ( string filename, vector<string> & MotifVec )
{
	ifstream infile(filename.c_str());
	if(!infile)
		return 0;
	string str;
	while ( getline( infile, str ) ) {
		vector<string> str_vec = string_tokenize(str);
		str = str_vec[0];
		ToUpperString( str );
		MotifVec.push_back ( str );
	}
	infile.close();
	return 1;
}

string AntiSense( string &seq ) {
	string AntiSeq = seq;
	for ( string::size_type i = 0; i < AntiSeq.length(); i++ ) {
		switch( seq[i] ) {
		case 'A':
		case 'a':
			AntiSeq[i] = 'T';
			break;
		case 'C':
		case 'c':
			AntiSeq[i] = 'G';
			break;
		case 'G':
		case 'g':
			AntiSeq[i] = 'C';
			break;
		case 'T':
		case 't':
			AntiSeq[i] = 'A';
			break;
		}
	}
	reverse(AntiSeq.begin(), AntiSeq.end());
	return AntiSeq;
}

bool SameSeq ( string strA, string strB) {
	if ( strA.size() != strB.size() )
		return false;
	int n = (int) strA.size();
	for ( int i = 0; i < n; i++ ) {
		if( strA[i] != '-' && strA[i] != 'N' && strB[i] != '-' && strB[i] != 'N'
			&& strA[i] != strB[i])
			return false;
	}
	return true;
}

int MapBack ( vector<Seq> &SeqVec ){
	ofstream of(param.outfile.c_str());
	vector<Seq>::iterator iter;
	string tmp = "";
	bool mode = false;
	bool p  = false; // Indicate whether to output the position information
	if( param.pos.length() != 0 ) p = true; 
	int start = 0;
	ofstream pos;
	if(p) {
		pos.open( param.pos.c_str() );
	}
	vector<string>::iterator MV_end;
	MV_end = unique( MotifVec.begin(), MotifVec.end() );

	for ( iter = SeqVec.begin(); iter != SeqVec.end(); iter++ ) {
		int n = (int)iter->content.size();
		int i = 0, end = 0;
		vector<Tread> Reads;
		Tread rd;

		bool state = false;	// whether in matched
		while ( i < n ) {
			if( state ) {
				if( rd.start + rd.length - param.Overlap < i ) {
					Reads.push_back(rd);
					state = false;
				} 
			}
			bool matched = false;
			int end = 0;
			for ( vector<string>::iterator Miter = MotifVec.begin(); Miter != MV_end; Miter++ ) {
				int ml = (int) Miter->size();
				if ( i + ml <= n && SameSeq(iter->content.substr(i,ml), (*Miter)) ) {
					matched = true;
					end = ( end < i + ml ) ? ( i + ml ) : end;
				} 
			} // for the motif vector
			if ( matched && !state ) {
				state = true;
				rd.start = i; rd.length = end - i;
			}
			i++;
		}
		for ( vector<Tread>::iterator riter = Reads.begin(); riter != Reads.end(); riter++ ) {
			of << ">" << iter->name << ":" << riter->start << endl;
			of << iter->content.substr(riter->start, riter->length) << endl;
			if( param.H == 0 ) {
				pos << riter->start + riter->length/2 -param.c << endl;
			} else if ( param.H == 1) {
				for ( int k = 0; k < riter->length; k++ )
					pos << riter->start + k -param.c << endl;
			}
		}
	}
	if(p) pos.close();
	return 0;
}

void init(){
	// default values
	param.motif="";
	param.infile="";
	param.outfile="";
	param.Anti = false;
	param.m	= 1;
	param.c	= 0;
	param.pos  = "";
	param.H	= 0;
	param.Overlap = 3;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.motif = string("motif.txt");
	param.Anti = true;
	param.outfile = string("output.txt");
	param.m	= 0;
	param.c	= 50;
	param.pos  = string("pos.txt");
	param.H = 1;
	param.Overlap = 3;
#else
	parse_command_line(argc,argv);
#endif
	vector<Seq> SeqVec;	// Vector to store the sequences
	if( !ReadSequence( param.infile, SeqVec ) ){
		cout << "Reading sequence file error.\n";
		return 0;
	}
	if( !ReadMotif( param.motif, MotifVec ) ){
		cout << "Reading motif file error.\n";
		return 0;
	}
	MapBack( SeqVec );
	return 1;
}
