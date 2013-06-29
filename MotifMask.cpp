/* Usage: MotifMask -s <filename.fa> -M <motif.txt> -m 1 -c _ -A -o <output_file>
 * Weilong GUO 2010-06-07
 * *********************
 * Modification:
 * 2010-07-30: add -c parameter; correct antisense function
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
	string motif;	// -M
	string infile;	// -s
	string outfile;	// -o
	char ch;
	bool Anti;	// -A Antisense Check
};

struct Seq{
	string name;
	string content;
	string status;
};

parameter param;

void exit_with_help( void )
{
	printf(
	"Usage: MotifMask -s <filename.fa> -M <motif.txt> -m 1 -c _ -A -o <output_file>\n"
	"options:\n"
	"-s Source file : file in fasta format\n"
	"-o outfile : file to store result\n"
	"-m <int> : to specify the mismatch number\n"
	"-A : The Antisense of the motif will also be checked\n"
	"-c : A char to take the place of the matched places\n"
	);
	exit(1);
}

void ToUpperString(string &str)
{
    transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	// MotifMask -s <filename.fa> -M <motif.txt> -m 1 -c _ -A -o <output_file>
	int i;

	for(i=2;i<argc;i++)
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
			case 'A':
				param.Anti = true;
				break;
			case 'm':
				param.m = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'c':
				param.ch = argv[i][0];
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
	for ( int i = 0; i < AntiSeq.length(); i++ ) {
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

int MaskSeq( vector<Seq> &SeqVec, string &motif ) {
	string temp;
	vector<Seq>::iterator iter;
	int l = motif.length();
	for ( iter = SeqVec.begin(); iter != SeqVec.end(); iter++ ) {
		int n = iter->content.length();
		bool mstate = false;
		string seq = iter->content;
		for ( int i = 0; i < n - l + 1; i++ ) {
			int mismatch = 0;
			int j = 0;
		//	cout << seq[i];
			while ( j < l && mismatch <= param.m ) {
				if( motif[j] != 'N' && motif[j]!= 'n') {
					if ( seq[i+j] != motif[j] )
						mismatch++;
				}
				j++;
			}
			if ( mismatch <= param.m ) {
				if ( mstate ) {
					iter->status[i+l-1] = param.ch;	// 需要减1,6mer的最后一位比第一位的下标大（l-1）
				} else {
					for ( int k = i; k < i + l; k++ )
						iter->status[k] = param.ch;
					mstate = true;
				}
			} else {
				mstate = false;
			}
		}
	}
	return 1;
}

int OutputSeq ( vector<Seq> &SeqVec ){
	ofstream outfile(param.outfile.c_str());
	vector<Seq>::iterator iter;
	for ( iter = SeqVec.begin(); iter != SeqVec.end(); iter++ ) {
		outfile << ">" ;
		outfile << iter->name << endl;
#ifdef _DEBUG
		//ToUpperString(iter->content);
		outfile << iter->content << endl;
#else
#endif
		//ToUpperString(iter->status);
		outfile << iter->status << endl;
	}
	outfile.close();
	return 1;
}

void init(){
	// default values
	param.motif="";
	param.infile="";
	param.outfile="";
	param.Anti = false;
	param.m	= 1;
	param.ch = 'N';
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.motif = string("motif.txt");
	param.Anti = true;
	param.outfile = string("output.txt");
	param.m	= 1;
	param.ch = 'N';
#else
	parse_command_line(argc,argv);
#endif
	vector<Seq> SeqVec;	// Vector to store the sequences
	ReadSequence( param.infile, SeqVec );

	vector<string> MotifVec;	// Vector to store the motifs
	ReadMotif ( param.motif, MotifVec );

	vector<string>::iterator MIter;	// iterator for the motifs
	for ( MIter = MotifVec.begin(); MIter != MotifVec.end(); MIter++ ) {
		int l = MIter->length();
		MaskSeq( SeqVec, *MIter );
		if( param.Anti ) {
			string AntiSeq = AntiSense( *MIter );
			MaskSeq( SeqVec, AntiSeq );
		}
	}
	OutputSeq( SeqVec );
	return 1;
}
