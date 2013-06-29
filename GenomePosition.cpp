/* Usage: GenomePosition -f <filename.fa> -s TACAG -c 2 -A -o <output.fa>
 * Weilong GUO 2010-08-09
 * *********************
 * Modification:
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;
#include "math.h"

struct parameter
{
	string motif;	// -s
	string infile;	// -f
	string outfile;	// -o
	bool Anti;		// -A Antisense Check
	int c;
};

struct Seq{
	string name;
	string content;
};

parameter param;

void exit_with_help( void )
{
	printf(
		"Usage: GenomePosition -f <filename.fa> -s TACAG -c 2 -A -o <output.fa>\n"
		"options:\n"
		"-s Source file : Fasta file to be read in\n"
		"-o outfile : write the result\n"
		"-c <int> : to specify the center of the fasta sequence,\n"
		"			EX: TACAG [2] for C\n"
		"-A : The Antisense of the motif will also be checked\n"
	);
	exit(1);
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

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);
}

void ToLowerString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))tolower);
}

map<string, string> chr2seq;
int GetGenome( const char* genomefile )
{
	ifstream geneNome(genomefile);
	if(!geneNome) {
		printf("cannot open geneNome file %s\n", genomefile);
		return -1;
	}
	cerr<<"reading genome ..."<<endl;
	string chr="";
	string chr_long_strings = "";
	while(!geneNome.eof()){
		char buffer[10000 + 1];
		geneNome.getline(buffer, 10000);
		if(buffer[strlen(buffer) - 1] == '\r'){
			buffer[strlen(buffer) - 1] = '\0';
		}
		string tmp = buffer;
		if(buffer[0] == '>'){
			if(chr == ""){
				chr_long_strings = "";
				vector<string> tokens = string_tokenize(tmp.substr(1));
				chr = tokens[0];	//specify which chromosome it is
				continue;
			}
			chr2seq[chr] = chr_long_strings;
			chr_long_strings = "";
			vector<string> tokens = string_tokenize(tmp.substr(1));
			chr = tokens[0];
		}else{
			chr_long_strings += tmp;
		}
	}
	chr2seq[chr] = chr_long_strings;
	// store the sequence in chr2seq[chr]
	cerr<<"reading genome done"<<endl;
	return 1;
}

int ReadSequence ( string filename, vector<Seq> & SeqVec )
{
	ifstream infile(filename.c_str());
	string str;
	vector<Seq>::iterator iter = SeqVec.begin();
	Seq seq;
	seq.name = "";
	cout << "Reading the FASTA file...\n";
	while ( getline( infile, str ) ) {
		if( str[0] == '>' ) {
			if( !seq.name.empty() ){
				SeqVec.push_back(seq);
			}
			seq.name = str.substr(1,str.length());
			seq.content = "";
		} else {
			ToUpperString(str); //only toupper the sequence line
			seq.content += str;
		}
	}
	
	if( !seq.name.empty() ){
		SeqVec.push_back( seq );
	} // push back the last element
	infile.close();
	cout << "Finished reading fasta file\n";
	return 1;
}

void parse_command_line(int argc, char **argv)
{
	// GenomePosition -f <filename.fa> -s TACAG -A -o <output.fa>
	int i;

	for(i=2;i<argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 's':
				param.motif = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'A':
				param.Anti = true;
				break;
			case 'f':
				param.infile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'c':
				param.c = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.infile.length() || !param.motif.length() || !param.outfile.length() ){
		exit_with_help();
	}
}


string Antisense( string &seq ) {
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

int ScanSeq( vector<Seq> &SeqVec ) {
	string temp;
	vector<Seq>::iterator iter;
	vector<Seq>::size_type l = param.motif.length();
	string motif = param.motif; ToUpperString(motif);
	string anti = Antisense(motif); ToUpperString(anti);
	ofstream of(param.outfile.c_str());
	if( !of ) {
		cerr << "Open output file error\n";
		exit(-1);
	}
	for ( iter = SeqVec.begin(); iter != SeqVec.end(); iter++ ) {
		vector<Seq>::size_type n = iter->content.length();
	//	ToUpperString(iter->content);
		bool mstate = false;
		for ( vector<Seq>::size_type i = 0; i < n - l + 1; i++ ) {
			string sub = iter->content.substr(i,l);
			if( sub == motif ) {
				of << iter->name << "\t" << i+param.c << "\t+\t" << motif << endl;
#ifdef _DEBUG
				of << iter->content << endl;
#endif
			} else if ( param.Anti && sub == anti) {
				of << iter->name << "\t" << i+l-1-param.c << "\t-\t" << motif << endl;
#ifdef _DEBUG
				of << iter->content << endl;
#endif
			}
		}
	}
	of.close();
	return 1;
}


void init(){
	// default values
	param.motif="";
	param.infile="";
	param.outfile="";
	param.Anti = false;
	param.c = 0;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.motif = string("TACAG");
	param.Anti = true;
	param.outfile = string("output.txt");
	param.c	= 2;
#else
	parse_command_line(argc,argv);
#endif
	//GetGenome( param.infile.c_str() );

	vector<Seq> SeqVec;	// Vector to store the sequences
	ReadSequence( param.infile, SeqVec );

	ScanSeq( SeqVec );

	return 1;
}

