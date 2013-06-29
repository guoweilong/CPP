/* Usage: FastaMask -f <filename.fa> -s 20 -l 25 -o <output_file>
 * Weilong GUO 2010-07-30
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
#include "Fasta.h"

struct parameter
{
	int s;		// start site
	int len;		// length to mask	
	string infile;	// -s
	string outfile;	// -o
};

parameter param;

void exit_with_help( void )
{
	printf(
	"Usage: FastaMask -f <filename.fa> -s 20 -l 5 -o <output_file>\n"
	"options:\n"
	"\t-f Source file : file in fasta format\n"
	"\t-o Outfile : file to store result\n"
	"\t-s  Start site : from which site to mask\n"
	"\t-l  length to be masked\n"
	);
	exit(1);
}

void ToUpperString(string &str)
{
    transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	// FastaMask -f <filename.fa> -s 20 -l 25 -o <output_file>
	int i;
	// parse options
	for(i=2;i<argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 's':
				param.s = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'l':
				param.len = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'f':
				param.infile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.infile.length() || !param.outfile.length()
		|| param.len < 0 || param.s < 0 ){
		exit_with_help();
	}
}

void init() {
	// default values
	param.infile = "";
	param.outfile = "";
	param.s = -1;
	param.len = 0;
}

int main(int argc, char* argv[])
{
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.outfile = string("output.txt");
	param.s	= 0;
	param.len = 10
#else
	parse_command_line(argc,argv);
#endif
	//Read in the FASTA file
	vector<seqline> SeqLine;
	ReadFasta(param.infile, SeqLine);
//	if( param.i ) cout << "Finish reading input FASTA file.\n";

	vector<seqline>::iterator iter;
	for( iter=SeqLine.begin(); iter!=SeqLine.end(); iter++){
		string seq = iter->seq; 
		iter->seq = seq.substr(0,param.s) + string(param.len,'N')
			+ seq.substr(param.s+param.len);
	}
	
	OutputFasta(param.outfile, SeqLine); 
	return 1;
}
