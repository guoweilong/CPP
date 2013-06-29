/* Usage:CountKmers -s <foreground.fa> -b <background.fa> -K 4 -r -o [output_file]
 * Weilong GUO 2010-06-19
 * ***********************
 * add the parameter for background
 */
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
#include "math.h"

struct parameter
{
	int K;
	string infile;
	string backfile;
	string outfile;
	bool RC; //RecursiveCheck
};

parameter param;


void exit_with_help( void )
{
	printf(
	"Usage: CountKmers -s <foreground.fa> -b <background.fa> -K 4 -r -o [output_file]\n"
	"options:\n"
	"-s Source file : forground file in fasta format\n"
	"-b Background file : file in fasta format\n"
	"-o outfile : file to store result\n"
	"-K <int> : to specify the length to evaluate\n"
	"-r : when specified, the program will run in a model considering\n"
	"	the recursive case, for example:\n"
	"		AAAAAAAAAAAAAAAAAA\n"
	"	will only be recognize as once for 6-mer \"AAAAAA\"\n"
	);
	exit(1);
}

void parse_command_line(int argc, char **argv)
{
	int i;

	// default values
	param.infile.clear();
	param.K	= 6;
	param.outfile.clear();
	param.backfile.clear();
	param.RC = false;

	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 's':
				param.infile = string(argv[i]);
				break;
			case 'b':
				param.backfile = string(argv[i]);
				break;
			case 'o':
				param.outfile = string(argv[i]);
				break;
			case 'r':
				param.RC = true;
				i--;
				break;
			case 'K':
				param.K = atoi(argv[i]);
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if( !param.infile.length() )
        exit_with_help();
}

int CountFromFile( string filename, int* &Count, int K )
{
	int i;
	ifstream infile;

	infile.open(filename.c_str());
	long MAX = (long) pow( 4.0, K );
	Count = new int[MAX];
// Initialization for counting result of each pattern
	for( i = 0; i < MAX; i++ ) 
		Count[i]=0;
	char ch;
	int pos = 0;
	int v = 0;
	char temp[100];
	int* Record = new int[K];	// used when considering the recursive sequence
    while( infile.get(ch) ){
		switch ( ch ) {
		case '>':
			infile.getline( temp, 100 );
			for ( i = 0; i < K; i++)
				Record[i] = -1;
			pos = 0;
			v = 0;
			continue;
			break;
		case 'A':
		case 'a':
			pos++;
			v = ( v * 4 + 0 ) % MAX;
			break;
		case 'C':
		case 'c':
			pos++;
			v = ( v * 4 + 1 ) % MAX;
			break;
		case 'G':
		case 'g':
			pos++;
			v = ( v * 4 + 2 ) % MAX;
			break;
		case 'T':
		case 't':
			pos++;
			v = ( v * 4 + 3 ) % MAX;
			break;
		default:
			continue;
			break;
		}
		if( pos > K - 1 ) {
			if(!param.RC) {
				Count[v]++;
			} else {
				int j;
				for ( j = 0; j < K && v != Record[j]; j++ );
				if ( j == K )
					Count[v]++;
				Record[pos%K] = v;
			}
		}
	}
	infile.close();
	return 1;
}

int main(int argc, char* argv[])
{
	int *Count, *BCount;
	int n;
#ifdef _DEBUG
	param.K = 4;
	param.infile = string("input.fa");
	param.backfile = string("background.fa");
	param.RC = true;
#else
	parse_command_line(argc,argv);
#endif

	// Read the foreground file
	CountFromFile( param.infile, Count, param.K );

	// Read the background file
	bool AgainstBack = false;
	if ( param.backfile.length() )
		AgainstBack = true;
	if ( AgainstBack ) {
		CountFromFile( param.backfile, BCount, param.K );
	}

	// Output the results
	char* read = new char[param.K+1] ;
	read[param.K] = '\0';
		long MAX = (long) pow( 4.0, param.K );

	
	for( int i = 0; i < MAX; i++ ){
		if( Count[i] > 0 ) {
			int n = i;
			for( int j = param.K - 1; j >= 0; j-- ){
				switch( n % 4 ) {
				case 0:	read[j] = 'A'; break;
				case 1:	read[j] = 'C'; break;
				case 2:	read[j] = 'G'; break;
				case 3:	read[j] = 'T'; break;
				default: cerr << "Error! " << endl;
				}
				n /= 4;
			}
			if ( AgainstBack ) {
				cout << setprecision(4) << read << "\t" << (double)Count[i] / BCount[i] << "\t" << Count[i] <<"\t" << BCount[i]  << endl;
			} else {
				cout << read << "\t" << Count[i] << endl;
			} 
		}
	}
	return 1;
}

