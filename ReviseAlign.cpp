/* Usage: ReviseAlign -f <AlignmentResult> -o <output.txt>
 * Weilong GUO 2010-09-12
 * *********************
 * Description:
 *     Revise the result of alignment, cutting the sites with
 *     weak signal from both sides.
 * Modification:
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

struct parameter
{
	string infile;	// -f
	string outfile;	// -o
};
parameter param;

void exit_with_help( void )
{
	printf(
		"Usage:  ReviseAlign -f <AlignmentResult> -o <output.txt>\n"
		"options:\n"
		"-f Source file : Alignment resulted file to be read in\n"
		"-o outfile : write the result\n"
	);
	exit(1);
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
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'f':
				param.infile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.infile.length() || !param.outfile.length() ){
		exit_with_help();
	}
}

void init(){
	// default values
	param.infile="";
	param.outfile="";
}
// Definition of the arrays to store the countings
int * A, *C, *G, *T;

int Count(int index, char ch ) {
	switch(ch) {
	case 'A':
	case 'a':
		A[index]++; 
		break;
	case 'C':
	case 'c':
		C[index]++; 
		break;
	case 'G':
	case 'g':
		G[index]++; 
		break;
	case 'T':
	case 't':
		T[index]++; 
		break;
	default:
		break;
	}
	return 0;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = string("Align.txt");
	param.outfile = string("output.txt");
#else
	parse_command_line(argc,argv);
#endif

	double count = 0;
	const double BOTTOM = 0.25;

	ifstream infile;
	infile.open(param.infile.c_str());
	if( infile.fail() ){ 
		cerr << "Input file can not be found.\n";
		return 3;
	}
	vector<string> vstr;
	char ch[100];

	infile.getline(ch, 100);	count++;
	string str(ch); vstr.push_back(str);
	// Dydamicly alloc the spaces
	const int len = (int) str.size();
	A = new int[len];
	C = new int[len];
	G = new int[len];
	T = new int[len];
	char * state = new char[len];
	// Initiation the arrays
	for( int i=0; i<len; i++ ) {
		A[i] = C[i] = G[i] = T[i] = 0; 
	}
	for( int i=0; i<len; i++ ) {
		 Count(i,str[i]);
	}

	while ( infile.getline(ch, 100) ) {
		count++;
		str = string(ch); vstr.push_back(str);
		if( (int)str.size() != len ) {
			cerr << "Lengths are not the same..\n";
			return 2;
		}
		for( int i=0; i<len; i++ ) {
			 Count(i,str[i]);
		}
	}
	infile.close();
	// valuate the sites
	for( int i=0; i<len; i++ ) {
		char s = '-';
		if( A[i]+C[i]+G[i]+T[i] <= BOTTOM * count ) {
			state[i] = s;
		} else { 
			// Compare the counting results among the four kinds of nucleotides
			if( A[i] >= C[i] && A[i] >= G[i] && A[i] >= T[i] ) {
				s = 'A';
			} else if ( C[i] >= A[i] && C[i] >= G[i] && C[i] >= T[i] ) {
				s = 'C';
			} else if ( G[i] >= A[i] && G[i] >= C[i] && G[i] >= T[i] ) {
				s = 'G';
			} else {
				s = 'T';
			}
			state[i] = s;
		}
	}
	// Get the start and end positions which are meaningful
	int k;
	for ( k = 0; k < len && state[k] == '-' ; k++ ) ;
	int head = k;
	for ( k = len -1; k >= 0 && state[k] == '-' ; k-- ) ;
	int tail = k;

	ofstream outfile( param.outfile.c_str() );
	for ( vector<string>::iterator iter =  vstr.begin(); iter != vstr.end(); iter++ ) {
		outfile << iter->substr(head, tail-head+1) << endl;
	}
	outfile.close();
	// Output the pattern that have been recognized.
	cout << string(state).substr(head, tail-head+1) << endl;

/*** for test  ***
	cout << "\tA\tC\tG\tT\t*\n";
	for( int i = 0; i < len; i++ ) {
		 cout << i << "\t" << A[i] << "\t" << C[i] << "\t" << G[i] << "\t" << T[i] << "\t" << state[i] << endl;
	}
***/
	return 0;
}


