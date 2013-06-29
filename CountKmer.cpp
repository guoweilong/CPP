/* Usage: CountKmers -s <filename.fa> -K 4 -g 5 -r -i -A[nti] -o [output_file]
 * Weilong GUO 2010-07-24
 * -----------------------
 * Modified on 2010-07-25 add parameter -i -r
 * still can not be transfer to Unix directly
 */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
#include "math.h"
#include "Fasta.h"

// In order to be compatible with Linux and Windows
#if defined(__GNUC__)
# if __GNUC__ < 3 && __GNUC__ >= 2 && __GNUC_MINOR__ >= 95
#     include <hash_map>
# elif __GNUC__ == 3
#     include <ext/hash_map>
      using namespace __gnu_cxx;
# elif __GNUC__ >= 4
#include <ext/hash_map>
using namespace __gnu_cxx;
namespace __gnu_cxx
{
  template<> struct hash<string>
  {
	size_t operator()(const string& s) const
	{
	  return __stl_hash_string(s.c_str());
	}
  };
}
# endif
#elif defined(__MSVC_VER__)
# if __MSVC_VER__ >= 7
#     include <hash_map>
# else
#     error "std::hash_map is not available with this compiler"
# endif
#elif defined(__sgi__)
# include <hash_map>
#else
	#include <hash_set>
	#include <hash_map>
	using namespace stdext;
#endif 

struct parameter
{
	string infile;
	string outfile;
	int K;
	int g;
	bool RC; //RecursiveCheck
	bool i;  // whether to give the progress information
	bool Anti;
};

static parameter param;

void exit_with_help( void )
{
	printf(
		"Usage: CountKmers -s <filename.fa> -K 4 -g 5 -r -i -A[nti] -o [output_file]\n"
		"options:\n"
		"-s Source file : file in fasta format\n"
		"-o outfile : file to store result\n"
		"-A/-Anti : to indicate the antisense is equal to sense\n"
		"-K <int> : to specify the length to evaluate\n"
		"-g <int> : to specify the max length of gaps that are allowed in the middle\n"
		"-i <int> : to tell the progress of this program\n"
		"-r : when specified, the program will run in a model considering\n"
		"    the recursive case, for example:\n"
		"        AAAAAAAAAAAAAAAAAA\n"
		"    will only be recognize as once for 6-mer \"AAAAAA\"\n"
	);
	exit(1);
}

void parse_command_line(int argc, char **argv)
{
	int i;
	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 's':
				param.infile = argv[i];
				break;
			case 'o':
				param.outfile = argv[i];
				break;
			case 'r':
				param.RC = true;
				i--;
				break;
			case 'i':
				param.i = true;
				i--;
				break;
			case 'A':
				param.Anti = true;
				i--;
				break;
			case 'K':
				param.K = atoi(argv[i]);
				break;
			case 'g':
				param.g = atoi(argv[i]);
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if( param.infile.empty() && param.outfile.empty() )
        exit_with_help();
}

void ToUpper(string &str)
{
    transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);
}

int toanti(int ch){
	switch(ch){
	case 'A':
	case 'a':
		return 'T';
	case 'C':
	case 'c':
		return 'G';
	case 'G':
	case 'g':
		return 'C';
	case 'T':
	case 't':
		return 'A';
	default:
		return ch;
	}
}

/** Return the antisence of input sequence
 *  Example: AANNTTGGCC -> GGCCAANNTT
 */
string ToAntisense(string seq)
{
	string str = seq;
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toanti);
	reverse(str.begin(), str.end());
	return str;
}

void init() {
	// default values
	param.infile = "";
	param.outfile = "";
	param.K	= 4;
	param.g = 0;
	param.RC = false;
	param.i = false;
	param.Anti = false;
}

struct CountEle{
	string seq;
	int num;
};

// Used when doing the sorting of the CountEles
bool moreThan( const CountEle &a, const CountEle &b)
//bool moreThan(CountEle &a, CountEle &b)
//The second way complied wrong under g++
{
	if(a.num == b.num) {
		return (a.seq <= b.seq);
	} else {
		return (a.num > b.num);
	}
}

int main(int argc, char* argv[])
{
	//Initiation
	init();
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.outfile = string("out.txt");
	param.K = 16;
	param.g = 8;
	param.RC = true;
	param.i  = true;
	param.Anti = true;
#else
	parse_command_line(argc,argv);
#endif

	//Read in the FASTA file
	vector<seqline> SeqLine;
	ReadFasta(param.infile, SeqLine);
	if( param.i ) cout << "Finish reading input FASTA file.\n";

	//OutputFasta(param.outfile, SeqLine); // for test

	//Word counting
	hash_map<string, int>* HT = new hash_map<string,int>[param.g+1];
	string* last = new string[param.g+1];
	int k = param.K; int fh = (k+1)/2; int sh = k/2; // fh: first half; sh: second half 
	vector<seqline>::iterator iter;
	for( int j = 0; j <= param.g; j++ ) {
		if( param.i ) cout << "\rWord counting: "<< j << "/" << param.g << flush;
		for( iter = SeqLine.begin(); iter != SeqLine.end(); iter++ ) {
			last[j] = "";
			string seq = iter->seq;
			int n = seq.length();
			for( int s = 0; s < n-k-j; s++ ) {
				string kmer = seq.substr(s,fh)+seq.substr(s+fh+j,sh); ToUpper(kmer);
				if ( kmer.find('N') == string::npos ) {
					// If 'N' is found in the kmer, then the kmer is not valid
					if( param.Anti ){
						string Anti = ToAntisense(kmer);
						kmer = (kmer<=Anti)?kmer:Anti;
					}
					if( last[j] != kmer ) {
						HT[j][kmer] = ++HT[j][kmer];
						last[j] = kmer;
					}
				}
			}
		}
	}
	if( param.i ) cout << "\nFinish word counting.\n";


	//Store the result in a vector
	vector<CountEle> Result;
	Result.clear();
	for( int j = 0; j <= param.g; j++ ) {
		if( param.i ) printf("\rProcessing result: %d/%d", j, param.g);
		hash_map<string,int>::iterator hash_iter;
		for( hash_iter = HT[j].begin(); hash_iter != HT[j].end(); hash_iter++ ) {
			CountEle CE;
			string first = hash_iter->first; 
			CE.seq = first.substr(0,fh) + string(j,'N') + first.substr(fh,sh);
			CE.num = hash_iter->second;
			Result.push_back(CE);
		}
	}
	//Sort the result in order
	if( param.i ) printf("\nSorting the result...\n");
	sort(Result.begin(), Result.end(), moreThan);

	//Save the result in outfile
	ofstream of( param.outfile.c_str() );
	for( vector<CountEle>::iterator CE_Iter = Result.begin(); CE_Iter != Result.end(); CE_Iter++ ) {
		of << CE_Iter->seq << "\t" << CE_Iter->num << endl;
	}
	if( param.i ) cout << "Finish saving the result.\n";
	of.close();
	
	return 1;
}
