/* Usage: MotifCluster -M <motif.txt> -m 1 -A -o <output_file> -n <network.txt> -p <prefix>
 * Weilong GUO 2010-07-30
 * *********************
 * Modification:
 * 2010-08-01 Weilong GUO | add the parameter '-n'
 * 2010-08-06 Weilong GUO | modify the clustering algorithm
 * 2010-09-04 Weilong GUO | modify the OUT format, add the prefix to list motif clusters one-by-one
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

struct parameter
{
	int m;		// -m
	string motif;	// -M
	string outfile;	// -o
	string netfile;	// -n
	string prefix;	// -p
	bool Anti;	// -A Antisense Check
};

parameter param;

typedef struct 	//Element of Array
{
	int weight;	//weight at this element
	char dir; //'U' for up, 'L' for left, 'S' for slop. 
}ELE;

void exit_with_help ( void )
{
	printf(
		"Usage: MotifCluster -M <motif.txt> -m 1 -A -o <output.txt> -n <network.txt> -p <prefix>\n"
		"options:\n"
		"-M Motif file : first word in each line is considered to be a motif\n"
		"-o outfile : file to store result\n"
		"-m <int> : to specify the mismatch number\n"
		"-A : The Antisense of the motif will also be checked\n"
		"-n netfile : file to store the network of motifs\n"
		"-p prefix :  prefix of the filenames to show each motif cluster\n"
	);
	exit(1);
}

void ToUpperString (string &str)
{
    transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line (int argc, char **argv)
{
	// MotifCluster -M <motif.txt> -m 1 -A -o <output_file> -n <network.txt>
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
			case 'm':
				param.m = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'n':
				param.netfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'A':
				param.Anti = true;
				break;
			case 'p':
				param.prefix = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.motif.length() ){
		exit_with_help();
	}
}

vector<string> string_tokenize (const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
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

string AntiSense ( string &seq ) {
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

void init (){
	param.motif = "";
	param.outfile = "";
	param.netfile = "";
	param.Anti = false;
	param.m	= 1;
	param.prefix = "";
}

string Clean( string a ) {
	string str = a;
	string::size_type n = str.size();
	string::size_type i,j;
	for( i = 0, j = 0; j < n; j++ ) {
		if( str[j] != 'N' && str[j] != 'n' ) {
			str[i] = str[j];
			i++;
		}
	} // Move the non-N character forward
	str.resize(i); // resize the string before return
	return str;
}

int MotifDistance( string a, string b){
	string A = Clean(a);
	string B = Clean(b);
	string::size_type al = A.size()+1;
	string::size_type bl = B.size()+1;
	// Find the longest common string
	string::size_type i, j;
	ELE* Array = new ELE[al*bl];
	for( i = 0; i < al; i++) 
	{
		Array[i*bl].weight = 0;
		Array[i*bl].dir = '\0';
	}
	for( j = 0; j < bl; j++) 
	{
		Array[j].weight = 0;
		Array[j].dir = '\0';
	}
	// Mark the direction for each element
	for( i = 1; i < al; i++){
		for( j = 1; j < bl; j++){
			int a = Array[(i-1)*bl+j].weight;
			int b = Array[i*bl+(j-1)].weight;
			int c = (A[i-1] == B[j-1]) ? (Array[(i-1)*bl+(j-1)].weight+1) : -100000;
			if(c > a && c > b){
				Array[i*bl+j].weight = c;
				Array[i*bl+j].dir = 'S';
			}
			else if(a > b){
				Array[i*bl+j].weight = a;
				Array[i*bl+j].dir = 'U';
			} 
			else{
				Array[i*bl+j].weight = b;
				Array[i*bl+j].dir = 'L';
			}
		}
	}
	int k = Array[al*bl-1].weight; // The value of last cell is number of max matched sites
	return (int) ((((al-k)>=(bl-k)) ? (al-k) : (bl-k) )-1); // compare with the longer one is more strict 10/08/07
	//"al" and "bl" are 1 larger than the lengths of the two strings
}

// Static variables definition
char* state;	// An array to record the state of each motifs
char* mat;	// An matrix to record the
list<list<int> > Components; // A list to store the information of components // ">>" should be "> >"

// Recursive function to get the components
void FindComponent( vector<string>::size_type k, vector<string>::size_type l, list<int>& com ) {
	state[k] = '1';  // set state to '1' to label it as used
	com.push_back((int)k);
	for(vector<string>::size_type i = 0; i < l; i++ ) {
		if( mat[k*l+i] == '1' && state[i] == '0' ) {
			FindComponent(i,l,com);
		}
	}
}

void ComponentAnalysis( vector<string> &MotifVec ) {
	vector<string>::size_type len = MotifVec.size();
	// Initialization
	state = new char[len];
	memset( state, '0', len ); // Shouldn't be memset(state,'0',sizeof(state))
	for( vector<string>::size_type i=0; i<len; i++ ) {
		if( state[i] == '0' ){
			list<int> com;
			com.empty();
			FindComponent(i,len,com);
			Components.push_back(com);
		}
	}
}

int main (int argc, char* argv[])
{
	init();

#ifdef _DEBUG
	param.motif = string("13allmotif.txt");
	param.Anti = true;
	param.outfile = string("output.txt");
	param.netfile = string("network.txt");
	param.m	= 1;
	param.prefix = string("");
#else
	parse_command_line(argc,argv);
#endif

	// Get a vector of motifs
	vector<string> MotifVec;	// Vector to store the motifs
	ReadMotif ( param.motif, MotifVec );

	// Get the pairwise distance
	vector<string>::size_type l = MotifVec.size();
	mat = new char[l*l];
	memset(mat, '0', l*l);
	// Initialize a matrix, one dimension instead of two dimesion
	for( vector<string>::size_type i=0; i<l; i++ ) {
		for( vector<string>::size_type j=i+1; j<l; j++ ) {
			if( MotifDistance(MotifVec[i],MotifVec[j]) <= param.m ) {
				mat[i*l+j] = '1'; mat[j*l+i] = '1';
			} else if( param.Anti ) {
				if( MotifDistance(MotifVec[i],AntiSense(MotifVec[j])) <= 1 ) {
					mat[i*l+j] = '1'; mat[j*l+i] = '1';
				}
			}
		}
	}

	// Doing the component analysis
	ComponentAnalysis( MotifVec );

	// Output the list of the components
	ofstream of(param.outfile.c_str());
	int k = 0;
	int ks = 0; // Counter for seperate motifs
	cout << "Components size:" << Components.size() << endl;
	cout << "No.\tCount\n---\t-----\n";
	for( list<list<int> >::iterator iter = Components.begin(); iter != Components.end(); iter++ ) {
		cout << "#" << (k+1) << "\t" << iter->size() << endl;
		of << "Component #" << ++k << ":\n";
		for( list<int>::iterator liter = iter->begin(); liter != iter->end(); liter++ ) {
			of << "\t" << MotifVec[*liter] << endl;
		}
		if( param.prefix.size() &&  iter->size() >= 9 ) {
			ks++;
			stringstream fname; fname << param.prefix << ks << ".txt";
			ofstream sof(fname.str().c_str());
			for( list<int>::iterator liter = iter->begin(); liter != iter->end(); liter++ ) {
				sof << MotifVec[*liter] << endl;
			}
			sof.close();
		}
	}
	of.close();

	if(!param.netfile.length()) 
		return 1;
	// If netfile is specified, write network into the file
	of.open(param.netfile.c_str());
	for( vector<string>::size_type i=0; i<l; i++ ) {
		bool noedge = true;
		for( vector<string>::size_type j=i+1; j<l; j++ ) {
			if( mat[i*l+j]=='1' ) {
				noedge = false;

#ifdef _DEBUG
				of << MotifVec[i] << "\n" << MotifVec[j] << "\n\n"<< endl;
#else
				of << MotifVec[i] << "\t" << MotifVec[j] << endl;
#endif
			}
		}
		if (noedge){
				noedge = false;

#ifdef _DEBUG
				of << MotifVec[i] << "\n\n"<< endl;
#else
				of << MotifVec[i] << endl;
#endif
		}
	}
	of.close();

	return 1;
}

