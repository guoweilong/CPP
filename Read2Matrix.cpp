/* Usage: Read2Matrix -f <Reads file> -o <matrix file>
 * Detail:
 *    This program will get the weight matrix from a file,
 *    each line of which can be a sequence. Gap (or other 
 *    characters) will not be counted.
 * Author: GUO Weilong
 * 2010-12-10
 */
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
using namespace std;

struct parameter
{
	string infile;	// -s
	string outfile;	// -o
};

parameter param;

void exit_with_help( void )
{
	printf(
	"Usage: Read2Matrix -f <Reads file> -o <matrix file>\n"
	"options:\n"
	"\t-f Source file : file in fasta format\n"
	"\t-o Outfile : file to store result\n"
	);
	exit(1);
}

void ToUpperString(string &str)
{
    transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	// 
	int i;
	// parse options
	for(i=2;i<argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
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
	if ( !param.infile.length() || !param.outfile.length() ){
		exit_with_help();

      }
}

void init() {
	// default values
	param.infile = "";
	param.outfile = "";
}

int main(int argc, char* argv[])
{
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.outfile = string("output.txt");
#else
	parse_command_line(argc,argv);
#endif
	ifstream infile(param.infile.c_str());
      if (!infile ) {
            cerr << "Can't open file" << param.infile << "!" << endl;
            return -1;
      }
      char buffer[100];
      infile.getline(buffer, 100);
      int len = strlen(buffer);
      int *A, *C, *G, *T;
      A = new int[len]; C = new int[len];  G = new int[len]; T = new int[len];
      int i;
      for (i = 0; i < len; i++) {
            A[i] = C[i] = G[i] = T[i] = 0;
            switch (buffer[i]) {
            case 'a':
            case 'A':
                  A[i]++; break;
            case 'c':
            case 'C':
                  C[i]++; break;
            case 'g':
            case 'G':
                  G[i]++; break;
            case 't':
            case 'T':
                  T[i]++; break;
            default: break;
            }
      }
      while (!infile.eof()) {
            infile.getline(buffer, 100);
            //cout << buffer << endl;
            for (i = 0; i < len; i++) {
                  switch (buffer[i]) {
                  case 'a':
                  case 'A':
                        A[i]++; break;
                  case 'c':
                  case 'C':
                        C[i]++; break;
                  case 'g':
                  case 'G':
                        G[i]++; break;
                  case 't':
                  case 'T':
                        T[i]++; break;
                  default: break;
                  }
            }
       }
      infile.close();
      
      ofstream of(param.outfile.c_str());
      if(!of) {
            cerr << "Can't create the file" << param.outfile << endl; return -1;
      }
      of << "A [ ";  for (i = 0; i < len; i++) { of << A[i] << " "; }; of << "]\n";
      of << "C [ ";  for (i = 0; i < len; i++) { of << C[i] << " "; }; of << "]\n";
      of << "G [ ";  for (i = 0; i < len; i++) { of << G[i] << " "; }; of << "]\n";
      of << "T [ ";  for (i = 0; i < len; i++) { of << T[i] << " "; }; of << "]\n";
      of.close();

      return 1;
}


