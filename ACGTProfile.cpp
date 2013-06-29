/* ACGTProfile 400 <inputfile> <output_file>
 * count the number that A,C,G,T appear on each position
 * Author: Weilong Guo
 * last modification: 2010/08/08
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[])
{
	int i;
	int section = 1;
	if (argc < 4){
		cout << "Usage: ACGTProfile [reads_length] <input.fa> <result.txt> [option:sections number]" << endl;
		cout << "Author: Guo, Weilong   Last Modified: 2011-06-14" << endl;
		cout << "Options:" << endl;
		cout << "   reads_length    <int>   Length of the whole flanking region" << endl;
		cout << "   sections number <int>   Length of bin to get the profile" << endl;
		cout << "Example: " << endl;
		cout << "  \t ACGTProfile 100 11H1R1_100chr1+.fa result.txt"	<< endl;
		cout << "or\t ACGTProfile 100 11H1R1_100chr1+.fa result.txt 10"	<< endl;
		return 0;
	}
	int len = atoi(argv[1]);
	if(argc == 5) {
		section = atoi( argv[4] );
	}

	ifstream infile;
	infile.open(argv[2]);
	ofstream outfile;
	outfile.open(argv[3]);
	
	// Initiation
	int *A, *C, *G, *T;
	A = new int[len];
	C = new int[len];
	G = new int[len];
	T = new int[len];
	for ( i = 0; i < len; i++ ) {
		A[i] = C[i] = G[i] = T[i] = 0;
	}
	
	string seq = "";
	string tmp;

        while(!infile.eof()){
                char buffer[10000 + 1];
                infile.getline(buffer, 10000);
                // For UNIX system, the end for line might be '\n\r'
                if(buffer[strlen(buffer) - 1] == '\r'){
                        buffer[strlen(buffer) - 1] = '\0';
                }
                tmp = buffer;
		if( tmp[0] == '>' ){
			for ( i = 0; i < len && seq[i]; i++ ){
				switch( seq[i] ){
				case 'A': 
				case 'a':
					A[i]++; break;
				case 'C':
				case 'c':
					C[i]++; break;
				case 'G':
				case 'g':
					G[i]++; break;
				case 'T':
				case 't':
					T[i]++; break;
				default:
					break;
				}
			}
			seq = "";
		} else {
			seq += tmp;
		}
	}
	// After while, for the last sequence
	for ( i = 0; i < len && seq[i]; i++ ){
		switch( seq[i] ){
		case 'A': 
		case 'a':
			A[i]++; break;
		case 'C':
		case 'c':
			C[i]++; break;
		case 'G':
		case 'g':
			G[i]++; break;
		case 'T':
		case 't':
			T[i]++; break;
		default:
			break;
		}
	}

	infile.close();
	outfile << "P0\tA\tC\tG\tT" << endl;
	int a, c, g, t;
	for ( i = 0; i < len/section; i++ ) {
		a = c = g = t = 0;
		for ( int j = 0; j < section; j++ ){
			a += A[ section * i + j ];
			c += C[ section * i + j ];
			g += G[ section * i + j ];
			t += T[ section * i + j ];
		}
		outfile << i*section << "\t" << a << "\t" << c << "\t" << g << "\t" << t << endl;
	}
	outfile.close();
	return 1;
}

