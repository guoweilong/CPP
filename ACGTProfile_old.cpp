/* ACGTProfile 400 <inputfile> <output_file>
 * count the number that A,C,G,T appear on each position
 */
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	int i;
	int section = 1;
	if (argc < 4){
		cout << "Too few parameters.\n";
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
	char line[100];
	int *A, *C, *G, *T;
	A = new int[len];
	C = new int[len];
	G = new int[len];
	T = new int[len];
	for ( i = 0; i < len; i++ ) {
		A[i] = C[i] = G[i] = T[i] = 0;
	}
	int index = 0;
	while(infile.getline( line, 100, '\n')){
		if( line[0] == '>' ){
			index = 0;
			continue;
		}
		for ( i = 0; i < 100 && line[i]; i++ ){
			switch( line[i] ){
				case 'A':   A[index]++; index++;  break;
				case 'a':   A[index]++; index++;  break;
				case 'C':   C[index]++; index++;  break;
				case 'c':   C[index]++; index++;  break;
				case 'G':   G[index]++; index++;  break;
				case 'g':   G[index]++; index++;  break;
				case 'T':   T[index]++; index++;  break;
				case 't':   T[index]++; index++;  break;
				default:	break;
			}
		}
	}
	infile.close();
	outfile << "position\tA\tC\tG\tT" << endl;
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

