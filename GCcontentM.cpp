#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	int i;
	if (argc < 2){
		cout << "Too few parameters.\n";
		return 0;
	}
	long A,C,G,T;
	A = C = G = T = 0;
	ifstream infile;
	infile.open(argv[1]);
	char line[100];
	while(infile.getline( line, 100, '\n')){
		if( line[0] == '>' ){
                  C--;
			continue;
		}
		for ( i = 0; i < 100 && line[i]; i++ ){
			switch( line[i] ){
				case 'A':   A++;  break;
				case 'a':   A++;  break;
				case 'C':   C++;  break;
				case 'c':   C++;  break;
				case 'G':   G++;  break;
				case 'g':   G++;  break;
				case 'T':   T++;  break;
				case 't':   T++;  break;
				default:	break;
			}
		}
	}
	double total = A + C + G + T;
	cout << "A:\t" << A/total << endl << "C:\t" << C/total << endl << "G:\t" << G/total << endl << "T:\t" << T/total << endl;
	infile.close();
	return 1;
}

