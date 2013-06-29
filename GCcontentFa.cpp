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
			double total = A + C + G + T;
			if(total)
				cout << A/total << "\t" << C/total << "\t" << G/total << "\t" << T/total << "\t" << (A+G)/total << endl;
			A = C = G = T = 0;
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
	if(total)
		cout << A/total << C/total << G/total << T/total << (A+G)/total << endl;
	infile.close();
	return 1;
}

