#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	int i;
	if (argc < 3){
		cout << "Too few parameters.\n";
            cout << "Antisense input.fa output.fa\n";
		return 0;
	}
	ifstream infile;
	infile.open(argv[1]);
	ofstream outfile;
	outfile.open(argv[2]);
	char line[100];
	while(infile.getline( line, 100, '\n')){
		if( line[0] == '>' ){
			outfile << line << endl;
			continue;
		}
		for ( i = 0; i < 100 && line[i]; i++ ){
			switch( line[i] ){
				case 'A':   line[i] = 'T';  break;
				case 'a':   line[i] = 't';  break;
				case 'C':   line[i] = 'G';  break;
				case 'c':   line[i] = 'g';  break;
				case 'G':   line[i] = 'C';  break;
				case 'g':   line[i] = 'c';  break;
				case 'T':   line[i] = 'A';  break;
				case 't':   line[i] = 'a';  break;
				default:	break;
			}
		}
		outfile << line << endl;
	}
	infile.close();
	outfile.close();
	return 1;
}

