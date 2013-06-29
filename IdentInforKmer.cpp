/* Usage:IdentInforKmer <foreground.fa> <background.fa> 4 
 * Weilong GUO 2010-05-21
 */
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include "math.h"

double CountKmer(char* filename, int *Count, int MAX, int K) {
	// Initialization for counting result of each pattern
	fstream infile;
	infile.open(filename);
	double total = 0;
	for(int i = 0; i < MAX; i++ ) 
		Count[i]=0;
	char ch;
	int pos = 0;
	int v = 0;
	char temp[100];
	while( infile.get(ch) ){
		switch ( ch ) {
		case '>':
			infile.getline( temp, 100 );
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
			Count[v]++;
			total++;
		}
	}
	infile.close();
	return total;
}

int main(int argc, char* argv[])
{
#ifdef _DEBUG
	int K = 4;
#else
	if( argc < 4 ) {
		cout << "CountKmers foreground.fa background.fa 4\n";
		return 0;
	} 
	int K = atoi( argv[3] ); //length of reads to be counted
#endif

	long MAX = (long) pow( 4.0, K );
	int *ForeCount, *BackCount;
	ForeCount = new int[MAX];
	BackCount = new int[MAX];
	double ForeTotal, BackTotal;

#ifdef _DEBUG
	ForeTotal = CountKmer( "foreground.fa", ForeCount, MAX, K );
	BackTotal = CountKmer( "background.fa", BackCount, MAX, K );
#else
	ForeTotal = CountKmer( argv[1], ForeCount, MAX, K );
	BackTotal = CountKmer( argv[2], BackCount, MAX, K );
#endif

	double II = 0;
	double SS = 0;
	double f, b;
	for(int i = 0; i < MAX; i++ ){
		if( ForeCount[i] && BackCount[i] ) {
			f = ForeCount[i] / ForeTotal;
			b = BackCount[i] / BackTotal;
			II += f * log( f / b ) / log(2.0);
			SS += b * log( b / f ) / log(2.0);
		}
	}
	cout << II << endl;
	cout << SS << endl;
	return 1;
}
