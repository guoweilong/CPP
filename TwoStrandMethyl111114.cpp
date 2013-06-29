/* TwoStrandMethyl.cpp
 * Modification:
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cmath>

using namespace std;


struct parameter
{
	string sitefile;	// -s
	string refFlatfile;	// -r
	string IR;		// -I
};

parameter param;

class Site {
public:
	unsigned int	pos;
	char	strand;
	float	ML; // Methylation level
};


bool operator<(const Site &x, const Site &y)
{
	return x.pos < y.pos;
}

// For the same site in different strand, it's considered to be the same
bool operator==(const Site &x, const Site &y)
{
        return (x.pos == y.pos);
}

void exit_with_help( void )
{
	printf(
		"Usage:	TwoStrandMethyl -s <site.bed> -r <refFlat> \n" //[-I <region_of_Intron>]
		"Author:	Guo, Weilong	guoweilong@gmail.com 	2011-11-13\n"
		"Options:\n"
		"-s	Input, Source file, store sites information\n"
		"	12  19898777 +  XXXX XXXX \n"
		"-r	Input, refFlat file \n"

		);
	exit(1);
}

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	int i;

	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
		case 's':
			if(i == argc)	exit_with_help();
			param.sitefile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'r':
			if(i == argc)	exit_with_help();	
			param.refFlatfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.sitefile.length() || !param.refFlatfile.length() ){
		exit_with_help();
	}
}


vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
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

void init(){
	// default values
	param.sitefile = "";
	param.refFlatfile = "";
}


#define MAX 10000

int ReadSitesFile ( map< string, vector<Site> > & Genome ) {
	ifstream Sfile(param.sitefile.c_str());
	// Sitefile format: 1       5636    -       CHG     3       9
	if(!Sfile) {
		cout << "cannot open input file" << param.sitefile.c_str() << endl;
		return -1;
	}
	while ( !Sfile.eof() ) {
		char buffer[MAX+1];
		Sfile.getline(buffer, MAX);
		if(!strlen(buffer)) continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
		Site site;
		string chr = string("chr") + tokens[0];
		site.pos = atoi ( tokens[1].c_str() );
		site.strand = tokens[2][0];
		site.ML = atof ( tokens[4].c_str() ) / atof (tokens[5].c_str() ) ;
		if ( Genome.find(chr) == Genome.end() ) {
			vector<Site> tsv;
			Genome[chr] = tsv; // ?
		}
		Genome[chr].push_back(site);
	}
	Sfile.close();
	// Sort
        map<string, vector<Site> >::iterator giter;
        for ( giter = Genome.begin(); giter != Genome.end(); giter++ ) {
                string chr = giter->first;
                sort(Genome[chr].begin(), Genome[chr].end());
		unique(Genome[chr].begin(), Genome[chr].end());
        }

	return 1;
}

// To test function { ReadSitesFile }
int WriteSitesFile ( map< string, vector<Site> > & Genome ) {
        map< string, vector<Site> >::iterator giter;
        for ( giter = Genome.begin(); giter != Genome.end(); giter++ ) {
		string chr = giter->first;
		vector<Site>::iterator siter;
		for ( siter = Genome[chr].begin(); siter != Genome[chr].end(); siter++ ) {
			cout << chr << "\t" << siter->pos << "\t" << siter->strand << "\t"
			     << siter->ML << endl;
		}
	}
	return 1;
}

class RefFlat {
public:
	string 	gene;
	string 	trans;
	string 	chr;
	char	strand;
	unsigned int tStart;	// Start of transcript
	unsigned int tEnd;	// End of transcript
	unsigned int cStart;	// Start of coding region
	unsigned int cEnd;	// End of coding region
	short	nExon;		// Number of exons
	vector<unsigned int> eStart;	// Starts of exons
	vector<unsigned int> eEnd;	// Ends of exons
};

/*
 format of refFlat file
 DLX1	NM_178120	chr2	+	172658453	172662647	172658651	172661231	3	172658453,172659627,172660976,	172658964,172659827,172662647,
*/

int ReadRefFlatFile ( list<RefFlat> & Annotation ) {
	ifstream Rfile(param.refFlatfile.c_str());
	if(!Rfile) {
		cout << "cannot open input file" << param.refFlatfile.c_str() << endl;
		return -1;
	}
	while ( !Rfile.eof() ) {
		char buffer[MAX+1];
		Rfile.getline(buffer, MAX);
		if(!strlen(buffer)) continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);

		RefFlat rf;
		rf.gene = tokens[0];
		rf.trans = tokens[1];
		rf.chr = tokens[2];
		rf.strand = tokens[3][0];
		rf.tStart = atoi ( tokens[4].c_str() );
		rf.tEnd = atoi ( tokens[5].c_str() );
		rf.cStart = atoi ( tokens[6].c_str() );
		rf.cEnd = atoi ( tokens[7].c_str() );
		rf.nExon = atoi ( tokens[8].c_str() );
		vector<string> st = string_tokenize( tokens[9], "," );
		int n = st.size();
		for (int i = 0; i < n; i++ ) {
			rf.eStart.push_back ( atoi( st[i].c_str() ) );
		}
		vector<string> et = string_tokenize( tokens[10], "," );
		n = et.size();
		for (int i = 0; i < n; i++ ) {
			rf.eEnd.push_back ( atoi( et[i].c_str() ) );
		}
		Annotation.push_back(rf);
	}
	Rfile.close();
	
	return 1;
}


// To test function { ReadRefFlatFile }
int WriteRefFlatFile ( list<RefFlat> & Annotation ) {
	list<RefFlat>::iterator liter;
	for ( liter = Annotation.begin(); liter != Annotation.end(); liter++ ) {
		cout << liter->gene << "\t" << liter->trans << "\t" << liter->chr << "\t"
		     << liter->strand << "\t" << liter->tStart << "\t" << liter->tEnd << "\t"
		     << liter->cStart << "\t" << liter->cEnd << "\t"  << liter->nExon << "\t";
		int n = liter->eStart.size();
		for ( int i = 0; i < n; i++ ) {
			cout << liter->eStart[i] << ",";
		}
		cout << "\t"; n = liter->eEnd.size();
		for (int i = 0; i < n; i++ ) {
			cout << liter->eEnd[i] << ",";
		}
		cout << endl;
	}
	return 1;
}

// Find all the methylated sites in the region [IntronLeft, IntronRight],
// and put them in the lists of WStrand (Watson) and CStrand (Click)
int RegionMethylation ( vector<Site> & Chr, unsigned int IntronLeft, unsigned int IntronRight,
	list<Site> & WStrand, list<Site> & CStrand) {
	int low, high, mid;
	low = 0; 
	if (Chr.size() == 0) return 0;
	high = Chr.size() - 1;
	mid = ( low + high ) / 2;
	while ( low <= high ) {
		mid = ( low + high ) / 2;
		if ( IntronLeft  > Chr[mid].pos ) { 
		// Not >=, so that can get low as the lowest one which >= IntronLeft
			low = mid + 1;
		} else {
			high = mid - 1;
		}
	}
	for ( unsigned int i = low; i<Chr.size() && Chr[i].pos <= IntronRight; i++ ) {
		if ( Chr[i].strand == '+' ) 
			WStrand.push_back( Chr[i] );
		else
			CStrand.push_back( Chr[i] );
	}
	return 1;
}

float MethylMean ( list<Site> & mC ) {
	float sum = 0; unsigned int count = 0;
	for ( list<Site>::iterator iter = mC.begin(); iter !=mC.end(); iter++, count++ ) 
		sum += iter->ML;
	return sum / count; 
}

#define SMALLNUM 0.0001

int TwoStrandMethylation ( map< string, vector<Site> > & Genome, list<RefFlat> & Annotation ) {
	list<RefFlat>::iterator liter;
//	list<RefFlat>::iterator NewAnno;
	for ( liter = Annotation.begin(); liter != Annotation.end(); liter++ ) {
		string chr = liter->chr;
		char strand = liter->strand;
		unsigned int n = liter->nExon;
		
		cout << liter->trans << "\t";
		
		// Exon
		list<Site> EWatson, ECrick;
		for ( int i = 0; i < n; i++ ) {
			RegionMethylation ( Genome[chr], liter->eStart[i], liter->eEnd[i], EWatson, ECrick);
		}
		cout << EWatson.size() << "\t" << ECrick.size() << "\t";
		
		// Intron
		list<Site> LWatson, LCrick; // left SS
		list<Site> MWatson, MCrick; // middle intron
		list<Site> RWatson, RCrick; // right SS

		for ( int i = 0; i < n-1; i++ ) {
			if( liter->eStart[i+1] - liter->eEnd[i] >= 500 ) {
				// Left SS
				RegionMethylation ( Genome[chr], liter->eEnd[i] + 5, liter->eEnd[i] + 105, LWatson, LCrick);

				// MI
				unsigned int mpos = ( liter->eEnd[i] + liter->eStart[i+1] ) / 2;
				RegionMethylation ( Genome[chr], mpos - 50, mpos + 50, MWatson, MCrick);

				// Right SS
				RegionMethylation ( Genome[chr], liter->eStart[i+1] - 105, liter->eStart[i+1] - 5, RWatson, RCrick);
			}
		}
		cout << LWatson.size() << "\t" << LCrick.size() << "\t";
		cout << MWatson.size() << "\t" << MCrick.size() << "\t";
		cout << RWatson.size() << "\t" << RCrick.size() << endl;
	}	
	
	return 1;
}

/*
		for ( int i = 1; i < n; i++ ) {
			unsigned int IntronLeft, IntronRight;
			IntronLeft = liter->eEnd[i-1]; IntronRight = liter->eStart[i];
			if ( param.IR == "whole" ) {
				if (IntronRight -10 > IntronLeft)
					RegionBiasedMethylation ( Genome[chr], IntronLeft+5, IntronRight-5, WStrand, CStrand);
			} else { 
				cout << "[ERROR] Unknown parameter:\t" << param.IR << endl << "============" << endl;
				exit_with_help(); 
			}
		}
		
		if ( WStrand.size() != 0 && CStrand.size() != 0 ) {
			float AntiToSense;
			float methylW = MethylMean ( WStrand ); // Walson strand
			float methylC = MethylMean ( CStrand ); // Click strand
			if (methylC >= SMALLNUM && methylW >= SMALLNUM ) {
				if ( strand == '+' ) {
					AntiToSense = methylC / methylW;
				} else {
					AntiToSense = methylW / methylC;
				}
			}
		} 
*/


int main(int argc, char* argv[])
{
	init ();
	parse_command_line ( argc, argv );

	map< string, vector<Site> > Genome;
	ReadSitesFile ( Genome );
//	WriteSitesFile ( Genome );

	list<RefFlat> Annotation;
	ReadRefFlatFile ( Annotation );
//	WriteRefFlatFile ( Annotation );

	TwoStrandMethylation ( Genome, Annotation );

	return 1;
}



