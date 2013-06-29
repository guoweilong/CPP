// GetCoveredSiteNo.cpp

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;
//#include "math.h"

struct parameter
{
	string sitefile;	// -s
	string regionfile;	// -r
	bool freq;	// -f
};

parameter param;

class Region
{
public:
	unsigned int left;
	unsigned int right;
	char strand;
	string other;
//	unsigned int count; // used before the version 2011-03-29
	float count;
	Region(){ left = 0; right = 0; strand = '+'; other = ""; count  = 0; }
	Region(unsigned int l, unsigned int r, char s):left(l),right(r),strand(s){ other = ""; count = 0; }
	Region(unsigned int l, unsigned int r, char s, string o):left(l),right(r),strand(s),other(o){ count = 0; }
};

// For Struct Region 
bool operator<(const Region &x, const Region &y)
{
	if (x.left != y.left) 
		return x.left < y.left;
	if (x.right != y.right)
		return x.right < y.right;
      // The function should return under all the conditions
      return false;
}

// For the same site in different strand, it's considered to be the same
bool operator==(const Region &x, const Region &y)
{
        return (x.left == y.left && x.right == y.right);
}

class Site
{
public:
	unsigned int pos;
	char strand;
	unsigned short mC;
	unsigned short tC;
	Site(){ pos= 0; strand = '+'; tC = 1; mC = 1;}
	Site(unsigned int p, char s):pos(p),strand(s){ tC = 1; mC = 1;}
	Site(unsigned int p, char s, unsigned short m, unsigned short t):pos(p), strand(s), mC(m), tC(t){}
};

// For struct site
bool operator<(const Site &x, const Site &y){
	if (x.pos != y.pos)
		return x.pos < y.pos;
	if (x.strand != y.strand)
		return x.strand < y.strand;
	return false;
}

bool operator==(const Site &x, const Site &y){
	return (x.pos == y.pos && x.strand == y.strand);
}

void exit_with_help( void )
{
	printf(
		"Author: GUO Weilong, 2011-04-08, guoweilong@gmail.com\n"
		"Usage: GetCoveredSitesNo -s <site.txt> -r <region.txt> [-f] \n"
		"options:\n"
		"-s Input, Source file, store sites information\n"
		"     12  19898777 +  XXXX XXXX \n"
		"-r Input, Region file, BED file to store regions\n"
		"     chr12 19898766 19898966 + 0.35896\n"
		"-f, a symble switch, when specified, output are the sum of frequencies of\n"
		"     methylated sites; or else, [DEFAULT] output the number of these sites.\n"
		"Output, Distribution file, output the number of the sites lie in the regions\n"
		"     chr12 19898766 1989966 + 0.35896 121\n"
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
			param.sitefile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'r':
			param.regionfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'f':
			param.freq = true;
			break; 
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.sitefile.length() || !param.regionfile.length() ){
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
	param.regionfile = "";
	param.freq = false;
}

int ReadRegionFile (map<string, vector<Region> > & GenomeRegion) {
	ifstream Rfile(param.regionfile.c_str());
	if(!Rfile) {
		cout << "cannot open input file" << param.regionfile.c_str() << endl;
		return -1;
	}
	int k = 0;
	while (!Rfile.eof()) {
		char buffer[1000+1];
		Rfile.getline(buffer, 1000);
		if(!strlen(buffer))
			continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
		string chr = tokens[0];
		map< string,vector<Region> >::iterator giter;
		if ( (giter = GenomeRegion.find(chr)) == GenomeRegion.end() ) {
			vector<Region> trv;
			GenomeRegion[chr] = trv;
		}
		unsigned int left = atoi(tokens[1].c_str());
		unsigned int right = atoi(tokens[2].c_str());
		char strand = tokens[3][0];

		unsigned int n = tokens.size();
		string other = "";
		
		for ( int i = 4; i < n; i++ ) {
			other += "\t"; other += tokens[i];
		}

		Region region(left, right, strand, other);
		GenomeRegion[chr].push_back(region);
	}
	map<string, vector<Region> >::iterator giter;
	for ( giter = GenomeRegion.begin(); giter != GenomeRegion.end(); giter++ ) {
		string chr = giter->first;
		sort(GenomeRegion[chr].begin(), GenomeRegion[chr].end());
        	unique(GenomeRegion[chr].begin(), GenomeRegion[chr].end());
	}
	Rfile.close();
	return 0;
}

int ReadSiteFile (map<string, vector<Site> > & GenomeSite) {
	ifstream Sfile(param.sitefile.c_str());
	// Sitefile format: 1       5636    -       CHG     3       9
	if(!Sfile) {
		cout << "cannot open input file" << param.sitefile.c_str() << endl;
		return -1;
	}
//	cout << "ReadSiteFile{}\n";
	while (!Sfile.eof()) {
		char buffer[1000+1];
		Sfile.getline(buffer, 1000);
		if(!strlen(buffer))
			continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
//		cout << "tokens.size()" << tokens.size() << endl;
		if ( tokens.size() != 6) {
			cerr << "Input file have wrong number of columns\n" ;
			exit(1);
		}
		
		string chr = string("chr") + tokens[0];
		map< string,vector<Site> >::iterator giter;
		
//		cout << "Mark 1\n" ;
		if ( (giter = GenomeSite.find(chr)) == GenomeSite.end() ) {
			vector<Site> vecsite; vecsite.clear();
			GenomeSite[chr]	= vecsite;
		}
		
//		cout << "Mark 2\n";
		unsigned int pos = atoi(tokens[1].c_str());
		char strand = tokens[2][0];
		// Need to check for type casting
		unsigned short mC = (unsigned short) atoi(tokens[4].c_str());
		unsigned short tC = (unsigned short) atoi(tokens[5].c_str());
//		cout << "Mark 3\n" << mC << "\t" << tC << endl;
//		cout << "Site\t" << pos << "\tstrand\t" << strand << "\tmc\t" << mC << "\trC\t" << tC << endl ;
		GenomeSite[chr].push_back( Site(pos, strand, mC, tC) );
	}
//	cout << "Before sort\n";
	map< string, vector<Site> >::iterator miter;
	for ( miter = GenomeSite.begin(); miter != GenomeSite.end(); miter++ ) {
		string str = miter->first;
		sort(GenomeSite[str].begin(), GenomeSite[str].end());
		unique(GenomeSite[str].begin(), GenomeSite[str].end());
	}
	Sfile.close();
	return 0;
}

int CoverageAnalysis( map<string, vector<Site> > & GenomeSite, map<string, vector<Region> > & GenomeRegion)
{
	map<string, vector<Region> >::iterator GRiter;
	for ( GRiter = GenomeRegion.begin(); GRiter != GenomeRegion.end(); GRiter++ ) {
		string chr = GRiter->first;
		vector<Region>::iterator Riter;
		
		vector<Site>::iterator Siter, preSiter;
		Siter = GenomeSite[chr].begin();
		preSiter = Siter;

		for ( Riter = GenomeRegion[chr].begin(); Riter != GenomeRegion[chr].end() && Siter != GenomeSite[chr].end(); Riter++ ) {
		//	cout << "Riter" << Riter->left << "\t" << Riter->right << endl;
			if ( Riter->left < Siter->pos )
				Siter = preSiter;
			while ( Siter->pos < Riter->left ) 
				Siter++;
			preSiter = Siter;
			unsigned int count = 0;
			float freq = 0;
		//	cout << "Riter" << Riter->left << "\t" << Riter->right << "\t" << Siter->pos << endl;
			while ( Siter != GenomeSite[chr].end() && Siter->pos <= Riter->right ) {
		//		cout << Siter->pos << endl;
				float mC = Siter->mC;
				float tC = Siter->tC;
//				cout << "pos\t" << Siter->pos << "\tmc\t" << mC << "\ttC\t" << tC << endl;
				freq += ( mC / tC );
				count++;	Siter++;
			}
			if ( !param.freq ) {
				Riter->count = count;
//				cout << "Using count\n";
			} else {
//				cout << "freq\t" << freq << endl;
				Riter->count = freq;
			}
		//	cout << "hello" << endl;
		}
		if ( Siter == GenomeSite[chr].end() ) {
			while ( Riter != GenomeRegion[chr].end() ) {
				Riter->count = 0; Riter++;
			}
		}

	}
	return 1;
}

int OutputRegions ( map<string, vector<Region> > & GenomeRegion ) {
	map<string, vector<Region> >::iterator GRiter;
	for ( GRiter = GenomeRegion.begin(); GRiter != GenomeRegion.end(); GRiter++ ) {
		string chr = GRiter->first;
		vector<Region>::iterator Riter;
		for ( Riter = GenomeRegion[chr].begin(); Riter != GenomeRegion[chr].end(); Riter++ ) {
			cout << chr << "\t" << Riter->left << "\t" << Riter->right << "\t"
			<< Riter->strand /*<< "\t"*/ << Riter->other << "\t" << Riter->count << endl;
		}
	}
	return 1;
}

int OutputSites ( map<string, vector<Site> > & GenomeSite) {
	map<string, vector<Site> >::iterator GSiter;
	for ( GSiter = GenomeSite.begin(); GSiter != GenomeSite.end(); GSiter++ ) {
		string chr = GSiter->first;
		vector<Site>::iterator Siter;
		for ( Siter = GenomeSite[chr].begin(); Siter != GenomeSite[chr].end(); Siter++ ) {
			cout << chr << "\t" << Siter->pos << "\t" << Siter->strand << "\t" 
			<< Siter->mC << "\t" << Siter->tC << endl;
			
		}
	}
	
	return 1;
}

int main(int argc, char* argv[])
{
	init();
//	cout << "Welcome\n";
#ifdef _DEBUG
	param.sitefile = string("site.bed");
	param.regionfile = string("region.bed");
#else
	parse_command_line(argc,argv);
#endif
//	cout << "Parse_command_line\n";
	map< string, vector<Site> > GenomeSite;
	map< string, vector<Region> > GenomeRegion; 
	//store all regions in each chromosome
	ReadSiteFile ( GenomeSite );
//	cout << "ReadSiteFile\n";
//	OutputSites ( GenomeSite );

	ReadRegionFile ( GenomeRegion );
//	cout << "ReadRegionFile\n";

//	cout << "No:\t" << GenomeSite["chr1"].size() << "\t" << GenomeRegion["chr1"].size() << endl;

	CoverageAnalysis( GenomeSite, GenomeRegion );
//	cout << "CoverageAnalysis\n";
	OutputRegions ( GenomeRegion );
	return 0;
}


