/* CoverageSiteToRegion
 * Fuction: Find all the sites in site.bed file are or not covered by any region in the region.bed file
 * Author: GUO Weilong @ 2010-12-06
 * *********************
 * Modification:
 *	2010-12-20: Correct the way to show the positions of methylated sites, 
 *		make the '+' and '-' strands different
 *	2010-12-23: Find the sites by binary search
 *	2010-12-29: Add option for methods to search;
 *	2011-05-09: Add option '-d' to output relative location information for [D]ebug;
 *  2011-05-23: Add option '-m' to output the sequencing information of [M]ethyl-C sites;
 */

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
//	bool bruteforce;    // -b
	bool debug;			// -d
	bool methyl;		// -m
};

parameter param;

class Region
{
public:
	unsigned int left;
	unsigned int right;
	char strand;
	Region () { left = 0; right = 0; strand = '+'; }
	Region (unsigned int l, unsigned int r, char s):left(l),right(r),strand(s){
	}
};

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

void exit_with_help( void )
{
	printf(
		"Usage:	CoverageSiteToRegion -s <site.bed> -r <region.bed> [-d -m]\n"
		"Author:	Guo, Weilong	guoweilong@gmail.com 	2011-05-09\n"
		"Options:\n"
		"-s	Input, Source file, store sites information\n"
		"	12  19898777 +  XXXX XXXX \n"
		"-r	Input, Region file, BED file to store regions\n"
		"	chr12 19898766 19898966 + XXXXXXX XXXX\n"
//		"-b	(optional switch) Bruteforce searching if specified, or else\n"
//		"	binary searching is used by default\n"
		"-m	(optional switch) Output sequence information of mC\n"
		"-d	(optional switch) Output coverage information when specified\n"
		"	Output file, distribution file, output the distribution of the\n"
		"	sites lie in the regions. Format:\n"
		"	0\t34\n\t1\t12\n\t2\t20\n"
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
			param.regionfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
//		case 'b':
//			param.bruteforce = true;
//			break;
		case 'd':
			param.debug = true;
			break;
		case 'm':
			param.methyl = true;
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
//	param.bruteforce = false;
	param.debug = false;
	param.methyl = false;
}

int ReadRegionFile (map<string, vector<Region> > & Genome) {
	ifstream Rfile(param.regionfile.c_str());
	if(!Rfile) {
		cout << "cannot open input file \"" << param.regionfile.c_str() << "\"\n";
		exit(1);
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
		if ( (giter = Genome.find(chr)) == Genome.end() ) {
			vector<Region> trv;
			Genome[chr] = trv;
		}
		unsigned int left = atoi(tokens[1].c_str());
		unsigned int right = atoi(tokens[2].c_str());
		char strand = tokens[3][0];
		Region region(left, right, strand);
		Genome[chr].push_back(region);
	}
	map<string, vector<Region> >::iterator giter;
	for ( giter = Genome.begin(); giter != Genome.end(); giter++ ) {
		string chr = giter->first;
		sort(Genome[chr].begin(), Genome[chr].end());
            unique(Genome[chr].begin(), Genome[chr].end());
	}
	Rfile.close();
	return 0;
}

int OutputSite (int rpos, string chr, int pos, unsigned int left, unsigned int right, char strand, unsigned int MC, unsigned int ALLC) {
	// rpos = relative position
	if ( param.methyl ) {
		if ( param.debug ) {
			cout << rpos << "\t" << chr << "\t" << pos << "\t[" << left << ",\t" << right << "]\t" << strand << "\t" << MC << "\t" << ALLC << endl; 
		} else {
			cout << rpos << "\t" << MC << "\t" << ALLC << endl; 
		}
	} else {
		if ( param.debug ) {
			cout << rpos << "\t" << chr << "\t" << pos << "\t[" << left << ",\t" << right << "]\t" << strand << endl; 
		} else {
			cout << rpos << endl; 
		}
	}
	return 1;
}

int MapSiteFile (map<string, vector<Region> > & Genome) {
	ifstream Sfile(param.sitefile.c_str());
	// Sitefile format: 1       5636    -       CHG     3       9
	if(!Sfile) {
		cout << "cannot open input file" << param.sitefile.c_str() << endl;
		return -1;
	}
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
		string chr = string("chr") + tokens[0];
		map< string,vector<Region> >::iterator giter;
		if ( (giter = Genome.find(chr)) == Genome.end() )
			continue;

		unsigned int pos = atoi(tokens[1].c_str());
		size_t size = Genome[chr].size();
		if (size == 0)    continue;
		unsigned int MC = atoi(tokens[4].c_str());
		unsigned int ALLC = atoi(tokens[5].c_str());

		/*Binary search*/
		int start = 0;
		int end = (int)size - 1;
		int mid;
		while (start <= end) {
		      mid = (start + end) / 2;
			if (pos < Genome[chr][mid].left){
				end = mid - 1;
		#ifdef _DEBUG
				cout << Genome[chr][mid].left << endl;
		#endif
			      } else if (pos > Genome[chr][mid].right) {
				      start = mid + 1;
		#ifdef _DEBUG
				cout << Genome[chr][mid].right << endl;
		#endif
			} else {
				break;
			}
		}
		unsigned int left = Genome[chr][mid].left;
		unsigned int right = Genome[chr][mid].right;
		char strand = Genome[chr][mid].strand;

		if (pos >= left && pos <= right) {
			if (strand == '+') {
				OutputSite (pos - left, chr, pos, left, right, strand, MC, ALLC) ;
			} else {
				OutputSite (right - pos, chr, pos, left, right, strand, MC, ALLC) ;
			}
		}
		int p = mid-1;
		while (p >= 0 && pos >= Genome[chr][p].left && pos <= Genome[chr][p].right) {
			if (Genome[chr][p].strand == '+') {
				OutputSite (pos - Genome[chr][p].left, chr, pos, Genome[chr][p].left, Genome[chr][p].right, Genome[chr][p].strand, MC, ALLC) ;
			} else {
				OutputSite (Genome[chr][p].right - pos, chr, pos, Genome[chr][p].left, Genome[chr][p].right, Genome[chr][p].strand, MC, ALLC) ;
			}
			p--;
		}
		p = mid+1;
		while (p < size && pos >= Genome[chr][p].left && pos <= Genome[chr][p].right) {
			if (Genome[chr][p].strand == '+') {
				OutputSite (pos - Genome[chr][p].left, chr, pos, Genome[chr][p].left, Genome[chr][p].right, Genome[chr][p].strand, MC, ALLC) ;
			} else {
				OutputSite (Genome[chr][p].right - pos, chr, pos, Genome[chr][p].left, Genome[chr][p].right, Genome[chr][p].strand, MC, ALLC) ;
			}
			p++;
		}
	}
	
	Sfile.close();
	return 0;
}

int main(int argc, char* argv[])
{
	init();

#ifdef _DEBUG
	param.sitefile = string("site.bed");
	param.regionfile = string("region.bed");
#else
	parse_command_line(argc,argv);
#endif

	map< string, vector<Region> > Genome; 
	//store all regions in each chromosome

	//cout << "Before Read Region file" <<endl;
	ReadRegionFile(Genome);
	//cout << "After Read Region file" << endl;
	MapSiteFile(Genome);
	return 1;
}


