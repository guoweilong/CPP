/* 
 Guo, Weilong; guoweilong@gmail.com; 2013-07-01
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
};

parameter param;

class Region
{
public:
	unsigned int    left;
	unsigned int    right;
	char            strand;
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
	if (x.strand != y.strand)
		return x.strand < y.strand;
    return false;
}
// Region with same boundary but different strands would be considered as same region
bool operator==(const Region &x, const Region &y)
{
        return (x.left == y.left && x.right == y.right && x.strand == y.strand);
}


class Site {
public:
	unsigned int	pos;
    char            strand;
	unsigned int    MC;
	unsigned int    ALLC;
};


bool operator<(const Site &x, const Site &y)
{
	if (x.pos != y.pos)
		return x.pos < y.pos;
	if (x.strand != y.strand)
		return x.strand < y.strand;
    return false;
}

bool operator==(const Site &x, const Site &y)
{
    return (x.pos == y.pos && x.strand == y.strand);
}


void exit_with_help( void )
{
	printf(
		"Usage:	StrandSpecificMCinRegion -s <site.bed> -r <region.bed>\n"
		"Author: Guo, Weilong; guoweilong@gmail.com; 2013-07-01\n"
        "INPUT:\n"
		"-s	Input, Source file, store sites information\n"
		"	chr12  19898777    +   10  20 XXX\n"
		"-r	Input, Region file, BED file to store regions\n"
		"	chr12 19898766 19898966 + XXX\n"
        "OUTPUT:\n"
		"Ex:  chr12  19898777    +   10  20   19898766 19898966 +\n"
        "Des: (1)chr (2)site_pos (3)site_strand (4)mC (5)All_C (6) region_left\n"
        "     (7)region_right (8)region_strand\n"
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
}

// The first 4 columns should be ...
// chrom   start_pos   end_pos     strand
// chr12   19898766    19898966    +
//

int ReadRegionFile (map<string, vector<Region> > & AllRegions) {
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
		if ( (giter = AllRegions.find(chr)) == AllRegions.end() ) {
			vector<Region> trv;
			AllRegions[chr] = trv;
		}
		unsigned int left = atoi(tokens[1].c_str());
		unsigned int right = atoi(tokens[2].c_str());
		char strand = tokens[3][0];
		Region region(left, right, strand);
		AllRegions[chr].push_back(region);
	}
	Rfile.close();
	map<string, vector<Region> >::iterator giter;
	for ( giter = AllRegions.begin(); giter != AllRegions.end(); giter++ ) {
		string chr = giter->first;
		sort(AllRegions[chr].begin(), AllRegions[chr].end());
        unique(AllRegions[chr].begin(), AllRegions[chr].end());
	}
	return 1;
}

int OutputRegion (map<string, vector<Region> > & AllRegions) {
    map<string, vector<Region> >::iterator citer;
    for (citer = AllRegions.begin(); citer != AllRegions.end(); citer++) {
        string chr = citer->first;
        vector<Region>::iterator riter;
        for (riter = AllRegions[chr].begin(); riter != AllRegions[chr].end(); riter++) {
            cout << chr << "\t" << riter->left << "\t" << riter->right << "\t" << riter->strand << endl;
        }
    }
    return 1;
}



// Site file format:
// strand  pos strand  #mC  #all_C
// chr10   54718   -   8   10
// chr10   54720   -   0   10
//

int ReadSiteFile (map<string, vector<Site> > & AllSites) {
	ifstream Sfile(param.sitefile.c_str());
	if(!Sfile) {
		cout << "Can not open input file : " << param.sitefile.c_str() << endl;
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
		Site one_site;
		string chr = tokens[0];
		one_site.pos = atoi(tokens[1].c_str());
		one_site.strand = tokens[2][0];
		one_site.MC = atoi(tokens[3].c_str());
		one_site.ALLC = atoi(tokens[4].c_str());
        if ( AllSites.find(chr) == AllSites.end() )
            AllSites[chr] = vector<Site>(); // To check grammar right or not
        AllSites[chr].push_back(one_site);
	}
	
	Sfile.close();
    map<string, vector<Site> >::iterator siter;
	for ( siter = AllSites.begin(); siter != AllSites.end(); siter++ ) {
		string chr = siter->first;
		sort(AllSites[chr].begin(), AllSites[chr].end());
        unique(AllSites[chr].begin(), AllSites[chr].end());
	}
	return 1;

}


int OutputSite ( map<string, vector<Site> > & AllSites) {
    map<string, vector<Site> >::iterator citer;
    for (citer = AllSites.begin(); citer != AllSites.end(); citer++) {
        string chr = citer->first;
        vector<Site>::iterator siter;
        for (siter = AllSites[chr].begin(); siter != AllSites[chr].end(); siter++) {
            cout << chr << "\t" << siter->pos << "\t" << siter->strand << "\t" << siter->MC << "\t" << siter->ALLC << endl;
        }
    }
	return 1;
}

class MapInfo
{
public:
    string          chr;
    unsigned int    s_pos;
	char            s_strand;
    unsigned int    s_MC;
	unsigned int    s_ALLC;
	unsigned int    r_left;
	unsigned int    r_right;
	char            r_strand;
    
	MapInfo() {};
	MapInfo(string CHR, unsigned int S_POS, char S_STRAND, unsigned int S_MC, unsigned int S_ALLC, unsigned int R_LEFT, unsigned int R_RIGHT, char R_STRAND):chr(CHR),s_pos(S_POS),s_strand(S_STRAND),s_MC(S_MC),s_ALLC(S_ALLC),r_left(R_LEFT),r_right(R_RIGHT),r_strand(R_STRAND){};
};

int OutputMapInfo ( vector<MapInfo> & MappedSites ) {
    vector<MapInfo>::iterator miter;
    for (miter = MappedSites.begin(); miter != MappedSites.end(); miter++) {
        cout << miter->chr << "\t" << miter->s_pos << "\t" << miter->s_strand << "\t" << miter->s_MC << "\t"  << miter->s_ALLC << "\t" << miter->r_left << "\t" << miter->r_right << "\t" << miter->r_strand << endl;
    }
}

int GetMappedSites(string chr, vector<Site> & S, Region & R, vector<MapInfo> & MappedSites) {
    unsigned int size = S.size();
    //Binary search
    //  The 1st sites of which postion is no less than left position of Region
    int start = 0;
    int end = (int)size - 1;
    int mid;
    while (start <= end) {
        mid = (start + end) / 2;
        if (S[mid].pos > R.left){
            end = mid - 1;
        } else if (S[mid].pos < R.left) {
            start = mid + 1;
        } else {
            break;
        }
    }
    //Get the region
    int i = start;
    while (i < size && S[i].pos <= R.right ) {
        MappedSites.push_back( MapInfo(chr, S[i].pos, S[i].strand, S[i].MC, S[i].ALLC, R.left, R.right, R.strand) );
        i++;
    }
    
    return 1;
}

int StrandSpecificMCinRegion(map<string, vector<Site> > & AllSites, map<string, vector<Region> > & AllRegions){
    vector<MapInfo> MappedSites;
    map<string, vector<Region> >::iterator citer;
    for (citer = AllRegions.begin(); citer != AllRegions.end(); citer++) {
        string chr = citer->first;
        vector<Region>::iterator riter;
        for (riter = AllRegions[chr].begin(); riter != AllRegions[chr].end(); riter++) {
            GetMappedSites(chr, AllSites[chr], *riter, MappedSites);
        }
    }
    OutputMapInfo(MappedSites);
    return 1;
}


int main(int argc, char* argv[])
{
	init();

	parse_command_line(argc,argv);
    //cout << "Before Read Region file" <<endl;
	map< string, vector<Region> > AllRegions;
	ReadRegionFile(AllRegions);
    //OutputRegion(AllRegions);
    
	//cout << "Before Read Site file" << endl;
	map< string, vector<Site> > AllSites;
	ReadSiteFile(AllSites);
    //OutputSite(AllSites);

    StrandSpecificMCinRegion( AllSites, AllRegions);
    
	return 0;
}


