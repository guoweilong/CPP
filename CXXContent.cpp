/* CXXContent.cpp
* Author: GUO Weilong @ 2010-12-26
* Usage: CoverageSiteToRegion -s <site.bed> -r <region.bed>
* *********************
* Modification:
*/

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

struct parameter
{
	string fastafile;	      // -f
//	string outputfile;	// -o
	int    center;		// -c
      bool   singlestrand;    // -s
};

parameter param;


struct seqline{
	string name;
	string seq;
};

/** Read in a file which is in FASTA format, and store in
 *  a 'vector<seqline>' structure fasta.
 */
int ReadFasta( string filename, vector<seqline> & fasta )
{
	fasta.clear();
	ifstream fa(filename.c_str());
	if(!fa) {
		cerr << "Cannot open fasta file " << filename << "\n";
		return -1;
	}
	seqline sl;
	sl.name = "";
	sl.seq = "";
	while(!fa.eof()){
		char buffer[10000 + 1];
		fa.getline(buffer, 10000);
		// For UNIX system, the end for line might be '\n\r'
		if(buffer[strlen(buffer) - 1] == '\r'){
			buffer[strlen(buffer) - 1] = '\0';
		}
		string tmp = buffer;
		if(buffer[0] == '>'){
			if( sl.name.compare("") ) {
				fasta.push_back(sl);
			}
			sl.name = tmp.substr(1);
			sl.seq = "";
		}else{
			sl.seq += tmp;
		}
	}
	fasta.push_back(sl);
	fa.close();
	return fasta.size();
}

void exit_with_help( void )
{
	printf(
		"Usage: CXXContent -c <int> -f <input.fa> [-s]\n"
		"Author: GUO Weilong	Version:2010/12/29\n"
		"options:\n"
		"-f	Input FASTA file\n"
		"-c	where the center is in the FASTA file\n"
		"	[default: 3] The third position in FASTA file\n"
            "-s   Only count on single strand (+ strand) when specified\n"
            "     [default: both strands will be counted]\n"
		"Function: If the region of the input FASTA file is [-A,+B]\n"
		"	The output are counts of CG,CHG,CHH in the region of [-(A-2),+(B-2)]\n"
		//"-o	Output file\n"
		"Output format:\n"
		"	[not show] [position	CG	CHG	CHH]\n"
		"			0	15	5	3\n"
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
		case 'f':
			param.fastafile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		/*case 'o':
			param.outputfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;*/
		case 'c':
			param.center = atoi(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
            case 's':
                  param.singlestrand = true;
                  break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.fastafile.length() ){
		exit_with_help();
	}
}

int GetContent (vector<seqline> &Seq) {
	vector<seqline>::iterator iter;
	iter = Seq.begin();
	size_t len = iter->seq.length();
	if (len < 5) {
		cerr << "The lengths of the sequences are too short, should be >=5\n";
		exit(1);
	}
	unsigned int *CG = new unsigned int[len - 4]; 
	unsigned int *CHG = new unsigned int[len - 4];
	unsigned int *CHH = new unsigned int[len - 4]; 
	// Cut two basepair on each end
	for (int i = 0; i < len - 4; i++) {
		CG[i] = CHG[i] = CHH[i] = 0;
	}
	for (iter = Seq.begin(); iter != Seq.end(); iter++) {
		string read = iter->seq;
		//if (read.length() != len) {
		if (read.length() < len) {
			cerr << "The lengths of sequences are not the same!" << endl;
			exit(1);
		}
		ToUpperString (read);
		// From the + strand
		for (int i = 0; i < len - 4; i++) {
			if (read[i+2] == 'C') {
				if (read[i+3] == 'G')	CG[i]++;
				else if (read[i+4] == 'G') CHG[i]++;
				else	CHH[i]++;
			}
		}
            if (!param.singlestrand) {
		// From the - strand
		for (int j = len - 5; j >= 0; j--) {
			if (read[j+2] == 'G') {
				if (read[j+1] == 'C')	CG[j]++;
				else if (read[j] == 'C') CHG[j]++;
				else	CHH[j]++;
			}
		}
            }
	}
	for (int i = 0; i < len - 4; i++) {
		cout << i - param.center << "\t" << CG[i] << "\t" << CHG[i] << "\t" << CHH[i] << endl;
	}
	return 0;
}

void init(){
	// default values
	param.fastafile = "";
	//param.outputfile = "";
	param.center = 0;
      param.singlestrand = false;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.fastafile = string("test.fa");
	param.center = 102;
#else
	parse_command_line(argc,argv);
#endif
	vector<seqline> Seq; 
	ReadFasta (param.fastafile, Seq);
	GetContent (Seq);
	return 0;
}


