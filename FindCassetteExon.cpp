// FindCassetteExon.cpp

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
	string gtfFile;	// -s
};

parameter param;

//  Standard format of "gtf" format
//  <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

class Exon
{
public:
	unsigned int 	left;
	unsigned int 	right;
	unsigned short 	exon_number;
	Exon () { left = right = 0; exon_number = 0; }
	Exon ( unsigned int LEFT, unsigned int RIGHT, unsigned short EXON_NO):left(LEFT),right(RIGHT),exon_number(EXON_NO) {}
};

class Transcript{
public:
	unsigned int 	left;
	unsigned int 	right;
	double		fpkm;
	double		frac;
	vector<Exon>	exons;
	Transcript () {
		left = right = 0; fpkm = frac = 0;
		exons.clear();
	}
	Transcript ( unsigned int LEFT, unsigned int RIGHT, double FPKM, double FRAC):left(LEFT),right(RIGHT),fpkm(FPKM),frac(FRAC) {
		exons.clear();
	}
};

class Gene{
public:
	char			strand;
	map<string, Transcript>	isoforms;
	Gene () {};
	Gene (char STRAND) : strand (STRAND) { isoforms.clear(); }	
};

class Chromosome{
public:
	map<string, Gene>	genes;
};

// For Struct Exon

bool operator<(const Exon &x, const Exon &y)
{
	if (x.left != y.left) 
		return x.left < y.left;
	if (x.right != y.right)
		return x.right < y.right;
      // The function should return under all the conditions
      return false;
}

bool operator>(const Exon &x, const Exon &y)
{
	if (x.left != y.left)
		return x.left > y.left;
	if (x.right != y.right)
		return x.right > y.right;
	return false;
}

// For the same site in different strand, it's considered to be the same
bool operator==(const Exon &x, const Exon &y)
{
	return (x.left == y.left && x.right == y.right);
}

bool operator!=(const Exon &x, const Exon &y)
{
	return (x.left != y.left || x.right != y.right);		
}

void exit_with_help( void )
{
	printf(
		"Author: GUO Weilong, 2011-04-18, guoweilong@gmail.com\n"
		"Usage: FindCassetteExon -f <input.gtf> \n"
		"Options:\n"
		"-f Input, file in GTF format\n"
		"     chr1	Cufflinks	transcript	665900	666003	1000	-	.	gene_id \"CUFF.491\"; transcript_id \"CUFF.491.1\"; FPKM \"75.0190504609\"; frac \"1.000000\"; conf_lo \"57.696343\"; conf_hi \"92.341758\"; cov \"6.661290\";\n"
		"     chr1	Cufflinks	exon	665900	666003	1000	-	.	gene_id \"CUFF.491\"; transcript_id \"CUFF.491.1\"; exon_number \"1\"; FPKM \"75.0190504609\"; frac \"1.000000\"; conf_lo \"57.696343\"; conf_hi \"92.341758\"; cov \"6.661290\";\n"
		"History:	GTFParse 11-04-13\n\t\tFindCasseteExon 11-04-18\n"
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
			param.gtfFile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.gtfFile.length() ){
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
	param.gtfFile = "";
}

struct Transcript_Struct {
	string chr;
	string source;
	unsigned int left;
	unsigned int right;
	unsigned int score;
	char strand;
	char frame;
	string gene_id;
	string transcript_id;
	double fpkm;
	double frac;
	double conf_lo;
	double conf_hi;
	double cov;
};

int TransParser ( Transcript_Struct & trans, vector<string> & tokens ) {
	trans.chr = tokens[0];
	trans.source = tokens[1];
	trans.left = atoi(tokens[3].c_str());
	trans.right = atoi(tokens[4].c_str());
	trans.score = atoi(tokens[5].c_str());
	trans.strand = tokens[6][0];
	trans.frame = tokens[7][0];
	
	string tmp = tokens[8];
	
	vector<string> comments = string_tokenize(tmp, ";");
	vector<string>::iterator siter;
	for ( siter = comments.begin(); siter != comments.end(); siter++) {
		string title = ((*siter)[0] == ' ') ? siter->substr(1) : (*siter) ;  // In case some are begin with " " [gap]
		title = title.substr( 0, title.find(' ') ) ;
		string content = siter->substr( siter->find("\"") + 1, siter->rfind("\"") - siter->find("\"") - 1);

//		cout << "Title:\t\"" << title << "\"\tContent:\t" << content << endl;

		if ( title == "gene_id" )
			trans.gene_id = content;
		else if ( title == "transcript_id" )
			trans.transcript_id = content;
		else if ( title == "FPKM" )
			trans.fpkm = atof(content.c_str());
		else if ( title == "frac" )
			trans.frac = atof(content.c_str());
		else if ( title == "conf_lo" )
			trans.conf_lo = atof(content.c_str());
		else if ( title == "conf_hi" )
			trans.conf_hi = atof(content.c_str());
		else if ( title == "cov" )
			trans.cov = atof(content.c_str());
		else
			cout << "Unknow label:\t" << title << endl;
	}

	return 1;
}

struct Exon_Struct {
	string chr;
	string source;
	unsigned int left;
	unsigned int right;
	unsigned int score;
	char strand;
	char frame;
	string gene_id;
	string transcript_id;
	unsigned short exon_number;
	double fpkm;
	double frac;
	double conf_lo;
	double conf_hi;
	double cov;
};

int ExonParser ( Exon_Struct & exon, vector<string> & tokens) {
	exon.chr = tokens[0];
	exon.source = tokens[1];
	exon.left = atoi(tokens[3].c_str());
	exon.right = atoi(tokens[4].c_str());
	exon.score = atoi(tokens[5].c_str());
	exon.strand = tokens[6][0];
	exon.frame = tokens[7][0];
	
	string tmp = tokens[8];
	
	vector<string> comments = string_tokenize(tmp, ";");
	vector<string>::iterator siter;
	for ( siter = comments.begin(); siter != comments.end(); siter++) {
		string title = ( (*siter)[0] == ' ') ? siter->substr(1) : (*siter) ;
		title = title.substr( 0, title.find(' ') );
		string content = siter->substr( siter->find("\"") + 1, siter->rfind("\"") - siter->find("\"") - 1 );
		
//		cout << title << "\t:\t" << content << endl;

		if ( title == "gene_id" )
			exon.gene_id = content;
		else if ( title == "transcript_id" )
			exon.transcript_id = content;
		else if ( title == "FPKM" )
			exon.fpkm = atof(content.c_str());
		else if ( title == "frac" )
			exon.frac = atof(content.c_str());
		else if ( title == "conf_lo" )
			exon.conf_lo = atof(content.c_str());
		else if ( title == "cov" )
			exon.cov = atof(content.c_str());
		else if ( title == "exon_number" ) 
			exon.exon_number = atoi(content.c_str());
		else if ( title == "conf_hi")
			exon.conf_hi = atof(content.c_str());
		else
			cout << "Unknown label: \t" << title << endl;		
	}
	return 1;
}


int ReadGTFFile ( map<string, Chromosome> & Genome ) {
	ifstream gtfFile(param.gtfFile.c_str());
	if(!gtfFile) {
		cout << "cannot open input file" << param.gtfFile.c_str() << endl;
		return -1;
	}
	int k = 0;
	while (!gtfFile.eof()) {
		char buffer[1000+1];
		gtfFile.getline(buffer, 1000);
		if(!strlen(buffer))
			continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp, "\t");
		if (tokens.size() != 9) {
			cout << "The following line are not in RTF format:\n" << tmp << endl;
		}
		string feature = tokens[2];
		string chr = tokens[0];
		if ( feature == "transcript" ) {
			Transcript_Struct trans;
			TransParser ( trans, tokens );
			if ( Genome.find(chr) == Genome.end() ) {
				Chromosome NewChr;
				Genome[chr] = NewChr;
			}
			if ( Genome[chr].genes.find(trans.gene_id) == Genome[chr].genes.end() ) {
				Gene NewGene( trans.strand );
				Genome[chr].genes[trans.gene_id] = NewGene;			
			}
			if ( Genome[chr].genes[trans.gene_id].isoforms.find(trans.transcript_id) ==  Genome[chr].genes[trans.gene_id].isoforms.end() ) {
//				cout << "frac:\t" << trans.frac << endl;
				Transcript NewTrans( trans.left, trans.right, trans.fpkm, trans.frac );
				Genome[chr].genes[trans.gene_id].isoforms[trans.transcript_id] = NewTrans;
			} else {
				cerr << "Duplicate items: \t" << tmp << endl;
			}
			
		} else if ( feature == "exon") {
			Exon_Struct exon;
			ExonParser ( exon, tokens );
			if ( Genome.find(chr) == Genome.end() ) {
				cout << "Error logic happened.\n" << tmp << endl;
			} else if ( Genome[chr].genes.find(exon.gene_id) == Genome[chr].genes.end() ) {
				cout << "chr:\t" << chr << "\tGene_id:\t" << exon.gene_id << endl;
				cout << "No corresponding gene_id:\n" << tmp << endl;
			} else if ( Genome[chr].genes[exon.gene_id].isoforms.find(exon.transcript_id) == Genome[chr].genes[exon.gene_id].isoforms.end() ) {
				cout << "No corresponding transcript_id:\n" << tmp << endl;
			} else {
				Exon NewExon(exon.left, exon.right, exon.exon_number);
				Genome[chr].genes[exon.gene_id].isoforms[exon.transcript_id].exons.push_back( NewExon );
			}
			
		} else {
			cerr << "Unknown type of the new line!\n" << tmp << endl;
		}
	}
	gtfFile.close();
	return 0;
}

int OutputGTFFile ( map<string, Chromosome> & Genome ) {
	map<string, Chromosome>::iterator citer;
	string chr;	string type;	unsigned int left;	unsigned int right;
	char strand;	string gene_id;	string transcript_id;	double	FPKM;
	double	frac;

//	cout << "OutputGTFFile function.\n" << endl;

	for ( citer = Genome.begin(); citer != Genome.end(); citer++ ) {
		chr = citer->first;
		map<string, Gene>::iterator giter;
		for ( giter = Genome[chr].genes.begin(); giter != Genome[chr].genes.end(); giter++ ) {
			gene_id = giter->first;
			strand = giter->second.strand;
			map<string, Transcript>::iterator titer;
			for ( titer = giter->second.isoforms.begin(); titer != giter->second.isoforms.end(); titer++ ) {
				transcript_id = titer->first;
				cout << chr << "\tTranscript\t" << titer->second.left << "\t" << titer->second.right << "\t" << strand << "\tgene_id \"" << gene_id << "\"; trancript_id \"" << transcript_id << "\"; FPKM \"" << titer->second.fpkm << "\"; frac \"" << titer->second.fpkm << "\";"  << endl;
				
				vector<Exon>::iterator eiter;
				for ( eiter = titer->second.exons.begin(); eiter != titer->second.exons.end(); eiter++ ) {
					cout << chr << "\tExon\t" << eiter->left << "\t" << eiter->right << "\t" << strand << "\tgene_id \"" << gene_id << "\"; trancript_id \"" << transcript_id << "\"; exon_number \"" << eiter->exon_number << "\"; FPKM \"" << titer->second.fpkm << "\"; frac \"" << titer->second.frac << "\";"  << endl;	
					
				}
			}
		}
	}
	return 1;
}

int FindTwoIsoformsGene ( map< string, Chromosome > & Genome ) {
	int count = 0;

	map< string, Chromosome >::iterator citer;
	for ( citer = Genome.begin(); citer != Genome.end(); citer++ ) {
		map<string,Gene>::iterator giter;
		for ( giter = citer->second.genes.begin(); giter != citer->second.genes.end(); giter++ ) {
			if ( giter->second.isoforms.size() == 2 ) {
				count++;
				map<string,Transcript>::iterator titer, titer_1, titer_2;
				titer_1 = giter->second.isoforms.begin();
				titer_2 = titer_1; titer_2++;
				
				/*for ( titer = giter->second.isoforms.begin(); titer != giter->second.isoforms.end(); titer++ ) {
					cout << titer->second.exons.size() << "\t";
					
				}
				cout << endl;*/
				
				vector<Exon>::iterator eiter_1, eiter_2;
				eiter_1 = titer_1->second.exons.begin();
				eiter_2 = titer_2->second.exons.begin();
				
				while ( eiter_1 != titer_1->second.exons.end() && eiter_2 != titer_2->second.exons.end() && *eiter_1 == *eiter_2 ) {
					eiter_1++; eiter_2++;	
				}
				
				if ( eiter_1 != titer_1->second.exons.end() && eiter_2 != titer_2->second.exons.end() ) {
					if ( *eiter_1 > *eiter_2 && *eiter_1 == *(eiter_2 + 1) ) {
						/* Read 1:	==[  ]================[    ]======
						   Read 2:  ==[  ]====[    ]======[    ]======
						*/
					/*	// For read 1
						cout << "Read 1:\t[" << (eiter_1 - 1)->left << ",\t" << (eiter_1 - 1)->right << "]\t=>[\t\t\t]\t=> ["
							 << eiter_1->left << ",\t" << eiter_1->right << "]\n";
						// For read 2
						cout << "Read 2:\t[" << (eiter_2 - 1)->left << ",\t" << (eiter_2 - 1)->right << "]\t=>["
							 << eiter_2->left << ",\t" << eiter_2->right << "]\t=> ["
							 << (eiter_2 + 1)->left << ",\t" << (eiter_2 + 1)->right << "]\n"; 
						cout << endl;*/
						double ratio = titer_2->second.frac / (titer_1->second.frac + titer_2->second.frac);
						if ( giter->second.strand == '+' ) {
							cout << citer->first << "\t" << eiter_2->right << "\t" << (eiter_2 + 1)->left << "\t"
								 << giter->second.strand << "\t" << ratio << endl;
						} else {
							cout << citer->first << "\t" << (eiter_2 - 1)->right << "\t" << eiter_2->left << "\t"
								 << giter->second.strand << "\t" << ratio << endl;
						}
						
					} else if ( *eiter_1 < *eiter_2 && *(eiter_1 + 1) == *eiter_2 ) {
						/* Read 1:	==[  ]====[    ]======[    ]======
						   Read 2:  ==[  ]================[    ]======
						*/
					/*	// For read 1
						cout << "Read 1:\t[" << (eiter_1 - 1)->left << ",\t" << (eiter_1 - 1)->right << "]\t=>["
							 << eiter_1->left << ",\t" << eiter_1->right << "]\t=> ["
							 << (eiter_1 + 1)->left << ",\t" << (eiter_1 + 1)->right << "]\n"; 
						// For read 2
						cout << "Read 2:\t[" << (eiter_2 - 1)->left << ",\t" << (eiter_2 - 1)->right << "]\t=>[\t\t\t]\t=> ["
							 << eiter_2->left << ",\t" << eiter_2->right << "]\n";
						cout << endl;*/
						double ratio = titer_1->second.frac / (titer_1->second.frac + titer_2->second.frac);
//						cout << "1:\t" << titer_1->second.frac << "\t2:\t" << titer_2->second.frac << "\t3:\t" << ratio << endl;
						if ( giter->second.strand == '+' ) {
							cout << citer->first << "\t" << eiter_1->right << "\t" << (eiter_1 + 1)->left << "\t"
								 << giter->second.strand << "\t" << ratio << endl;
						} else {
							cout << citer->first << "\t" << (eiter_1 - 1)->right << "\t" << eiter_1->left << "\t"
								 << giter->second.strand << "\t" << ratio << endl;
						}
					}
				}
			}
		}
	} 
	
	return count;
}

int main(int argc, char* argv[])
{
	init();
//	cout << "Welcome\n";
#ifdef _DEBUG
	param.gtfFile = string("input.gtf");
#else
	parse_command_line(argc,argv);
#endif
	map< string, Chromosome > Genome;
	
	ReadGTFFile( Genome );
//	cout <<	"The number of genes have only two isoforms are:";
	FindTwoIsoformsGene ( Genome ) ;
//	cout << endl;
//	OutputGTFFile( Genome );
	return 1;
}


