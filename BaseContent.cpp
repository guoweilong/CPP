// BaseContent.cpp
// 2011-05-03
#include <iostream>
#include <fstream>

using namespace std;

struct parameter
{
        string	fastafile;	// -f
        bool	all;      	// -a
	bool	strand;		// -s
};

parameter param;

void exit_with_help( void )
{
        printf(
                "Author: GUO Weilong, 2011-05-03, guoweilong@gmail.com\n"
                "Usage: BaseContent -f <fasta.txt> [-a -s]\n"
                "options:\n"
                "-f Input, Fasta file, Sequences of the chromosome, each line\n"
                "     should be no longer than 1000. [ACGTNacgtn]\n"
                "-a switch, when specified, all the information will be output,\n"
                "     or else, only the frequency of the four bases\n"
		"-s switch, when specified, the information of the +/- strands are\n"
		"     showed, or else, only the + strand are showed\n"
                );
        exit(1);
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
                        if(++i>argc)    exit_with_help();
                        break;
                case 'a':
                        param.all = true;
                        break;
		case 's':
			param.strand = true;
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

void init(){
        // default values
        param.fastafile = "";
        param.all = false;
	param.strand = false;
}


int main(int argc, char* argv[])
{
        init();

#ifdef _DEBUG
        param.fastafile = string("hg18.fa");
#else
        parse_command_line(argc,argv);
#endif

	int i;
	long A, C, G, T, N, a, c, g, t, n, all;
	A = C = G = T = N = a = c = g = t = n = all = 0;
	ifstream infile;
	infile.open( param.fastafile.c_str() );
	int len = 1000;
	char* line = new char[len];
	while(infile.getline( line, len, '\n')){
		if(!strlen(line))
                        continue;
                if (line[strlen(line)-1] == '\r') {
                        line[strlen(line)-1] = '\0';
                }

		if( line[0] == '>' ){
			continue;
		}
		for ( i = 0; i < len && line[i]; i++ ){
			switch( line[i] ){
				case 'A':   A++;  break;
				case 'a':   a++;  break;
				case 'C':   C++;  break;
				case 'c':   c++;  break;
				case 'G':   G++;  break;
				case 'g':   g++;  break;
				case 'T':   T++;  break;
				case 't':   t++;  break;
				case 'N':   N++;  break;
				case 'n':   n++;  break;
				default:	break;
			}
			all++;
		}
	}
	infile.close();

	double total = A + C + G + T + a + c + g + t;
	if ( !param.all ) {
		cout << "A:\t" << (A+a)/total << endl << "C:\t" << (C+c)/total << endl << "G:\t" << (G+g)/total << endl << "T:\t" << (T+t)/total << endl;
	} else {
		if ( param.strand ) cout << "The + strand:\n";

		cout << "A:\t" << A << "\ta:\t" << a << "\tA/a:\t" << A+a << "\tRatio:\t" << (A+a)/total << endl;
		cout << "C:\t" << C << "\tc:\t" << c << "\tC/c:\t" << C+c << "\tRatio:\t" << (C+c)/total << endl;
		cout << "G:\t" << G << "\tg:\t" << g << "\tG/g:\t" << G+g << "\tRatio:\t" << (G+g)/total << endl;
		cout << "T:\t" << T << "\tt:\t" << t << "\tT/t:\t" << T+t << "\tRatio:\t" << (T+t)/total << endl;
		cout << "N:\t" << N << "\tn:\t" << n << "\tN/n:\t" << N+n << "\tN/known:\t" << (N+n)/total <<endl;
		cout << "Total(without N/n):\t" << (long) total << endl;
		cout << "Whole(with N/n):   \t" << (long) (total + N + n) << endl;
		
		if ( param.strand ) { 
			cout << "The - strand\n";
			cout << "A:\t" << T << "\ta:\t" << t << "\tA/a:\t" << T+t << "\tRatio:\t" << (T+t)/total << endl;
			cout << "C:\t" << G << "\tc:\t" << g << "\tC/c:\t" << G+g << "\tRatio:\t" << (G+g)/total << endl;
			cout << "G:\t" << C << "\tg:\t" << c << "\tG/g:\t" << C+c << "\tRatio:\t" << (C+c)/total << endl;
			cout << "T:\t" << A << "\tt:\t" << a << "\tT/t:\t" << A+a << "\tRatio:\t" << (A+a)/total << endl;
			cout << "N:\t" << N << "\tn:\t" << n << "\tN/n:\t" << N+n << "\tN/known:\t" << (N+n)/total << endl;
			cout << "Total(without N/n):\t" << (long) total << endl;
		}
	}
	return 1;
}



