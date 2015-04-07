#include "reads.h"
#include "smallapp.h"
#include "seq.h"

using namespace std;

int reads::countreads(char *inpfile, bool fastq)
{
	std::string s1;
        std::ifstream fp(inpfile);
        if(!fp) {
                cout<<"cannot open file "<<inpfile<<endl;
                exit(0);
        }

	int n = 0;
      	while(!fp.eof())     {
               	std::getline(fp, s1);
               	const char *str1 = s1.c_str();
                if(fastq || str1[0] == '>') n ++;
	}
	if(fastq && (n % 4 != 0)) {
		cout<<"wrong fastq format"<<endl;
	}
	if(fastq)	n = n / 4;
	fp.close();
	return(n);
}

int reads::get1readname(ifstream &fp, char *seqname)
{
	std::string s1;	
	while(!fp.eof()) {
		std::getline(fp, s1);
		const char *str = s1.c_str();
		if(str[0] == '>')	{
			sscanf(&str[1], "%s", seqname);
			return(1);
		}
	}
	return(0);
}

int reads::get1readfasta(ifstream &fp, char *seq, char *seqname)
{
	std::string s1;	
	int len_seq = 0;
	while(!fp.eof()) {
		std::getline(fp, s1);
		const char *str = s1.c_str();
		if(str[0] == '>')	{
			sscanf(&str[1], "%s", seqname);
			return(len_seq);
		} else	{
			for(int i = 0; i < s1.length(); i ++)	{
				if(str[i] <= 'z' && str[i] >= 'a')	{
					seq[len_seq ++] = seq::char2int(str[i]);
				} else if(str[i] <= 'Z' && str[i] >= 'A')	{
					seq[len_seq ++] = seq::char2int(str[i] - 'A' + 'a');
				}
			}
		}
	}
	return(len_seq);
}

int reads::get1readfastq(ifstream &fp, char *seq, char *seqname)
{
	std::string s1;	
	int len_seq = 0;
	if(!fp.eof()) {
		std::getline(fp, s1); //name
		const char *str = s1.c_str();
		sscanf(&str[1], "%s", seqname);

		std::getline(fp, s1); //sequence line
		str = s1.c_str();
		for(int i = 0; i < s1.length(); i ++)	{
			if(str[i] <= 'z' && str[i] >= 'a')	{
				seq[len_seq ++] = seq::char2int(str[i]);
			} else if(str[i] <= 'Z' && str[i] >= 'A')	{
				seq[len_seq ++] = seq::char2int(str[i] - 'A' + 'a');
			}
		}

		std::getline(fp, s1); 
		std::getline(fp, s1); //quality scores
	}
	return(len_seq);
}
