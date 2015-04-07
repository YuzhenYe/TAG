#include "seq.h"
#include "smallapp.h"
#include "contigapp.h"

int seq::readnum(char *inpfile, int *num, CONTIGNAMES *NameOfContigs, int nc, int num_contigs, int maxc)
{
	char	contig[100];
	std::string s;
	std::ifstream fp(inpfile);
	if(!fp) { cout<<"open file error "<<inpfile<<endl; exit(0); }
	int l = 0;
	while(!fp.eof()) {
		std::getline(fp, s); 
		sscanf(s.c_str(), "%*s%s%*d%*d", contig);
		int k = contigapp::findcontig(contig, NameOfContigs, nc, num_contigs, maxc);
		num[k] ++;
		l ++;
	}
	fp.close();
	return(l);
}

int seq::readfile(char *inpfile, SEGMENT **segment, int *num, CONTIGNAMES *NameOfContigs, int nc, int num_contigs, int maxc)
{
	int	pos[2];
	char	contig[100], name[100];
	int	l = 0;
	std::string s;
	std::ifstream fp(inpfile);
	if(!fp) { cout<<"open file error "<<inpfile<<endl; exit(0); }
	while(!fp.eof()) {
		std::getline(fp, s);
		if(s.empty()) continue;
		sscanf(s.c_str(), "%s%s%d%d", name, contig, &pos[0], &pos[1]);
		int k = contigapp::findcontig(contig, NameOfContigs, nc, num_contigs, maxc);
		strcpy(segment[k][num[k]].name, name);
		segment[k][num[k]].pos[0] = pos[0];
		segment[k][num[k]].pos[1] = pos[1];
		num[k] ++;
		l ++;
	}
	fp.close();
	return(l);
}

int seq::readpairfile(char *inpfile, char *pairfile, SEGMENT **segment, int *num, CONTIGNAMES *NameOfContigs, int nc, int num_contigs, int maxc)
{
	int	k1, k2, pos1[2], pos2[2];
	char	contig1[100], name1[100];
	char	contig2[100], name2[100];
	char	lastname1[100], lastname2[100];
	char	tag1, tag2;
	int	bit1, bit2;
	int	index1, index2;

	std::ifstream fp1(inpfile);
	std::ifstream fp2(pairfile);
	std::string s1, s2;
	std::getline(fp1, s1);
	std::getline(fp2, s2);
	const char* str1 = s1.c_str();
	const char* str2 = s2.c_str();
	tag1 = tag2 = 1;
	bit1 = bit2 = 0;
	strcpy(lastname1, "\0");
	strcpy(lastname2, "\0");
	int l = 0;
	while(!fp1.eof() || !fp2.eof())      {
		if(!s1.empty() && tag1)	{
			sscanf(str1, "%s%s%d%d", name1, contig1, &pos1[0], &pos1[1]);
			k1 = contigapp::findcontig(contig1, NameOfContigs, nc, num_contigs, maxc);
			index1 = parseindex(name1, lastname1, &bit1);
			strcpy(lastname1, name1);
			tag1 = 0;
		}
		if(!s2.empty() && tag2)	{
			sscanf(str2, "%s%s%d%d", name2, contig2, &pos2[0], &pos2[1]);
			k2 = contigapp::findcontig(contig2, NameOfContigs, nc, num_contigs, maxc);
			index2 = parseindex(name2, lastname2, &bit2);
			strcpy(lastname2, name2);
			tag2 = 0;
		}
		if(!s1.empty() && !s2.empty())	{
			if(index1 == index2)	{
				if(k1 == k2)	{
					strcpy(segment[k1][num[k1]].name, name1);
					segment[k1][num[k1]].pos[0] = min(pos1[0], pos2[0]);
					segment[k1][num[k1]].pos[1] = max(pos1[1], pos2[1]);
					num[k1] ++;
					l ++;
				} else	{
					strcpy(segment[k1][num[k1]].name, name1);
					segment[k1][num[k1]].pos[0] = pos1[0];
					segment[k1][num[k1]].pos[1] = pos1[1];
					num[k1] ++;
					l ++;
					strcpy(segment[k2][num[k2]].name, name2);
					segment[k2][num[k2]].pos[0] = pos2[0];
					segment[k2][num[k2]].pos[1] = pos2[1];
					num[k2] ++;
					l ++;
				}
				std::getline(fp1, s1);
				str1 = s1.c_str();
				std::getline(fp2, s2);
				str2 = s2.c_str();
				tag1 = tag2 = 1;
			} else if(index1 < index2)	{
				strcpy(segment[k1][num[k1]].name, name1);
				segment[k1][num[k1]].pos[0] = pos1[0];
				segment[k1][num[k1]].pos[1] = pos1[1];
				num[k1] ++;
				l ++;
				std::getline(fp1, s1);
				str1 = s1.c_str();
				tag1 = 1;
			} else {
				strcpy(segment[k2][num[k2]].name, name2);
				segment[k2][num[k2]].pos[0] = pos2[0];
				segment[k2][num[k2]].pos[1] = pos2[1];
				num[k2] ++;
				l ++;
				std::getline(fp2, s2);
				str2 = s2.c_str();
				tag2 = 1;
			}
		} else if(!fp1.eof())	{
			strcpy(segment[k1][num[k1]].name, name1);
			segment[k1][num[k1]].pos[0] = pos1[0];
			segment[k1][num[k1]].pos[1] = pos1[1];
			num[k1] ++;
			l ++;
			std::getline(fp1, s1);
			str1 = s1.c_str();
			tag1 = 1;
		} else	{
			strcpy(segment[k2][num[k2]].name, name2);
			segment[k2][num[k2]].pos[0] = pos2[0];
			segment[k2][num[k2]].pos[1] = pos2[1];
			num[k2] ++;
			l ++;
			std::getline(fp2, s2);
			str2 = s2.c_str();
			tag2 = 1;
		}
	}
	fp1.close();
	fp2.close();
	return(l);
}

int seq::parseindex(char *name, char *lastname, int *bit)
{
	int	i, l, k;
	int	index;

	l = strlen(name);
	for(i = 0; i < l; i ++)	{
		if(name[i] == '.')	{
			if(strncmp(name, lastname, i))	{
				*bit = *bit + 1;
			}
			index = atoi(&name[i + 1]) + *bit * LARGENUMBER;
			break;
		}
	}
	if(i == l)	{
		printf("name format error: %s\n", name);
		exit(-1);
	}
	return(index);
}

int seq::readseq_num(char *filename)
{
	std::ifstream fp(filename);
	if(!fp) { cout<<"open file error "<<filename<<endl; exit(0); }
	int n = 0;
	std::string s;
	while(!fp.eof()) {
		std::getline(fp, s);
		const char* str = s.c_str();
		if(str[0] == '>')	n ++;
	}
	fp.close();
	return(n);
}
				
int seq::readseq(char **src_seq, char **src_name, int *len_seq, char *filename)
{
	int	i, k, n, len;
        char   *buf;

	std::ifstream fp(filename);
	if(!fp) { cout<<"open file error "<<filename; exit(0); }
        buf = (char *)smallapp::ckalloc(sizeof(char) * 1000000);

	n = 0;
	k = -1;
	std::string s;
	while(!fp.eof()) {
		std::getline(fp, s);
		const char *str = s.c_str();
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			if(k >= 0)
                        {
                          src_seq[k] = (char *)smallapp::ckalloc(sizeof(char) * (n + 1));
                          memcpy((void *)src_seq[k], (void *)buf, sizeof(char) * (n + 1));
                          len_seq[k] = n;
                        }
			n = 0;
			k ++;
			sscanf(&str[1], "%s", buf);
                        src_name[k] = (char *)smallapp::ckalloc(sizeof(char) * (strlen(buf) + 1));
                        strcpy(src_name[k], buf);
		} else {
			for(i = 0; i < s.length(); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					buf[++ n] = char2int(str[i]);
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					buf[++ n] = char2int(str[i]
						 - 'A' + 'a');
				}
			}
		}
	}
	fp.close();
        
        if (k >= 0)
        {
          src_seq[k] = (char *)smallapp::ckalloc(sizeof(char) * (n + 1));
          memcpy((void *)src_seq[k], (void *)buf, sizeof(char) * (n + 1));
          len_seq[k] = n;
        }
	k ++;

        free((void *)buf);
	return(k);
}

char seq::char2int(char c)
{
	int	i, k;

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	int	idum = 1234; //to be checked YY
	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = smallapp::ran_number(4, &idum);
	}

	if(k > 3)	{
		k = smallapp::ran_number(4, &idum);
	}

	return(k);
}

char seq::char2intgen(char c)
{
	int	i, k;

	for(i = 0; i < total_nuc; i ++)		{
		if(c == na_name[i])	{
			k = i;
			break;
		}
	}

	int idum = 1234; //to-be-checked by YY
	if(i == total_nuc)	{
		printf("Not found %c\n", c);
		k = smallapp::ran_number(4, &idum);
	}

	return(k);
}


int seq::outputseq(ofstream &fp, char *seq, char *name, int pbeg, int pend)
{
	int i, j, l;
	char str[1000];

	sprintf(str, ">%s-%d-%d %d\n", name, pbeg, pend, pend-pbeg+1);
	fp<<str;
	l = 0; 
	for(i = pbeg; i <= pend; i ++)	{
		fp<<na_name[seq[i]];
		l ++;
		if(l == 60)	{
			l = 0;
			fp<<endl;
		}
	}
	if(l != 60)	fp<<endl;
	return(1);
}

int seq::writeseq(ofstream &fp, char *seq, char *name, int len_seq)
{
	int i, j, l;

	fp<<">"<<name<<endl;
	l = 0; 
	for(i = 0; i < len_seq; i ++)	{
		fp<<na_name[seq[i]];
		l ++;
		if(l == 60)	{
			l = 0;
			fp<<endl;
		}
	}
	if(l != 60)	fp<<endl;
	return(1);
}
