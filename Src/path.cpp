#include "path.h"
#include "smallapp.h"

EDGESEG* path::insert_single_path(EDGESEG *pathlist, EDGESEG *newpathlist)
{
	newpathlist -> length = 1;
	newpathlist -> num_reads = 1;
	newpathlist -> total_length += (newpathlist -> pos[1] - newpathlist -> pos[0] + 1);
	newpathlist -> next = pathlist;
	return(newpathlist);
}

int path::ReadSingleSegment(bool stranded, char *inpfile, EDGESEG **pathlist, EDGE *edge, int num_edge)
{
	int	i, j, k, l, index, bal_index, pos[2], cpos[2], cindex;
	char	contig[100], strand, strandtmp[10];
	ifstream fp(inpfile);
	if(!fp) { cout<<"open file error "<<inpfile<<endl; exit(0); }

	std::string s;
	k = 0;
	while(!fp.eof()) {
		std::getline(fp, s);
		if(s.empty()) continue;
		const char* str = s.c_str();
		sscanf(str, "%*s%d%s%d%d", &cindex, strandtmp, &(pos[0]), &(pos[1]));
		strand = smallapp::strand2char(strandtmp);
		pos[1] --; /* note: the end mapped position (pos[1]) is one larger than the last mapped position        */
		cindex --;
		bal_index = edge[cindex].bal_edge;	/* complementary edge	*/
		cpos[0] = edge[cindex].length + 1 - pos[0];
		cpos[1] = edge[cindex].length + 1 - pos[1];
		EDGESEG *newpathlist = smallapp::newEDGESEGpt();
		if(strand == 0)	{	/* add the segment where the read is mapped to positive strand	*/
			newpathlist -> pos[0] = pos[0];
			newpathlist -> pos[1] = pos[1];
			newpathlist -> index = cindex;
		} else	{
			newpathlist -> pos[0] = cpos[0];
			newpathlist -> pos[1] = cpos[1];
			newpathlist -> index = bal_index;
		}
		pathlist[newpathlist -> index] = insert_single_path(pathlist[newpathlist -> index], newpathlist); 
		k ++;
		if((!stranded) && bal_index != cindex)	{ // unstranded protocol: both forward and reverse strand will be added
			newpathlist = smallapp::newEDGESEGpt();
			if(strand == 0)	{// add the segment where the read is mapped to positive strand
				newpathlist -> pos[0] = cpos[0];
				newpathlist -> pos[1] = cpos[1];
				newpathlist -> index = bal_index;
			} else	{
				newpathlist -> pos[0] = pos[0];
				newpathlist -> pos[1] = pos[1];
				newpathlist -> index = cindex;
			}
			newpathlist -> index = bal_index;
			newpathlist -> pos[0] = edge[bal_index].length + 1 - pos[0];
			newpathlist -> pos[1] = edge[bal_index].length + 1 - pos[1];
			pathlist[bal_index] = insert_single_path(pathlist[bal_index], newpathlist); 
			k ++;
		}
	}
	fp.close();
	return(k);
}

int path::countseg(char *inpfile)
{
	int k = 0;
	std::string s;
	ifstream fp(inpfile);
	if(!fp) { cout<<"open file error "<<inpfile<<endl; exit(0); }
	while(!fp.eof()) {
		std::getline(fp, s);
		if(!s.empty()) k ++;
	}
	fp.close();
	return(k);
}

int path::ReadPairSegment(bool stranded, int insert_size, char *inpfile, char *pairfile, EDGESEG **pathlist, EDGE *edge, int num_edge)
{
	int	l, pos1[2], pos2[2], pt[2];
	int	distance;
	EDGESEG	*newpathlist;
	char	name1[100], name2[1000];
	char	lastname1[100], lastname2[100];
	char	tag1, tag2, strand1, strand2, jointstrand, strandtmp[10];
	int	bit1, bit2;
	int	index1, index2;
	int	cindex1, cindex2;
	int	bal_index;
	l = 0;
	std::ifstream fp1(inpfile);
	std::ifstream fp2(pairfile);
	if(!fp1) { cout<<"read file error "<<inpfile<<endl; exit(0); }
	if(!fp2) { cout<<"read file error "<<pairfile<<endl; exit(0); }
	std::string s1, s2;
	std::getline(fp1, s1);
	const char* str1 = s1.c_str();
	std::getline(fp2, s2);
	const char* str2 = s2.c_str();
	tag1 = tag2 = 1;
	bit1 = bit2 = 0;
	strcpy(lastname1, "\0");
	strcpy(lastname2, "\0");
	while(!s1.empty() || !s2.empty()) {
		if(!s1.empty() && tag1)	{
			sscanf(str1, "%s%d%s%d%d", name1, &cindex1, strandtmp, &pos1[0], &pos1[1]);
			strand1 = smallapp::strand2char(strandtmp);
			pos1[1] --; /* note: the end mapped position (pos1[1]) is one larger than the last mapped position        */
			index1 = parseindex(name1, lastname1, &bit1);	/*	read index	*/
			cindex1 --;
			strcpy(lastname1, name1);
			tag1 = 0;
		}
		if(!s2.empty() && tag2)	{
			sscanf(str2, "%s%d%s%d%d", name2, &cindex2, strandtmp, &pos2[0], &pos2[1]);
			strand2 = smallapp::strand2char(strandtmp);
			pos2[1] --; /* note: the end mapped position (pos2[1]) is one larger than the last mapped position        */
			index2 = parseindex(name2, lastname2, &bit2);	/*	read index	*/
			cindex2 --;
			strcpy(lastname2, name2);
			tag2 = 0;
		}
		if(!fp1.eof() && !fp2.eof())	{	/*	both reads file are not empty	*/
			if(index1 == index2)	{	/*	read pairs	*/
				if(cindex1 == cindex2 && strand1 != strand2)	{	/* must be in complement strand	*/
					pt[0] = min(pos1[0], pos2[0]);
					pt[1] = max(pos1[1], pos2[1]);
					if(pt[0] == pos1[0] && strand1 == 0 && strand2 == 1)	{
						jointstrand = 1;		/*	mate-pairs are reverse complementary	*/
						distance = pt[1] - pt[0] + 1;
					} else if(pt[0] == pos2[0] && strand2 == 0 && strand1 == 1)	{
						jointstrand = 0;		/*	mate-pairs are forward	*/
						distance = pt[1] - pt[0] + 1;
					} else	{
						distance = insert_size + 100;
					}
				} else	{
					distance = insert_size + 100;
				}
				if(distance <= insert_size)	{		
					//both reads mapped to the same edge and within a distance of insert_size: they will form a single segment
					l += ProcessOneAln(pathlist, edge, stranded, cindex1, jointstrand, pt);	
				} else	{	/*	Each read mapped to a different edge: they will form two separate segments	*/
					l += ProcessOneAln(pathlist, edge, stranded, cindex1, !strand1, pos1);
					l += ProcessOneAln(pathlist, edge, stranded, cindex2, strand2, pos2);
				}
				std::getline(fp1, s1);
				str1 = s1.c_str();
				std::getline(fp2, s2);
				str2 = s2.c_str();
				tag1 = tag2 = 1;
			} else if(index1 < index2)	{
				l += ProcessOneAln(pathlist, edge, stranded, cindex1, !strand1, pos1);
				std::getline(fp1, s1);
				str1 = s1.c_str();
				tag1 = 1;
			} else {
				l += ProcessOneAln(pathlist, edge, stranded, cindex2, strand2, pos2);
				std::getline(fp2, s2);
				str2 = s2.c_str();
				tag2 = 1;
			}
		} else if(!fp1.eof())	{
			l += ProcessOneAln(pathlist, edge, stranded, cindex1, !strand1, pos1);
			std::getline(fp1, s1);
			str1 = s1.c_str();
			tag1 = 1;
		} else	{
			l += ProcessOneAln(pathlist, edge, stranded, cindex2, strand2, pos2);
			std::getline(fp2, s2);
			str2 = s2.c_str();
			tag2 = 1;
		}
	}
	fp1.close();
	fp2.close();
	return(l);
}

/*
SAM FLAG
Bit    Description
0x1    template having multiple segments in sequencing
0x2    each segment properly aligned according to the aligner
0x4    segment unmapped
0x8    next segment in the template unmapped
0x10   SEQ being reverse complemented (5th bit)
and so on
Each flag is a bit, containg a 0 or a 1, representing no and yes respectively. 
how to check the strand?
If the 5th bit (0x10) is 1 then the sequence is reverse complemented; if it is 0 than it is not reverse complemented. 
A bitwise flag where only the 5th bit is set would have a binary string of:
binary 10000 = 1 * 24 + 0 * 23 + 0 * 22 + 0 * 21 + 0 * 20 = decimal 16
e.g., to check the strand: int(flag) & 16:
 */
bool path::ReadSamOneLine(const char *str, int *pos, int &cindex, char &strand, int &mate)
{
//printf("str = %s\n", str);
	char	qname[100], rname[100], cigar[100], flagstr[100], posstr[100], tmp[100];
	int	i, flag, d;
	sscanf(str, "%s%s%s%s%s%s", qname, flagstr, rname, posstr, tmp, cigar);
	sscanf(flagstr, "%d", &flag);
	sscanf(posstr, "%d", &pos[0]);
	sscanf(rname, "%d", &cindex); //index as the seq-id
//printf("cindex %d\n", cindex);
	int 	n0 = 0;
	int	mlen, ins, dele, clip;
	mlen = ins = dele = clip = 0;	
	for(i = 0; i < strlen(cigar); i ++) {
		char c = cigar[i];
		if(c=='M' || c=='I' || c=='D' || c=='N' || c=='S' || c=='H' || c=='P' || c=='X') {
			sscanf(cigar + n0, "%d", &d);
			if (c == 'M') mlen += d; 
			else if (c == 'I') ins += d;
			else if (c == 'D') dele += d;
			else if (c == 'S') clip += d;
			n0 = i + 1;
		}
	}

	int mismatch = 0;
	for(i = 0; i < strlen(str); i ++) {
		if(!strncmp(str + i, "NM:i:", 5)) {
			sscanf(str + i + 5, "%d", &mismatch);
		}
	}
        if(mismatch + dele + ins > MAXD) return false; 
	int ref_alnlen = mlen + dele;
	pos[1] = pos[0] + ref_alnlen - 1; //-1
	//bitwise operation applied to mask (16 decimal, in binary 0001 0000)
        if (int(flag) & 16) strand = 1; //strandtmp = '-'; check smallapp::strand2char(strandtmp); reverse complemented
	else strand = 0; //strandtmp = '+';
//printf("pos %d %d strand %d\n", pos[0], pos[1], strand);
//getchar();
	if(int(flag) & 64) mate = 1; //7th bit -- the first segment in the template
	else mate = 2;
	if(ref_alnlen < 30) return false;
	else return true;
}

//input: cindex: contig-index; strand: 0 for reverse, 1 for forward; alnpos: alignment-position
int path::ProcessOneAln(EDGESEG **pathlist, EDGE *edge, bool stranded, int cindex, char strand, int *alnpos)
{
	EDGESEG* newpathlist;
	newpathlist = smallapp::newEDGESEGpt();
	int bal_index = edge[cindex].bal_edge; //index of the complementary strand
	if(strand == 0)	{ //forward strand +
		newpathlist -> index = cindex;
		newpathlist -> pos[0] = alnpos[0];
		newpathlist -> pos[1] = alnpos[1];
		pathlist[cindex] = insert_single_path(pathlist[cindex], newpathlist); 
	} else	{ //reverse strand -
		newpathlist -> index = bal_index;
		newpathlist -> pos[0] = edge[bal_index].length + 1 - alnpos[1];
		newpathlist -> pos[1] = edge[bal_index].length + 1 - alnpos[0];
		pathlist[bal_index] = insert_single_path(pathlist[bal_index], newpathlist); 
	}

	if(!stranded)	{ //if RNAseq was unstranded -- consider the reverse complementary chain as well 
		newpathlist = smallapp::newEDGESEGpt();	
		if(strand != 0)	{
			newpathlist -> index = cindex;
			newpathlist -> pos[0] = alnpos[0];
			newpathlist -> pos[1] = alnpos[1];
			pathlist[cindex] = insert_single_path(pathlist[cindex], newpathlist); 
		} else	{
			newpathlist -> index = bal_index;
			newpathlist -> pos[0] = edge[bal_index].length + 1 - alnpos[1];
			newpathlist -> pos[1] = edge[bal_index].length + 1 - alnpos[0];
			pathlist[bal_index] = insert_single_path(pathlist[bal_index], newpathlist); 
		}
		return 2;
	}
	else return 1;
}

int path::ReadPairedSam(bool stranded, int insert_size, char *samfile, EDGESEG **pathlist, EDGE *edge, int num_edge)
{
	int	i, j, k, k1, pos1[2], pos2[2], pt[2];
	int	distance;
	char	name1[1000], name2[1000];
	char	strand1, strand2, jointstrand, strandtmp[10];
	int	index1, index2, fl1, fl2, mate1, mate2;
	int	cindex1, cindex2;
	bool	merge, r1, r2;
	std::string s1, s2;

	cout<<"\n>>>Load sam file "<<samfile<<endl;
	std::ifstream infile(samfile);
	if(!infile) {
		cout<<"cannot open file "<<samfile<<endl; 
		exit(0);
	}
	int paired = 0;
	int pair_con = 0;
	int pair_con_merge = 0;
	int pair_noncon = 0;
	int single = 0;
	int considered = 0;
	int tot = 0;
	int wrongmate = 0;

	//important: for strand-specific protocol, the 2nd read comes
        //from the forward strand, and the 1st read comes from the rt strand;
        //as such, if mate #1 is mapped to the forward strand, it is actually from the rt strand (-)
        //don't change the direction for mate #2
	while(!infile.eof())     {
		std::getline(infile, s1);
		const char *str1 = s1.c_str();
		if(str1[0] == '@') continue;
		sscanf(str1, "%s%d", name1, &fl1);
		r1 = ReadSamOneLine(str1, pos1, cindex1, strand1, mate1);
		if(fl1==73 || fl1==153 || fl1==137 || fl1==89 || fl1==129 || fl1==147 || fl1==161 || fl1==177 || fl1==145 || fl1==163) { //only one mate mapped
			single += 1;
			if(!r1) continue;
			considered += 1;
			//mate #1
			if (mate1 == 1) tot += ProcessOneAln(pathlist, edge, stranded, cindex1, !strand1, pos1); //mate #1, note !strand1 (change the direction)
			else tot += ProcessOneAln(pathlist, edge, stranded, cindex1, strand1, pos1); //mate #2, don't change the direction
			//mate #2
		}
		else if(fl1==65||fl1==99||fl1==81||fl1==113||fl1==97||fl1==83) {//both mates mapped, concordantly or non-concordantly
			std::getline(infile, s2);
			const char *str2 = s2.c_str();
			sscanf(str2, "%s%d", name2, &fl2);
			r2 = ReadSamOneLine(str2, pos2, cindex2, strand2, mate2);
			if(!((mate1 == 1) && (mate2 == 2))) wrongmate += 1; //unknown reason!!!
			//printf("Read1 %s mate %d strand %d Read2 %s mate %d strand %d\n", name1, mate1, strand1, name2, mate2, strand2);
			//getchar();
			paired += 1;
			if(fl1 == 99 || fl1 == 83) pair_con += 1;
			merge = false;
			if(r1 && r2 && (mate1 == 1) && (mate2 == 2) && ((fl1==99 || fl1==83)))	{ //concordant
				if(strand1 != strand2)	{//must be in complement strand
					pt[0] = min(pos1[0], pos2[0]);
					pt[1] = max(pos1[1], pos2[1]);
					if(pt[0] == pos1[0] && strand1 == 0 && strand2 == 1)	{
						jointstrand = 1;	//mate-pairs are reverse complementary
						distance = pt[1] - pt[0] + 1;
						merge = true;
					} else if(pt[0] == pos2[0] && strand2 == 0 && strand1 == 1)	{
						jointstrand = 0;	//mate-pairs are forward
						distance = pt[1] - pt[0] + 1;
						merge = true;
					}
				}
			}
			if(merge) {
				if(distance > 1000) {
					cout<<"Warning: mates too far apart "<<distance<<endl;
					getchar();
				}
				pair_con_merge += 1;
				//paired reads mapped concordantly to the same edge: they will form a single segment
				tot += ProcessOneAln(pathlist, edge, stranded, cindex1, jointstrand, pt); //ok
			}
			else	{ // don't merge
				//mate #1
				if(r1 && (mate1 == 1)) {
					considered += 1;
					tot += ProcessOneAln(pathlist, edge, stranded, cindex1, !strand1, pos1); 
					//note: !strand1, change direction for mate 1
					//else unknown problem with sam (some cases mate1 == 2)
				}
				//mate #2
				if(r2 && (mate2 == 2)) {
					considered += 1;
					tot += ProcessOneAln(pathlist, edge, stranded, cindex2, strand2, pos2); 
					//note: just strand2, keep the direction for mate 2
					//else unknown problem with sam
				}
			}
		}
	}
	infile.close();
	cout<<"Total pairs "<<paired<<endl;
	cout<<"Total pairs-concordant "<<pair_con<<endl;
	cout<<"Total pairs-concordant-merged "<<pair_con_merge<<endl;
	cout<<"Total single mates "<<single<<endl;
	cout<<"Single mates or pairs but considered separately "<<considered<<endl;
	cout<<"Warning: wrong mate "<<wrongmate<<endl;
	return(tot);
}

int path::parseindex(char *name, char *lastname, int *bit)
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

EDGESEG* path::insert_path(EDGESEG *path, int len_path, EDGESEG *pathlist)
{
	int	i;
	//EDGESEG *newpathlist;
	//newpathlist = (EDGESEG *) smallapp::ckalloc(len_path * sizeof(EDGESEG));
	EDGESEG *newpathlist = new EDGESEG[len_path]; //pointer only!!
	for(i = 0; i < len_path; i ++)	{
		newpathlist[i] = path[i];
	}
	newpathlist -> length = len_path;
	newpathlist -> next = pathlist;

	return(newpathlist);
}

int path::reverse_path(EDGESEG *path, EDGESEG *rev_path, int len_path, EDGE *edge)
{
	int	i, j, k, l, bal_edge, length;

	for(i = 0; i < len_path; i ++)	{
		length = edge[path[i].index].length;
		bal_edge = edge[path[i].index].bal_edge;
		rev_path[len_path - i - 1].index = bal_edge;
		rev_path[len_path - i - 1].pos[0] = length - path[i].pos[1] + 1;
		rev_path[len_path - i - 1].pos[1] = length - path[i].pos[0] + 1;
		rev_path[len_path - i - 1].num_reads = path[i].num_reads;
		rev_path[len_path - i - 1].total_length = path[i].total_length;
	}
	rev_path -> length = len_path;
}

void path::delete_path(EDGESEG *pathlist)
{
	EDGESEG *pathlistnext;

	if(!pathlist)	return;

	pathlistnext = pathlist -> next;
	//free((void *) pathlist);
	delete[] pathlist;
	delete_path(pathlistnext);
}
