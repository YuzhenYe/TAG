//Package: tag
//Program Name: seqgraph (the main program for mapping reads to de Bruijn graph) 
//Latest modification on Oct 28th, 2014
//Yuzhen Ye <yye@indiana.edu> and Haixu Tang <hatang@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington
#include "seqgraph.h"
#include "smallapp.h"
#include "soap.h"
#include "reads.h"
#include "path.h"

int ifprint = 0;

seqgraph::seqgraph(void) {
	init();
}

//initialize
void seqgraph::init() {
	inppos = 0;
	ngap = 3;
	maxmut = 3;
	kmersize = 27;
	overlaplen = 20;
	stranded = true;

	hashw = 12;
	hashtablesize = smallapp::hashvalue(hashw);
	set_maskv();
}

void seqgraph::set_maskv()
{
        maskv = 1;
        maskv = maskv << (2 * kmersize);
        maskv = maskv - 1;
}

void seqgraph::cleanup(void)
{
	delete[] vindex;

	int i;
	for(i = 0; i < num_edge; i ++)	{
		path::delete_path(pathlist[i]);
	}
	delete[] pathlist;

	for(i = 0; i < num_edge; i ++)	{
		delete[] edge[i].seq;
	}
	delete[] edge;
	for(i = 0; i < num_vertex; i ++)	{
		soap::delete_edgelist(vertex[i].lastedge);
		soap::delete_edgelist(vertex[i].nextedge);
	}
	delete[] vertex;
}

//load the graph (from SOAPdenovo)
void seqgraph::loadsoap(bool SOAP2_tag, char *edgefile, char *edgeseqfile) 
{
	int i;
	//get num of edges
	cout<<"\n>>>Now read graph from file "<<edgefile<<endl;
	num_edge = soap::countedge(edgefile);
	cout<<"num edge to be allocated: "<<num_edge<<endl;

	//allocate vertex
	vertex = new VERTEX[2 * (num_edge + 1)];
	//initialize the variables -- it is important!
	cout<<"num_vertex to be allocated: "<<2 * (num_edge + 1)<<endl;
	for(i = 0; i < 2 * (num_edge + 1); i ++) {
		vertex[i].lastedge = vertex[i].nextedge = NULL;
		vertex[i].indegree = vertex[i].outdegree = 0;
	}
	num_vertex = soap::readvertex(SOAP2_tag, edgefile, vertex);
	cout<<"num_vertex actually extracted:  "<<num_vertex<<endl;

	//Initialize class attributes (global) -- edges for contigs
	edge = new EDGE[num_edge]; 
	for(i = 0; i < num_edge; i ++) { 
		edge[i].tag = edge[i].begtag = edge[i].endtag = edge[i].parttag = 0;
		edge[i].firstseg = edge[i].lastseg = edge[i].firstseg_num_reads = edge[i].lastseg_num_reads = 0;
		edge[i].firstseg_total_length = edge[i].lastseg_total_length = 0;
		edge[i].num_reads = edge[i].total_length = edge[i].index = edge[i].length = edge[i].bal_edge = 0;
		edge[i].seq = NULL;
		edge[i].nodeindex[0] = edge[i].nodeindex[1] = 0;
	}

	cout<<"SOAP2_tag "<<SOAP2_tag<<" kmersize "<<kmersize<<endl;
	num_edge = soap::readedge(SOAP2_tag, edgefile, edgeseqfile, edge, vertex, num_vertex, kmersize);
	cout<<"input edge extracted: "<<num_edge<<endl;

	for(i = 0; i < num_edge; i ++)	{
		edge[i].index = i;
		edge[i].lastseg = edge[i].length + 1;
	}
	cout<<"print_vertex_degree.."<<endl;
	print_vertex_degree();

	cout<<"hashtablesize "<<hashtablesize<<endl;
	//vindex = (int *) smallapp::ckalloc((hashtablesize + 1) * sizeof(int));
	vindex = new int[hashtablesize + 1];
	for(i = 0; i < hashtablesize + 1; i ++) vindex[i] = 0;
	vindex[hashtablesize] = num_vertex;

	cout<<"index vertex..."<<endl;
	index_vertex();

	for(i = hashtablesize - 1; i >= 0; i --)	{
		if(vindex[i] == 0)	{
			vindex[i] = vindex[i + 1];
		}
	}
}

void seqgraph::index_vertex(void)
{
	int shift = (kmersize - hashw) * 2;
	printf("kmer %d hashw %d shift %d\n", kmersize, hashw, shift);

	int hashv = vertex[0].index >> shift;
	vindex[hashv] = 1;
	//cout<<"index "<<vertex[0].index<<" hashv "<<hashv<<endl;
	printf("index %ld hashv %d\n", vertex[0].index, hashv);
	for(int i = 1; i < num_vertex; i ++)	{
		hashv = vertex[i].index >> shift;
		if(vindex[hashv] == 0)	{
			vindex[hashv] = i + 1;
		}
	}
}

void seqgraph::print_vertex_degree(void)
{
	int	i, j, n1, n2;
	int	vertex_in_out[5][5];

	for(i = 0; i < 5; i ++)	{
		for(j = 0; j < 5; j ++)	{
			vertex_in_out[i][j] = 0;
		}
	}
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i].indegree > 4)	{
			n1 = 4;
		} else	{
			n1 = vertex[i].indegree;
		}
		if(vertex[i].outdegree > 4)	{
			n2 = 4;
		} else	{
			n2 = vertex[i].outdegree;
		}
		vertex_in_out[n1][n2] ++;
	}
	printf("Degree distribution..\n");
	printf("-----------------------\n");
	for(i = 0; i < 5; i ++)	{
		for(j = 0; j < 5; j ++)	{
			printf("%d ", vertex_in_out[i][j]);
		}
		printf("\n");
	}
	printf("-----------------------\n");
}

//paired sam file 
void seqgraph::loadsam(char *samfile, int insert_size)
{
	int	i;
	pathlist = new EDGESEG*[num_edge];
	for(i = 0; i < num_edge; i ++) pathlist[i] = NULL;
	int tot_seg = path::ReadPairedSam(stranded, insert_size, samfile, pathlist, edge, num_edge);
	cout<<"Number of alignment regions considered "<<tot_seg<<endl;

	int n = 0;
	for(i = 0; i < num_edge; i ++)	{
		if(pathlist[i])	{
			n ++;
		}
	}
	cout<<n<<" out of "<<num_edge<<" edges covered by RNA-seq reads."<<endl;
}

//segfile & pairfile -- mapping information (of the reads against contigs/edges by the conventional mapper)
void seqgraph::loadmap(char *segfile, char *pairfile, int insert_size)
{
	int 	i, tot_seg;

	cout<<"input segments in the same edge\n";

	//pathlist records the regions (begin: pos[0], end: pos[1]) in each edge (contig), to which reads are mapped
	pathlist = new EDGESEG*[num_edge];

	for(i = 0; i < num_edge; i ++) pathlist[i] = NULL;
	if(pairfile[0] != 0)	{
		cout<<"paired file of mapped reads"<<endl;
		tot_seg = path::ReadPairSegment(stranded, insert_size, segfile, pairfile, pathlist, edge, num_edge);
	} else	{
		tot_seg = path::ReadSingleSegment(stranded, segfile, pathlist, edge, num_edge);
	}
	cout<<"Number of segment input "<<tot_seg<<endl;

	int n = 0;
	for(i = 0; i < num_edge; i ++)	{
		if(pathlist[i])	{
			n ++;
		}
	}
	cout<<n<<" out of "<<num_edge<<" edges covered by RNA-seq reads."<<endl;
}

//map reads to the graph -- read2file will be empty if not paired
//read1file & read2file -- reads file
void seqgraph::read2graph(char *read1file, char *read2file, bool fastq)
{
printf("ngap %d\n", ngap);
	//Mapping to subgraph with edges of tag =  1
	cout<<"\n>>>Now map the reads to graph.."<<endl;
	//EDGESEG* readpath = (EDGESEG *) smallapp::ckalloc(MIN_LEG * sizeof(EDGESEG *));
	EDGESEG* readpath = new EDGESEG[MIN_LEG];
	//EDGESEG* rev_path = (EDGESEG *) smallapp::ckalloc(MIN_LEG * sizeof(EDGESEG *));
	EDGESEG* rev_path = new EDGESEG[MIN_LEG];
	seq = new char[MIN_LEG];

	int r, i, num_reads1, n1, np1;
	char	*readfile;
	for(i = 0; i <= 10; i ++)	{
		pathstats[i] = 0;
	}
	n1 = np1 = 0;
	int totreadfile = 1;
	if(read2file[0] != 0) totreadfile = 2;
	for(r = 0; r < totreadfile; r ++) {
		if(r == 0) readfile = read1file;
		else readfile = read2file;
		num_reads1 = reads::countreads(readfile, fastq);
		cout<<"processing 1st read file: "<<readfile<<endl;
		cout<<"total reads expected: "<<num_reads1<<endl;
        	std::ifstream fp_reads(readfile);
	        if(!fp_reads) {
       	         	cout<<"cannot open file "<<readfile<<endl;
                	exit(0);
        	}

		//Binary search of k-mers in reads against K-mers in the graph (with either indegree > 1 or outdegree > 1
		tot_reads_mapped = 0;
		if(!fastq) {
			int k = reads::get1readname(fp_reads, seqname);	// more sequence, k = 1; otherwise, k = 0 
			len_seq = reads::get1readfasta(fp_reads, seq, seqname_new);
		}
		else len_seq = reads::get1readfastq(fp_reads, seq, seqname);
		while(len_seq > 0)	{
			len_path = map1read(seq, len_seq, readpath);
			if(len_path > 0)	{
				if(len_path < 10)	{
					pathstats[len_path] ++;
				} else	{
					pathstats[10] ++;
				}
				tot_reads_mapped ++;
				readpath -> length = len_path; //#of edges a mapping spans
				path::reverse_path(readpath, rev_path, len_path, edge);
				if(r == 0) pathlist[rev_path[0].index] = path::insert_path(rev_path, len_path, pathlist[rev_path[0].index]);
				else pathlist[readpath[0].index] = path::insert_path(readpath, len_path, pathlist[readpath[0].index]);
				//For 2nd read, the actual path should be the one on the rp, when strand-specific protocol is used
				np1 ++;	
				if(stranded == 0)	{
					if(r == 0) pathlist[readpath[0].index] = path::insert_path(readpath, len_path, pathlist[readpath[0].index]);
					else pathlist[rev_path[0].index] = path::insert_path(rev_path, len_path, pathlist[rev_path[0].index]);
					//For 2nd read, the actual path should be the one on the rp, when strand-specific protocol is used
					np1 ++;	
				}
			}
			if(!fastq) {
				strcpy(seqname, seqname_new);
				len_seq = reads::get1readfasta(fp_reads, seq, seqname_new);
			}
			else len_seq = reads::get1readfastq(fp_reads, seq, seqname);
			n1 ++;
		}
		fp_reads.close();
		cout<<"reads file processed: "<<np1<<" paths identified from "<<n1<<" reads."<<endl;
	}

	delete[] readpath;
	delete[] rev_path;
	delete[] seq;

	printf("------------------------------------\n");
	printf("Path length	Number of paths\n");
	printf("------------------------------------\n");
	for(i = 1; i < 10; i ++)	{
		printf("%d	%d\n", i, pathstats[i]);
	}
	printf(">10	%d\n", pathstats[10]);
	printf("Total	%d (%.2f%%)\n", tot_reads_mapped, (double) tot_reads_mapped * 100 / n1);
	printf("------------------------------------\n");
}

int seqgraph::map1read(char *seq, int len_seq, EDGESEG *path)
{
	int	i, j, k, l, m, n;
	int	mut, tot_mut;
	int	hashv, hash0;
	__int128_t	index;
	int	shift;
	int	len_path, len_tmppath;
	int	MAX_LEG = 100;
	EDGESEG	tmppath[MAX_LEG];
	int	matchpos[MAX_LEG], matchvertex[MAX_LEG];

	shift = (kmersize - hashw) * 2;

	index = 0;
	for(i = 0; i < kmersize; i ++)	{
		index = (index << 2) + seq[i];
		if(ifprint) printf("seq %d %d\n", i, seq[i]);
	}
	if(ifprint) printf(" kmersize %d, shift %d, index %ld maskv %ld\n", kmersize, shift, index, maskv);

	k = 0;
	while(i < len_seq)	{
		hashv = index >> shift;
		for(j = vindex[hashv] - 1; j < vindex[hashv + 1] - 1; j ++)	{
			if(vertex[j].index == index)	{
				matchpos[k] = i;
				matchvertex[k] = j;
				k ++;
			}
		}
		i ++;
		if(i < len_seq)	{
			index = (maskv & (index << 2)) + seq[i];	/* shift two bits; mask out bits before k	*/
		}
	}
	if(ifprint) printf(" index %ld\n", index);
	if(ifprint) printf(" find seeds %d first matchpos %d matchvertex %d\n", k, matchpos[0], matchvertex[k]);

	if(k == 0)	{
		return(0);
	}

	tot_mut = 0;
	len_tmppath = backalign(seq, matchpos[0], &vertex[matchvertex[0]], tmppath, 0, &mut);
	if(ifprint) printf(" backalign, length = %d\n", len_tmppath);
	if(len_tmppath < 0)	{
		return(0);
	}
	tot_mut += mut;
	for(j = 0; j < len_tmppath; j ++)	{
		path[j] = tmppath[j];
	}
	len_path = len_tmppath;
	for(i = 0; i < k - 1; i ++)	{
		len_tmppath = intalign(&seq[matchpos[i]], matchpos[i + 1] - matchpos[i], &vertex[matchvertex[i]], &vertex[matchvertex[i + 1]], tmppath, 0, &mut);
		if(len_tmppath < 0)	{
			return(0);
		}
		tot_mut += mut;
		for(j = 0; j < len_tmppath; j ++)	{
			path[len_path ++] = tmppath[j];
		}
	}
	len_tmppath = forwardalign(&seq[matchpos[k - 1]], len_seq - matchpos[k - 1], &vertex[matchvertex[k - 1]], tmppath, 0, &mut);
	if(ifprint) printf(" forwardalign, length = %d\n", len_tmppath);
	if(len_tmppath < 0)	{
		return(0);
	}
	tot_mut += mut;
	for(j = 0; j < len_tmppath; j ++)	{
		path[len_path ++] = tmppath[j];
	}
	if(ifprint) printf("ngap %d, maxmut %d, tot_mut %d\n", ngap, maxmut, tot_mut);
	//if(tot_mut <= ngap)	{
	if(tot_mut <= maxmut)	{
		return(len_path);
	} else	{
		return(0);
	}
}

int seqgraph::intalign(char *seq, int length, VERTEX *vertex1, VERTEX *vertex2, EDGESEG *tmppath, int len_tmppath, int *rmut)
{
	int	i, j;
	int	loc, mut, pos, min_mut, min_loc;
	EDGE	*thisedge, *min_edge;
	EDGELIST	*edgelist;

	if(vertex1 == vertex2)	{
		if(length <= ngap)	{
			*rmut = 0;
			return(0);
		}
	} else if(vertex1 -> outdegree == 0)	{
		return(-1);
	}

	min_mut = ngap + 1;
	edgelist = vertex1 -> nextedge;
	while(edgelist)	{
		thisedge = edgelist -> edge;
		if(abs(thisedge -> length - length - kmersize) <= ngap && &vertex[thisedge -> nodeindex[1]] == vertex2)	{
			loc = for_align(seq, length, &(thisedge -> seq[kmersize]), thisedge -> length - kmersize, &mut, &pos);
			if(mut < min_mut)	{
				min_mut = mut;
				min_edge = thisedge;
				min_loc = loc;
			}
		}
		edgelist = edgelist -> next;
	}
	if(min_mut <= ngap)	{
		*rmut = min_mut;
		loc = min_loc;
		tmppath[len_tmppath].index = min_edge -> index; 
		tmppath[len_tmppath].num_reads = 1;
		tmppath[len_tmppath].total_length = min_edge -> length;
		tmppath[len_tmppath].pos[0] = 1;
		tmppath[len_tmppath].pos[1] = min_edge -> length;
		len_tmppath ++;
		return(len_tmppath);
	} else	{
		return(-1);
	}
}

int seqgraph::forwardalign(char *seq, int length, VERTEX *vertex0, EDGESEG *tmppath, int len_tmppath,int *rmut)
{
	int	i, j;
	int	loc, mut, pos, min_mut, min_loc, min_pos;
	EDGE	*thisedge, *min_edge;
	EDGELIST	*edgelist;

	if(length <= ngap)	{
		*rmut = 0;
		return(0);
	} else if(vertex0 -> outdegree == 0)	{
		return(-1);
	}
	min_mut = ngap + 1;
	edgelist = vertex0 -> nextedge;
	while(edgelist)	{
		thisedge = edgelist -> edge;
		if(thisedge -> nodeindex[0] != thisedge -> nodeindex[1] && thisedge -> length >= length - ngap)	{
			loc = for_align(seq, length, &(thisedge -> seq[kmersize]), thisedge -> length - kmersize, &mut, &pos);
			if(mut < min_mut)	{
				min_mut = mut;
				min_edge = thisedge;
				min_loc = loc;
				min_pos = pos;
			}
		}
		edgelist = edgelist -> next;
	}
	if(min_mut <= ngap)	{
		*rmut = min_mut;
		tmppath[len_tmppath].index = min_edge -> index; 
		tmppath[len_tmppath].num_reads = 1;
		tmppath[len_tmppath].total_length = min_pos + kmersize;
		tmppath[len_tmppath].pos[0] = 1;
		tmppath[len_tmppath].pos[1] = min(min_edge -> length, min_pos + kmersize + 1);
		len_tmppath ++;
		return(len_tmppath);
	} else	{
		return(-1);
	}
}

int seqgraph::for_align(char *seq1, int length1, char *seq2, int length2, int *mut, int *pos)
{
	int	i, j, k, l, diff;
	int	*score, s, sl, min_score;
	int	loc;

	diff = length2 - length1;
	score = (int *) smallapp::ckalloc((2 * ngap + 1) * sizeof(int));
	for(j = 0; j <= 2 * ngap; j ++)	{
		score[j] = abs(j - ngap);
	}
	for(i = 0, k = 0; i < length1 && k < length2; i ++, k ++)	{
		j = max(-i, -ngap);
		s = score[ngap + j];
		sl = ngap;
		while(j <= min(length2 - i - 1, ngap))	{
			if(seq1[i] == seq2[i + j])	{
				score[ngap + j] = min(score[ngap + j], s + 1);
			} else	{
				score[ngap + j] = min(score[ngap + j], s) + 1;
			}
			score[ngap + j] = min(score[ngap + j], sl + 1);
			sl = score[ngap + j];
			j ++;
		}
	}
	min_score = ngap * 10 + 1;
	l = 0;
	for(j = 0; j <= 2 * ngap; j ++)	{
		if(score[j] < min_score)	{
			min_score = score[j];
			l = j - ngap;
		}
	}
	*mut = min_score;
	if(i == length1)	{
		loc = i - 1;
		*pos = length1 + l - 1;
	} else if(k == length2) {
		loc = length2 - l - 1;
		*pos = k - 1;
	}

	free((void *) score);
	return(loc);
}

int seqgraph::backalign(char *seq, int length, VERTEX *vertex0, EDGESEG *tmppath, int len_tmppath, int *rmut)
{
	int	i, j;
	int	loc, mut, pos, min_mut, min_loc, min_pos;
	EDGE	*thisedge, *min_edge;
	EDGELIST	*edgelist;

	if(length - kmersize < ngap)	{
		*rmut = 0;
		return(0);
	} else if(vertex0 -> indegree == 0)	{
		return(-1);
	}
	min_mut = ngap + 1;
	edgelist = vertex0 -> lastedge;
	while(edgelist)	{
		thisedge = edgelist -> edge;
		if(thisedge -> nodeindex[1] != thisedge -> nodeindex[0] && thisedge -> length >= length - ngap)	{
			loc = rev_align(seq, length - kmersize, thisedge -> seq, thisedge -> length - kmersize, &mut, &pos);
			if(mut < min_mut)	{
				min_mut = mut;
				min_edge = thisedge;
				min_loc = loc;
				min_pos = pos;
			}
		}
		edgelist = edgelist -> next;
	}
/*
printf("min_mut %d\n", min_mut);
*/
	if(min_mut <= ngap)	{
		*rmut = min_mut;
		tmppath[len_tmppath].index = min_edge -> index; 
		tmppath[len_tmppath].num_reads = 1;
		tmppath[len_tmppath].total_length = min_edge -> length - pos + 1;
		tmppath[len_tmppath].pos[0] = max(1, min_pos);
		tmppath[len_tmppath].pos[1] = min_edge -> length;
		len_tmppath ++;
		return(len_tmppath);
	} else	{
		return(-1);
	}
}

int seqgraph::rev_align(char *seq1, int length1, char *seq2, int length2, int *mut, int *pos)
{
	int	i, j, k, l, diff;
	int	*score, s, sl, min_score;
	int	loc;

	diff = length2 - length1;
	score = (int *) smallapp::ckalloc((2 * ngap + 1) * sizeof(int));
	for(j = 0; j <= 2 * ngap; j ++)	{
		score[j] = abs(j - ngap);
	}
	for(i = length1 - 1, k = length2 - 1; i >= 0 && k >= 0; i --, k --)	{
		j = max(-k, -ngap);
		s = score[ngap + j];
		sl = ngap;
		while(j <= min(length2 - 1 - k, ngap))	{
			if(seq1[i] == seq2[k + j])	{
				score[ngap + j] = min(score[ngap + j], s + 1);
			} else	{
				score[ngap + j] = min(score[ngap + j], s) + 1;
			}
			score[ngap + j] = min(score[ngap + j], sl + 1);
			sl = score[ngap + j];
			j ++;
		}
	}
	min_score = ngap * 10 + 1;
	l = 0;
	for(j = 0; j <= 2 * ngap; j ++)	{
		if(score[j] < min_score)	{
			min_score = score[j];
			l = j - ngap;
		}
	}
	*mut = min_score;
	if(i == -1)	{
		loc = i + 1;
		*pos = diff - l;
	} else if(k == -1) {
		loc = -diff + l;
		*pos = k + 1;
	}

	free((void *) score);
	return(loc);
}
