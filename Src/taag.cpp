#include "taag.h"
#include "path.h"
#include "smallapp.h"
#include "soap.h"

void taag::assembly(void)
{
	int	i, j, k, l, n, nc, n0, n1;
	int	num_part;
	int	index, edgeindex;
	int	pos[2];
	int	eidx;
	int	segcurr, segbeg;
	char	*label;
	EDGELIST	*edgelist;
	EDGESEG	*pathseg, *newpathlist;
	EDGESEG	*pathbeg, *pathcurr, *pathlast;

	cout<<"\n>>>Now assemble transcripts..."<<endl;
	//i = 610685;
	printf("Check at the beginning of assembly\n");
        //printf("Edge i = %d tag %d, begtag %d endtag %d length %d pathseg-length %d\n", i, edge[i].tag, edge[i].begtag, edge[i].endtag, edge[i].length, pathlist[i] -> length);

	for(i = 0; i < num_edge; i ++)	edge[i].tag = 0; //inilize the tag -- set to 0 

	for(i = 0; i < num_edge; i ++)	{
		pathseg = pathlist[i];
		//cout<<"check i "<<i<<endl;
		while(pathseg)	{
			//cout<<" pathseg -> length "<<pathseg->length<<endl;
			for(j = 0; j < pathseg -> length; j ++)	{
				eidx = pathseg[j].index;
				edge[eidx].tag = 3;	
				//3: edges covered by reads;
				//0: not in graph
				//1: traversed in graph
				//2: complement edge traversed in graph
				edge[eidx].num_reads ++;
				edge[eidx].total_length += (pathseg[j].pos[1] - pathseg[j].pos[0] + 1);
			}
			if(pathseg -> length > 1)	{ //why add this??
				eidx = pathseg[pathseg -> length - 1].index;
				newpathlist = smallapp::newEDGESEGpt();
				newpathlist -> index = eidx;
				newpathlist -> pos[0] = pathseg[pathseg -> length - 1].pos[0];
				newpathlist -> pos[1] = pathseg[pathseg -> length - 1].pos[1];
				pathlist[eidx] = path::insert_single_path(pathlist[eidx], newpathlist); 
			}
			pathseg = pathseg -> next;
		}
	}

	//i = 610685;
	//printf("Check 222\n");
        //printf("Edge i = %d tag %d, begtag %d endtag %d length %d pathseg-length %d\n", i, edge[i].tag, edge[i].begtag, edge[i].endtag, edge[i].length, pathlist[i] -> length);

	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].tag == 0)	continue;

		//create pathindex (pathlist[i]) & then qsort pathindex
		pathseg = pathlist[i];
		n = get_path_num(pathseg);
		if(n == 0)	continue;
		label = new char[n];
		for(j = 0; j < n; j ++) label[j] = 0;
		EDGESEG **pathindex = new EDGESEG*[n];
		n = get_path_list(pathindex, pathseg);

		//qsort to check , YY
		std::sort(pathindex, pathindex + n, sorter());

		k = 0;
		n0 = 1;
		pathcurr = pathindex[0];
		segbeg = segcurr = 0;

		pathlast = (EDGESEG *) NULL;
		for(j = 1; j < n; j ++)	{
			if(pathcurr -> pos[1] - pathindex[j] -> pos[0] >= overlaplen)	{
				//pathcurr & pathindex[j] overlaps
				if(pathindex[j] -> length == 1)	{
					//delete pathindex[j]; update pathcurr info (pos[1], num_reads, total_length)
					if(pathindex[j] -> pos[1] > pathcurr -> pos[1])	{
						pathcurr -> pos[1] = pathindex[j] -> pos[1];
					}
					pathcurr -> num_reads += pathindex[j] -> num_reads;
					pathcurr -> total_length += pathindex[j] -> total_length;
					label[j] = 1;	//mark segment j deleted
					smallapp::delEDGESEGpt(pathindex[j]);
				} else	{
					//delete pathcurr and update pathindex[j] (pos[0], num_reads, and total_length)
					if(pathindex[j] -> pos[0] > pathcurr -> pos[0])	{
						pathindex[j] -> pos[0] = pathcurr -> pos[0];
					}
					pathindex[j] -> num_reads += pathcurr -> num_reads;
					pathindex[j] -> total_length += pathcurr -> total_length;
					if(segbeg == segcurr)	segbeg = j;
					smallapp::delEDGESEGpt(pathcurr);
					label[segcurr] = 1;	// mark segment segcurr deleted
					segcurr = j;
					pathcurr = pathindex[j];
					if(pathlast)	pathlast -> next = pathcurr;
					pathlast = pathcurr;
				}
			} else	{
				//no overlap, re-start the checking of overlaps from pathindex[j] (set as pathcurr)
				pathlast = pathcurr;
				pathcurr = pathindex[j];
				segcurr = j;
				pathlast -> next = pathcurr;
				n0 ++;
			}
		}
		//cout<<"Edge "<<i<<" after simple merge, n0 "<<n0<<" original: "<<n<<endl;

		//chain together all mapped segments (& not deleted) for each edge i (by setting next)
		//so these segments can be explored from pathlist[i] following next
		pathlist[i] = pathindex[segbeg]; 
		pathbeg = pathlist[i]; //pathbeg is temp
		for(j = segbeg + 1; j < n; j ++)	{
			if(!label[j])	{ //pathindex[j] remained 
				pathbeg -> next = pathindex[j];
				pathbeg = pathindex[j];
			}
		}
		pathbeg -> next = (EDGESEG *) NULL;
		delete[] pathindex; //delete the pointer holder
		delete[] label;

		n1 = get_path_num(pathlist[i]);
		if(n1 != n0)	{
			printf("Path list not match: %d %d.\n", n1, n0);
			exit(-1);
		}

		//update edge (data structure for edges/contigs) based on pathlist
		pathbeg = pathlist[i];
		n = 0;
		while(pathbeg)	{
			//printf("Edge %d, path %d, length %d\n", i, n, pathbeg -> length);
			if(pathbeg -> length > 1) { //a mapping that spans multiple edges/contigs
				eidx = pathbeg[pathbeg -> length - 1].index;
				//update firstseg: first mapped segment for edge edgeedge
				edge[eidx].firstseg = pathbeg[pathbeg -> length - 1].pos[1];
				edge[eidx].firstseg_num_reads += pathbeg[pathbeg -> length - 1].num_reads;
				edge[eidx].firstseg_total_length += pathbeg[pathbeg -> length - 1].total_length;
				//update lastseg: last mapped segment for edge i
				edge[i].lastseg = pathbeg -> pos[0]; //last mapped segment for edge i
				edge[i].lastseg_num_reads += pathbeg -> num_reads;
				edge[i].lastseg_total_length += pathbeg -> total_length;
			}
			pathbeg = pathbeg -> next;
			n ++;
		}
	}
	cout<<" Merge done"<<endl;

	//i = 22266;
	//i = 610685;
	//printf("Check after merge\n");
        //printf("Edge i = %d tag %d, begtag %d endtag %d length %d pathseg-length %d\n", i, edge[i].tag, edge[i].begtag, edge[i].endtag, edge[i].length, pathlist[i] -> length);

	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].tag == 0)	continue;

		pathbeg = pathlist[i];
		while(pathbeg)	{
			if(pathbeg -> length > 1)	{
				eidx = pathbeg[pathbeg -> length - 1].index;
				pathseg = pathlist[eidx];
				if(pathseg && pathbeg[pathbeg -> length - 1].pos[1] - pathseg -> pos[0] >= overlaplen &&
				   pathseg -> pos[1] > pathbeg[pathbeg -> length - 1].pos[1])	{
					pathbeg[pathbeg -> length - 1].pos[1] = pathseg -> pos[1]; 
					// extend the end position of the last segment in a path to the end position 
					// of the first segment in the edge where the segment is from if these two segments overlap
					edge[eidx].firstseg = pathseg -> pos[1];
					edge[eidx].firstseg_num_reads += pathseg -> num_reads;
					edge[eidx].firstseg_total_length += pathseg -> total_length;
				} else	{
					edge[eidx].firstseg = pathbeg[pathbeg -> length - 1].pos[1];
				}
			}
			pathbeg = pathbeg -> next;
		}
	}
	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].tag == 0)	continue;

		pathbeg = pathlist[i];
		if(!pathbeg)	continue;
		if(pathbeg -> pos[0] == 1)	{
			if(pathbeg -> pos[1] > edge[i].firstseg)	{
				edge[i].firstseg = pathbeg -> pos[1];
			}
		}
	}
	cout<<" Extension done."<<endl;

	for(i = 0; i < num_edge; i ++)	{
		pathseg = pathlist[i];
		while(pathseg)	{
			if(pathseg -> length > 1)	{
				if(edge[i].begtag == 0)	{
					edge[i].begtag = 1; // tag the edges always at the beginning of a path for starting the search
				}
				for(j = 1; j < pathseg -> length - 1; j ++)	{
					eidx = pathseg[j].index;
					edge[eidx].begtag = edge[eidx].endtag = 3; // edges in the middle of a path
				}
				eidx = pathseg[j].index;
				if(edge[eidx].endtag == 0)	{
					edge[eidx].endtag = 1; // tag the edges always at the end of some paths
				}
			}
			pathseg = pathseg -> next;
		}
	}
}

void taag::writefile(char *outfile, bool sepmult, bool savegraph, int minlen)
{
	//import initialization
	transcriptlen_min = minlen;
	multtranscriptlen = 0;
	multtranscriptlen_long = 0;
	multtranscriptlen_max = 0;
	multtranscript = multtranscript_long = 0;
	multtranscriptlen = multtranscriptlen_long = 0;

	EDGELIST	*edgelist;
	EDGESEG	*pathseg, *newpathlist;
	EDGESEG	*pathbeg, *pathcurr, *pathlast;
	EDGESEG **pathindex;

	char	transcriptfile[1000];
	sprintf(transcriptfile, "%s.fa", outfile);
	std::ofstream fs(transcriptfile);
	if(!fs) { cout<<"cannot open file to write "<<transcriptfile<<endl; exit(0); }

	cout<<"\n>>>Output transcripts to file(s)"<<"..."<<endl;
	cout<<"Only transcripts of >= "<<transcriptlen_min<<" bp are considered"<<endl;

	//Output single-edge transcripts
	cout<<"\n>>>Output single-edge transcripts to "<<outfile<<" .."<<endl;
	int	i, l;
	int 	index = 0;
	int	index_long = 0;
	double 	avecov;
	long	totalbp = 0;
	long 	totalbp_long = 0;
	int	maxlen = 0;
	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].tag != 3)	continue;
		pathbeg = pathlist[i];
		while(pathbeg)	{
			if(pathbeg -> length == 1)	{
				if((edge[i].begtag == 0 || (edge[i].begtag == 1 && pathbeg -> pos[1] <= edge[i].length)) &&
				   (edge[i].endtag == 0 || (edge[i].endtag == 1 && pathbeg -> pos[0] > 1)))	{
					l = pathbeg -> pos[1] - pathbeg -> pos[0] + 1;
					avecov = ((double) pathbeg -> total_length) / l;
					index ++;
					totalbp += l;
					if(l >= transcriptlen_min) {
						output_transcript(fs, edge[i].seq, index_long, pathbeg -> pos[0] - 1, pathbeg -> pos[1] - 1, pathbeg -> num_reads, avecov);
						totalbp_long += l;
						index_long ++;
					}
					if(l > maxlen) maxlen = l;
				}
			}
			pathbeg = pathbeg -> next;
		}
		if(edge[i].begtag == 0 && edge[i].endtag == 0)	{
			edge[i].tag = 0;
			if(stranded == 0 && edge[i].bal_edge != i)	{
				edge[edge[i].bal_edge].tag = 0;
			}
		}
	}
	int single = index_long;
	long singlebp = totalbp_long;
	//cout<<" Transcripts from single edges: "<<index<<" total-bp "<<totalbp<<" average-length "<<int(totalbp/index)<<" maxlen "<<maxlen<<endl;
	cout<<" Transcripts from single edges: "<<index_long<<" total-bp "<<totalbp_long<<" average-length "<<int(totalbp_long/index_long)<<" maxlen "<<maxlen<<endl;

	int num_part = 0;
	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].tag == 3)	{
			num_part ++;
		}
	}
	cout<<" Now traverse edges to recover transcripts spanning edges"<<endl;
	cout<<"  Edges to be traversed: "<<num_part<<endl;

	int *SourceEdgeIndex = new int[num_edge];

	int nc = 0;
	for(i = 0; i < num_edge; i ++)	{
		if(edge[i].begtag == 1)	{
			SourceEdgeIndex[nc ++] = i;
		}
	}
	cout<<"  Source edges: "<<nc<<endl;

	//Output transcript subgraphs
	std::ofstream *fm; 
	if(sepmult) {
		char	multedgefile[1000];
		sprintf(multedgefile, "%s.mult", outfile);
		fm = new std::ofstream(multedgefile);
		if(!fm) { cout<<"cannot open file to write "<<multedgefile<<endl; exit(0); }
	}
	else fm = &fs;


	int num_updatededge = 0;
	EDGE* updatededge = new EDGE[num_edge];
	for(i = 0; i < nc; i ++)	{ //for each source segment
		index = SourceEdgeIndex[i];
//printf("source edge %d, %d, tag %d n %d\n", i, index, edge[index].tag, n);
//getchar();
		if(edge[index].tag != 0)	{
			if(edge[index].endtag == 1 && edge[index].firstseg > 0 && edge[index].firstseg < edge[index].length && edge[index].parttag == 0)	{
				avecov = edge[index].firstseg_total_length / edge[index].firstseg;
				outputmultitranscript(*fm, edge[index].seq, edge[index].firstseg, edge[index].firstseg_num_reads, avecov);
				edge[index].parttag = 1;
			}
			if(edge[index].tag == 3)	{
				num_updatededge = traversegraph(&edge[index], updatededge, num_updatededge, *fm, 0);
			}
		}
	}
	int avelen = multtranscriptlen / multtranscript;
	int alen_long = multtranscriptlen_long / multtranscript_long; 
	//cout<<" Transcripts spanning edges: "<<multtranscript<<" total-bp "<<multtranscriptlen<<" average-length "<<avelen<<" maxlen "<<multtranscriptlen_max<<endl;
	cout<<" Transcripts spanning edges: "<<multtranscript_long<<" total-bp "<<multtranscriptlen_long<<" average-length "<<alen_long<<" maxlen "<<multtranscriptlen_max<<endl;

	delete[] SourceEdgeIndex;

	//output transcripts from graph (those cannot be resolved)
	cout<<">>>Output partial transcripts to "<<transcriptfile<<" .."<<endl;

	if(savegraph) {
		char	edgefile[1000];
		sprintf(edgefile, "%s.edge", outfile);
		std::ofstream fe(edgefile);
		if(!fe) { cout<<"cannot open file to write "<<edgefile<<endl; exit(0); }
		fe<<"EDGEs "<<num_updatededge<<endl;
		maxlen = 0;
		for(i = 0; i < num_updatededge; i ++)	{
			vertex[updatededge[i].nodeindex[0]].tag = vertex[updatededge[i].nodeindex[1]].tag = 1;
			output_updatededge(fe, updatededge[i]);
			if(updatededge[i].length > maxlen) maxlen = updatededge[i].length;
		}
		fe.close();
	}

	totalbp = totalbp_long = 0;
	long totalbp_m = 0;
	long totalbp_m_long = 0;
	int n = 0;
	int n_long = 0;
	for(i = 0; i < num_updatededge; i ++) {
		if(updatededge[i].tag == 1)	{
			if(!updatededge[i].multitag) {
				if(updatededge[i].length == 0) avecov = 0; //????
				else avecov = ((double) updatededge[i].firstseg_total_length) / updatededge[i].length;
				n ++;
				totalbp += updatededge[i].length;
				if(updatededge[i].length >= transcriptlen_min) { 
					outputpartialtranscript(fs, updatededge[i].seq, n, updatededge[i].length, updatededge[i].firstseg_num_reads, avecov, updatededge[i].multitag);
					n_long ++;
					totalbp_long += updatededge[i].length;
				}
			}
		}
	}
	int m = 0;
	int m_long = 0;
	for(i = 0; i < num_updatededge; i ++)	{
		if(updatededge[i].tag == 1)	{
			if(updatededge[i].multitag) {
				if(updatededge[i].length == 0) avecov = 0; //????
				else avecov = ((double) updatededge[i].firstseg_total_length) / updatededge[i].length;
				m ++;
				totalbp_m += updatededge[i].length;
				if(updatededge[i].length >= transcriptlen_min) {
					outputpartialtranscript(*fm, updatededge[i].seq, m, updatededge[i].length, updatededge[i].firstseg_num_reads, avecov, updatededge[i].multitag);
					m_long ++;
					totalbp_m_long += updatededge[i].length;
				}
			}
			delete[] updatededge[i].seq;
		}
	}
	fs.close();
	if(sepmult) delete fm;

	int graph_s = n_long;
	int graph_m = m_long;
	int graph = graph_s + graph_m;
	long graph_s_bp = totalbp_long;
	long graph_m_bp = totalbp_m_long;
	long graph_bp = graph_s_bp + graph_m_bp;
	cout<<"Transcripts extracted from transcript graph "<<graph<<" total-bp "<<graph_bp<<" average-length "<<graph_bp/graph<<" maxlen "<<maxlen<<endl;
	cout<<"    single edges "<<graph_s<<" total-bp "<<graph_s_bp<<endl;
	cout<<"    spanning multiple edges "<<graph_m<<" total-bp "<<graph_m_bp<<endl;

	//cout<<" "<<num_updatededge<<" updated edges"<<endl;

	if(savegraph) {
		char	vertexfile[1000];
		sprintf(vertexfile, "%s.vertex", outfile);
		std::ofstream fv(vertexfile);
		if(!fv) { cout<<"cannot open file to write "<<vertexfile<<endl; exit(0); }
		n = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i].tag)	{
				output_vertex(fv, vertex[i]);
				n ++;
				if(n % 8 == 0)	{
					fv<<endl;
				}
			}
		}
		fv.close();
	}
	
	//k  single single-bp multi multi-bp graph graph-bp graph-s graph-s-bp graph-m graph-m-bp
	cout<<"TAG "<<kmersize<<" "<<single<<" "<<singlebp<<" "<<multtranscript_long<<" "<<multtranscriptlen_long<<" ";
	cout<<graph<<" "<<graph_bp<<" "<<graph_s<<" "<<graph_s_bp<<" "<<graph_m<<" "<<graph_m_bp<<endl;

	delete[] updatededge;
}

int taag::get_path_num(EDGESEG *pathlist_t)
{
	int	n = 0;
	EDGESEG	*pathlist0 = pathlist_t;

	while(pathlist0)	{
		n ++;
		pathlist0 = pathlist0 -> next;
	}
	return(n);
}

//int taag::get_path_list(vector<EDGESEG*> &pathindex, EDGESEG *pathlist_t)
int taag::get_path_list(EDGESEG** pathindex, EDGESEG *pathlist_t)
{
	EDGESEG	*pathlist0 = pathlist_t;

	int n = 0;
	while(pathlist0)	{
		//pathindex.insert(pathindex.end(), pathlist0);
		pathindex[n] = pathlist0;
		pathlist0 = pathlist0 -> next;
		n ++;
	}
	//return(pathindex.size());
	return n;
}

bool taag::comparseg(const EDGESEG* s1, const EDGESEG* s2)
{
	return (s1 -> pos[0] > s2 -> pos[0]);
}

void taag::output_vertex(ofstream &fp, VERTEX vertex)
{
	char	str[100];

	soap::index2code(vertex.index, str);
	fp<<str;
}

void taag::output_updatededge(ofstream &fp, EDGE updatededge)
{
	int	dir, unk;
	char	str1[100], str2[100], str[1000];

	soap::index2code(vertex[updatededge.nodeindex[0]].index, str1);
	soap::index2code(vertex[updatededge.nodeindex[1]].index, str2);
	if(updatededge.tag == 1)	{
		dir = 0;
	} else if(updatededge.tag == 2)	{
		dir = 1;
	}
	unk = 0;
	sprintf(str, ">length %d,%s,%s,%d,%d", updatededge.length - kmersize, str1, str2, dir, unk);
	fp<<str<<endl;
}

void taag::output_transcript(ofstream &fp, char *seq, int index, int pos1, int pos2, int num_reads, double avecov)
{
	int	i, j, k, l, m, n;
	int	length;
	char	str[1000];	

	length = pos2 - pos1 + 1;
	sprintf(str, ">transcript%d %d %d %.2f", index, length, num_reads, avecov);
	fp<<str<<endl;
	for(i = pos1, j = 0; i <= pos2; i ++, j++)	{
		sprintf(str, "%c", na_name[seq[i]] - 'a' + 'A');
		fp<<str;
		if(j % 60 == 59)	{
			fp<<endl;
		}
	}
	if(j % 60 != 0)	{
		fp<<endl;	
	}
}

void taag::outputmultitranscript(ofstream &fp, char *seq, int length, int num_reads, double avecov)
{
	int	i, j, k, l, m, n;
	char	str[1000];

	multtranscript ++;
	multtranscriptlen += length;
	if(multtranscriptlen_max < length) multtranscriptlen_max = length;

	if(length < transcriptlen_min) return;

	multtranscript_long ++;
	multtranscriptlen_long += length;
	sprintf(str, ">multiedgetranscript%d %d %d %.2f", multtranscript_long, length, num_reads, avecov);
	fp<<str<<endl;
	for(i = 0; i < length; i ++)	{
		sprintf(str, "%c", na_name[seq[i]] - 'a' + 'A');
		fp<<str;	
		if(i % 60 == 59)	{
			fp<<endl;
		}
	}
	if(i % 60 != 0)	{
		fp<<endl;
	}
}

void taag::outputpartialtranscript(ofstream &fp, char *seq, int index, int length, int num_reads, double avecov, bool multitag)
{
	int	i, j, k, l, m, n;
	char	str[1000];
	if(multitag) 
		sprintf(str, ">partialmultitranscript%d %d %d %.2f multi %d", index, length, num_reads, avecov, multitag);
	else
		sprintf(str, ">partialtranscript%d %d %d %.2f multi %d", index, length, num_reads, avecov, multitag);
	fp<<str<<endl;
	for(i = 0; i < length; i ++)	{
		sprintf(str, "%c", na_name[seq[i]] - 'a' + 'A');
		fp<<str;
		if(i % 60 == 59)	{
			fp<<endl;
		}
	}
	if(i % 60 != 0)	{
		fp<<endl;
	}
}

int taag::traversegraph(EDGE *thisedge, EDGE *updatededge, int num_updatededge, ofstream &fp, char tag)
{
	int	i, j, k, l, m, n;
	int	bal_edge;
	int	num_reads, total_length;
	int	endnodeindex = -1;
	int	outdegree = 0;
	int	*nextedgeindex;
	double	avecov;
	EDGELIST	*edgelist;
	int	length;
	char	*seq;

	//seq = (char *) ckalloc(MIN_CON_LEN * sizeof(char));
	seq = new char[MIN_CON_LEN];

	updatededge[num_updatededge].nodeindex[0] = thisedge -> nodeindex[0];
	length = 0;
	num_reads = total_length = 0;

	length = extendseq(seq, length, thisedge, &endnodeindex, &num_reads, &total_length, 0);
	if(endnodeindex >= 0)	{
		nextedgeindex = new int[vertex[endnodeindex].outdegree + 10];
		outdegree = chk_nextedge(&vertex[endnodeindex], nextedgeindex);
	}

	if(outdegree == 0 && tag == 0)	{//when tag = 1, the edge should be part of the transcript graph	*/
		if(thisedge -> tag != 2 && length > 0)	{
			// output only the transcripts whose reverse complements are not traversed (rev_edge -> tag != 2)
			avecov = ((double) total_length) / length;
			outputmultitranscript(fp, seq, length, num_reads, avecov);
		}
	} else if(length > overlaplen)	{
		if(!(length == thisedge -> length && thisedge -> parttag == 1 && thisedge -> firstseg >= thisedge -> length))	{
			updatededge[num_updatededge].nodeindex[1] = endnodeindex;
			if(endnodeindex != thisedge -> nodeindex[1])	{
				updatededge[num_updatededge].multitag = true; // this updated edge spans multiple edges in the de Bruijn graph
			} else	{
				updatededge[num_updatededge].multitag = false;
			}
			updatededge[num_updatededge].seq = new char[length];
			for(i = 0; i < length; i ++)	{
				updatededge[num_updatededge].seq[i] = seq[i];
			}
			updatededge[num_updatededge].length = length;
			updatededge[num_updatededge].firstseg_num_reads = num_reads;
			updatededge[num_updatededge].firstseg_total_length = total_length;
			if(thisedge -> tag == 2 && stranded == 0)	{
				updatededge[num_updatededge].tag = 2;//edge -> tag == 2 => revserse complement of the edge is traversed
			} else	{
				updatededge[num_updatededge].tag = 1;
			}
			num_updatededge ++;
		}
	}
	if(endnodeindex >= 0)	{
		delete[] nextedgeindex;
	}
	delete[] seq;

	edgelist = vertex[endnodeindex].nextedge;
	while(edgelist)	{
		if(edgelist -> edge -> parttag == 0 && edgelist -> edge -> tag != 0)	{
			length = edgelist -> edge -> firstseg;
			if(!(edgelist -> edge -> lastseg == 1 && edgelist -> edge -> begtag == 1) && edgelist -> edge -> lastseg != 1 && length > overlaplen)	{
				avecov = ((double) edgelist -> edge -> firstseg_total_length) / length;
				outputmultitranscript(fp, edgelist -> edge -> seq, length, thisedge -> firstseg_num_reads, avecov);
				edgelist -> edge -> parttag = 1;
			}
		}
		if(edgelist -> edge -> begtag != 1 && edgelist -> edge -> tag >= 2)	{
			num_updatededge = traversegraph(edgelist -> edge, updatededge, num_updatededge, fp, 1);
		}
		edgelist = edgelist -> next;
	}
	return num_updatededge;
}

int taag::extendseq(char *seq, int length, EDGE *thisedge, int *endnodeindex, int *num_reads, int *total_length, char tag)
{
	int	i, j, k, l, m, n, m1;
	int	index;
	int	indegree, outdegree;
	int	*nextedgeindex, *lastedgeindex;
	char	termtag;

	if(tag == 0)	{
		*num_reads += thisedge -> lastseg_num_reads;
		*total_length += thisedge -> lastseg_total_length;
	} else	{
		*num_reads += thisedge -> num_reads;
		*total_length += thisedge -> total_length;
	}
	termtag = 0;
	if(length == 0)	{
		if(thisedge -> begtag == 3)	{
			m = 0;
		} else if(thisedge -> lastseg < thisedge -> length)	{
			m = thisedge -> lastseg - 1;
		} else	{
			thisedge -> tag = 0;
			return(length);
		}
		m1 = thisedge -> length - kmersize;
	} else	{
		m = 0;
		if(thisedge -> endtag == 1 && thisedge -> firstseg > 0)	{
			termtag = 1;
			if(thisedge -> parttag == 0)	{
				m1 = thisedge -> firstseg;
			} else	{
				m1 = kmersize;
			}
		} else	{
			m1 = thisedge -> length - kmersize;
		}
	}
	for(i = m; i < m1; i ++)	{
		seq[length ++] = thisedge -> seq[i];
	}
	if(termtag == 1)	{
		thisedge -> parttag = 1;
		if(thisedge -> firstseg >= thisedge -> length)	{
			thisedge -> lastseg = thisedge -> length - kmersize;
		}
		return(length);
	}

	thisedge -> tag = 1;
	//edge -- all edge
	if(&(edge[thisedge -> bal_edge]) != thisedge)	{
		if(edge[thisedge -> bal_edge].tag == 3 && stranded == 0)	{
			edge[thisedge -> bal_edge].tag = 2;
		}
	}

	//lastedgeindex = (int *) ckalloc((vertex[edge -> nodeindex[1]].indegree + 10) * sizeof(int));
	lastedgeindex = new int[vertex[thisedge -> nodeindex[1]].indegree + 10];
	indegree = chk_lastedge(&vertex[thisedge -> nodeindex[1]], lastedgeindex);
	//nextedgeindex = (int *) ckalloc((vertex[thisedge -> nodeindex[1]].outdegree + 10) * sizeof(int));
	nextedgeindex = new int[vertex[thisedge -> nodeindex[1]].outdegree + 10];
	outdegree = chk_nextedge(&vertex[thisedge -> nodeindex[1]], nextedgeindex);

	if(outdegree != 1 || indegree > 1)	{
		for(i = thisedge -> length - kmersize; i < thisedge -> length; i ++)	{
			seq[length ++] = thisedge -> seq[i];
		}
		*endnodeindex = thisedge -> nodeindex[1];
	} else	{
		//edge -- all edge
		if(edge[nextedgeindex[0]].tag != 0)	{
			length = extendseq(seq, length, &(edge[nextedgeindex[0]]), endnodeindex, num_reads, total_length, 1);
		} else	{
			for(i = thisedge -> length - kmersize; i < thisedge -> length; i ++)	{
				seq[length ++] = thisedge -> seq[i];
			}
		}
	}
	delete[] lastedgeindex;
	delete[] nextedgeindex;
	return(length);
}

int taag::chk_lastedge(VERTEX *thisvertex, int *lastedgeindex)
{
	int	indegree;
	EDGELIST *edgelist;

	indegree = 0;
	edgelist = thisvertex -> lastedge;
	while(edgelist)	{
		if(edgelist -> edge -> tag != 0)	{
			lastedgeindex[indegree ++] = edgelist -> edge -> index;
		}
		edgelist = edgelist -> next;
	}
	return(indegree);
}

int taag::chk_nextedge(VERTEX *thisvertex, int *nextedgeindex)
{
	int	outdegree;
	EDGELIST *edgelist;

	outdegree = 0;
	edgelist = thisvertex -> nextedge;
	while(edgelist)	{
		if(edgelist -> edge -> tag != 0 && (edgelist -> edge -> tag == 3 || (edgelist -> edge -> parttag == 0 && edgelist -> edge -> endtag == 1)))	{
			nextedgeindex[outdegree ++] = edgelist -> edge -> index;
		}
		edgelist = edgelist -> next;
	}
	return(outdegree);
}
