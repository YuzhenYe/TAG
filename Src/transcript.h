#ifndef __TRANSCRIPT_H_
#define __TRANSCRIPT_H_

#include "lib.h"

namespace transcript
{
	int output_vertex(FILE *fp, VERTEX vertex);
	int getnum(EDGESEG *pathlist);
	int getlist(EDGESEG **pathindex, EDGESEG *pathlist);
	int comparseg(EDGESEG **s1, EDGESEG **s2);
	int chk_nextedge(VERTEX *vertex, int *nextedgeindex);
	int output_transcript(FILE *fp, char *seq, int index, int pos1, int pos2, int num_reads, double avecov);
	int traversegraph(EDGE *alledge, EDGE *edge, EDGE *updatededge, int num_updatededge, VERTEX *vertex, FILE *fp, int *num_singleedge, int kmer, char tag);
	int extendseq(char *seq, int length, EDGE *edge, EDGE *alledge, VERTEX *vertex, int *endnodeindex, int *num_reads, int *total_length, int kmer, char tag);
	int output_transcriptgraph(char *outfile, char *singleoutfile, EDGE *edge, int num_edge, VERTEX *vertex, int num_vertex, EDGESEG **pathlist, int kmer);
};

#endif
