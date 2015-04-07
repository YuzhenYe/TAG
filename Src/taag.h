#ifndef __TAAG_H_
#define __TAAG_H_

#include "seqgraph.h"

struct sorter {
  	bool operator() (const EDGESEG* s1, const EDGESEG* s2) { 
		if(s1->pos[0] == s2->pos[0]) {
			return s1->pos[1] < s2->pos[1];
		}
		else {
			return s1->pos[0] < s2->pos[0];
		}
		//return (s1->pos[0] < s2->pos[0]);
   	}
};

class taag:public seqgraph
{
	private:
		int	multtranscript;  //number of multi-edge transcripts
		int	multtranscript_long; //number of multi-edge transcripts of >= minlen
		long	multtranscriptlen; //total length of multi-edge transcripts	
		long	multtranscriptlen_long;		
		int	transcriptlen_min; //only output transcripts of this length
		int	multtranscriptlen_max;
	public:	
		void assembly(void);
		void writefile(char *outfile, bool sepmult, bool savegraph, int minlen);
		int get_path_num(EDGESEG *pathlist_t);
		//int getlist(EDGESEG **pathindex, EDGESEG *pathlist_t);
		//int get_path_list(vector<EDGESEG*> &pathindex, EDGESEG *pathlist_t);
		int get_path_list(EDGESEG** pathindex, EDGESEG *pathlist_t);
		//int comparseg(const void * p1, const void * p2);
		bool comparseg(const EDGESEG *s1, const EDGESEG *s2);
		void output_vertex(ofstream &fp, VERTEX vertex);
		void output_updatededge(ofstream &fp, EDGE updatededge);
		void output_transcript(ofstream &fp, char *seq, int index, int pos1, int pos2, int num_reads, double avecov);
		void outputmultitranscript(ofstream &fp, char *seq, int length, int num_reads, double avecov);
		void outputpartialtranscript(ofstream &fp, char *seq, int index, int length, int num_reads, double avecov, bool multitag);
		int traversegraph(EDGE *thisedge, EDGE *updatededge, int num_updatededge, ofstream &fp, char tag);
		int extendseq(char *seq, int length, EDGE *thisedge, int *endnodeindex, int *num_reads, int *total_length, char tag);
		int chk_lastedge(VERTEX *thisvertex, int *lastedgeindex);
		int chk_nextedge(VERTEX *thisvertex, int *nextedgeindex);
		void set_minlen(int len);
};


#endif
