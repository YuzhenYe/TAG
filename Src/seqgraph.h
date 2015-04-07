#ifndef __MAP2GRAPH_H_
#define __MAP2GRAPH_H_

#include "lib.h"

class seqgraph {
	protected: 
        	int inppos;
        	int ngap; //number of gaps
		int maxmut;
        	int kmersize;
        	int overlaplen;
		__int128_t maskv;
		bool stranded;
	
		int len_seq;
		char *seq;
		char seqname[100];
		char seqname_new[100];

		int num_vertex;
		int *vindex;
		VERTEX *vertex; //??

		int num_edge;
		EDGE *edge; //contigs -- contain seq etc
		int num_path;

		int len_path;
		EDGESEG **pathlist;
		int pathstats[15];

		int tot_reads_mapped;

		int hashw;
		int hashtablesize;

	public:
		seqgraph(void);
		~seqgraph(void) { cleanup(); }
		void init();
		void set_maskv(void);
		void cleanup(void);
		void set_kmersize(int inp) { kmersize = inp; set_maskv(); }
		void set_overlaplen(int inp) { overlaplen = inp; }
		void set_ngap(int inp) { ngap = inp; }
		void set_maxmut(int inp) { maxmut = inp; }
		void set_stranded_false(void) { stranded = false; }

		void loadsoap(bool SOAP2, char *edgefile, char *edgeseqfile);
		void index_vertex(void);
		void print_vertex_degree(void);

		void loadsam(char *samfile, int insert_size);
		void loadmap(char *segfile, char *pairfile, int insert_size);
		void read2graph(char *inpfile, char *pairinpfile, bool fastq);
		int  map1read(char *seq, int len_seq, EDGESEG* path);
		int  intalign(char *seq, int length, VERTEX *vertex1, VERTEX *vertex2, EDGESEG *tmppath, int len_tmppath, int *rmut);
		int  forwardalign(char *seq, int length, VERTEX *vertex0, EDGESEG *tmppath, int len_tmppath,int *rmut);
		int  for_align(char *seq1, int length1, char *seq2, int length2, int *mut, int *pos);
		int  backalign(char *seq, int length, VERTEX *vertex0, EDGESEG *tmppath, int len_tmppath, int *rmut);
		int  rev_align(char *seq1, int length1, char *seq2, int length2, int *mut, int *pos);
};

#endif
