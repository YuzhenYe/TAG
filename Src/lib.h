#ifndef __LIB_H_
#define __LIB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#define pow1(x) ((x) * (x))
#define rev(x) ((2 + (x)) % 4)
/*
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))
*/
#define numc(i, j) ((i) < (j)? (j)*((j)+1)/2+(i): (i)*((i)+1)/2+(j))
#define LARGENUMBER 100000000
#define MIN_CON_LEN 2000000
#define MIN_OVERLAP 5
#define MIN_LEG 1000
#define MBIG 1000000000
#define MSEED 16183398
#define MZ 0
#define FAC (1.0/MBIG)
#define MAXD 5

#define total_nuc 16
#define na_name "actgnrywsmkhbvdx"
#define code_name "ACTG"

typedef struct POSsegment {
	char	name[30];
	int	pos[2];
} SEGMENT;

typedef struct EDGEseg {
	int	index; //which edge
	int	pos[2];
	int	length; //number of edges a mapping spans
	int	num_reads; //number of reads supporing the mapping 
	int	total_length; //length of the mapping
	struct EDGEseg *next;
} EDGESEG;

typedef struct ContigNAMES {
	int	index;
	int	rank;
	char	name[100];
} CONTIGNAMES;

typedef struct edge {
	unsigned int	tag : 4 ; 	//4 bits for tag edges to be removed for various reasons
	unsigned int	begtag : 2 ; 	//for tag edges at the begining of a path in transcript mapping
	unsigned int	endtag : 2 ; 	//for tag edges at the end of a path in transcript mapping
	unsigned int	parttag : 2 ; 	// for tag edges in that the first segment has been output
	bool	multitag;
	int	firstseg, lastseg;	// the end position of the first segment and the first position of the last segment spanning multiple edges
	int	firstseg_num_reads, lastseg_num_reads;
	int	firstseg_total_length, lastseg_total_length;
	int	num_reads, total_length;
	int	index;		/* index of the edge	*/
	char	*seq;
	char	covertag;
	int	length;
	int	bal_edge;
	int	nodeindex[2];
} EDGE;

typedef struct edgelist {
	EDGE	*edge;
	struct edgelist	*next;
} EDGELIST;

typedef struct vertex {
	unsigned int	tag : 1; // tagging visited vertices
	__int128_t	index;
	int	indegree, outdegree;
	EDGELIST	*lastedge;
	EDGELIST	*nextedge;
} VERTEX;

typedef struct hashpos {
	int readindex, pos;
	struct hashpos *next;
} HASH;

using namespace std;

#endif
