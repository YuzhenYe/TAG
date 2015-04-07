#ifndef __READS_H_
#define __READS_H_
#include "lib.h"

namespace reads
{
	int countreads(char *inpfile, bool fastq);
	int get1readname(ifstream &fp, char *seqname);
	int get1readfasta(ifstream &fp, char *seq, char *seqname);
	int get1readfastq(ifstream &fp, char *seq, char *seqname);
};

#endif
