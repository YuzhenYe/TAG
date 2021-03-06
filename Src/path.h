#ifndef __PATH_H_
#define __PATH_H_

#include "lib.h"

namespace path
{
	int reverse_path(EDGESEG *path, EDGESEG *rev_path, int len_path, EDGE *edge);
	EDGESEG *insert_single_path(EDGESEG *pathlist, EDGESEG *newpathlist);
	int ReadSingleSegment(bool stranded, char *inpfile, EDGESEG **pathlist, EDGE *edge, int num_edge);
	int ReadPairSegment(bool stranded, int insert_size, char *inpfile, char *pairfile, EDGESEG **pathlist, EDGE *edge, int num_edge);
	int ReadPairedSam(bool stranded, int insert_size, char *samfile, EDGESEG **pathlist, EDGE *edge, int num_edge);
	bool ReadSamOneLine(const char *str, int *pos, int &cindex, char &strand, int &mate);
	int ProcessOneAln(EDGESEG **pathlist, EDGE *edge, bool stranded, int cindex, char strand, int *alnpos);
	int countseg(char *inpfile);
	
	int parseindex(char *name, char *lastname, int *bit);
	EDGESEG* insert_path(EDGESEG *path, int len_path, EDGESEG *pathlist);
	int reverse_path(EDGESEG *path, EDGESEG *rev_path, int len_path, EDGE *edge);
	void delete_path(EDGESEG *pathlist);
}

#endif
