#ifndef __CONTIGAPP_H_
#define __CONTIGAPP_H_

#include "lib.h"

namespace contigapp
{
	int findcontig(char *contig, CONTIGNAMES *NameOfContigs, int nc, int num_contigs, int maxc);
	int binarysearchindex(CONTIGNAMES *NameOfContigs, int beg, int end, int pos);
}

#endif
