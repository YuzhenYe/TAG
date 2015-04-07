#include "contigapp.h"

int contigapp::findcontig(char *contig, CONTIGNAMES *NameOfContigs, int nc, int num_contigs, int maxc)
{
	int	i, j, k;

	if(contig[0] == 'C')	{
		k = atoi(&contig[1]);
		j = binarysearchindex(NameOfContigs, 0, nc - 1, k);
		return(NameOfContigs[j].index);
	} else	{
		k = atoi(&contig[8]);
		j = binarysearchindex(NameOfContigs, nc, num_contigs - 1, k + maxc + 1);
		return(NameOfContigs[j].index);
	}
}

int contigapp::binarysearchindex(CONTIGNAMES *NameOfContigs, int beg, int end, int pos)
{
	int	i, j, k;
	int	mid;

	if(end < beg)	{
		printf("position not found: %d\n", pos);
		exit(-1);
	} 
	else if(end == beg)	mid = end;
	else mid = beg + (end - beg) / 2 + 1;
	if(NameOfContigs[mid].rank == pos)	{
		return(mid);
	} else if(NameOfContigs[mid].rank < pos)	{
		binarysearchindex(NameOfContigs, mid + 1, end, pos);
	} else	{
		binarysearchindex(NameOfContigs, beg, mid - 1, pos);
	}
}
