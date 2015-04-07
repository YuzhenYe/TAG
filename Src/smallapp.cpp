#include "smallapp.h"


FILE* smallapp::ckopen(char *name, char *mode)
{
	FILE *fp;

	if ((fp = fopen(name, mode)) == NULL)	{
		printf("Cannot open %s.\n", name);
		exit(-1);
	}
	return(fp);
}


/* ckalloc - allocate space; check for success */

void* smallapp::ckalloc(int amount)
{
	void *p;

	if(amount == 0)	{
		amount = (unsigned) 100;
	}
	if ((p = (void *) calloc( (unsigned) amount, 1)) == NULL)	{
		printf("Ran out of memory.\n");
                printf("There may be errors as follows:\n");
                printf("1) Not enough memory.\n");
                printf("2) The ARRAY may be overrode.\n");
                printf("3) The wild pointers.\n");
                exit(-1);
	}
	return(p);
}

double smallapp::random1(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj, mk;
	int i, ii,k;

	if(*idum < 0 || iff == 0)	{ /* initialization */
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for(i = 1; i < 54; i ++)	{
			ii = (21 + i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if(mk < MZ)	mk += MBIG;
			mj = ma[ii];
		}
		for(k = 1; k <= 4; k ++)	{
			for(i = 1; i <= 55; i ++)	{
				ma[i] -= ma[1 + (i + 30) % 55];
				if(ma[i] < MZ)	ma[i] += MBIG;
			}
		}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if( ++inext == 56) inext = 1;
	if( ++inextp == 56) inextp = 1;
	mj = ma[inext] - ma[inextp];
	if(mj < MZ) mj += MBIG;
	ma[inext] = mj;
	return(mj * FAC);
}

int  smallapp::ran_number(int n, int *idum)
{
	double	     t;
	int          p;

	t = random1(idum);
	if(t == 1.0)	{
		t = t - 1.0e-10;
	}
	p = ( int ) (n * t);
	return( p );
}

int smallapp::hashvalue(int hashw)
{
        int     i, k;
        k = 1;
        for(i = 0; i < hashw; i ++)     {
                k *= 4;
        }
        return(k);
}

int smallapp::comparcontig(const void * p1, const void * p2)
{
        CONTIGNAMES *s1 = (CONTIGNAMES *)p1;
        CONTIGNAMES *s2 = (CONTIGNAMES *)p2;
        if(s1 ->rank > s2 ->rank)       {
                return(1);
        } else if(s1 -> rank < s2 -> rank)      {
                return(-1);
        } else  {
                return(0);
        }
}

int smallapp::comparint(const void * p1, const void * p2)
{
        SEGMENT *s1 = (SEGMENT *)p1;
        SEGMENT *s2 = (SEGMENT *)p2;
        if(s1 -> pos[0] > s2 -> pos[0])   {
                return(1);
        } else if(s1 -> pos[0] < s2 -> pos[0])    {
                return(-1);
        } else  {
                return(0);
        }
}

//sort vertex based on their index (for fast retrival later)
int smallapp::comparvertex(const void * p1, const void * p2)
{
        VERTEX *v1 = (VERTEX *)p1;
        VERTEX *v2 = (VERTEX *)p2;
        if(v1 -> index > v2 -> index)   {
                return(1);
        } else if(v1 -> index < v2 -> index)    {
                return(-1);
        } else  {
                return(0);
        }
}

char smallapp::strand2char(char *strandtmp)
{
        if(strandtmp[0] == '+') {
                return(0);
        } else  {
                return(1);
        }
}

EDGESEG* smallapp::newEDGESEGpt(void)
{
	EDGESEG *pt = new EDGESEG;
	pt -> index = 0;
	pt -> pos[0] = pt -> pos[1] = 0;
	pt -> length = 0;
	pt -> num_reads = 0;
	pt -> total_length = 0;
	pt -> next = NULL;
}

EDGESEG* smallapp::newEDGESEGpt(int index, int pos[2], int length, int num_reads, int total_length, EDGESEG* next)
{
	EDGESEG *pt = new EDGESEG;
	pt -> index = index;
	pt -> pos[0] = pos[0]; pt -> pos[1] = pos[1];
	pt -> length = length;
	pt -> num_reads = num_reads;
	pt -> total_length = total_length;
	pt -> next = next;
}

void smallapp::delEDGESEGpt(EDGESEG *pt)
{
	delete pt;	
	pt = NULL;
}
