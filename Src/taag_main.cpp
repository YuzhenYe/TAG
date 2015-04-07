//Package: tag (version 0.91)
//Program Name: taag_main (the main program) 
//Latest modification on Jan 22th, 2014
//Yuzhen Ye <yye@indiana.edu> and Haixu Tang <hatang@indiana.edu>
//Affiliation: School of Informatics and Computing, Indiana University, Bloomington
#include "taag.h"
#include <time.h>

void printusage(char *error);

int main(int argc, char **argv)
{
	int copt, inp;
	bool SOAP2 = true;
	bool stranded = true;
	int insert_size = 600;
	bool fastq = false;
	int output_minlen = 100;
	bool sepmult = false;
	bool savegraph = false;

	time_t t0 = time(NULL);
	taag *eng = new taag();
	char samfile[1000], read1file[1000], read2file[1000], edgefile[1000], edgeseqfile[1000];
	samfile[0] = read1file[0] = read2file[0] = edgefile[0] = edgeseqfile[0] = 0; 
	char outfile[1000], singleoutfile[1000];
	strcpy(outfile, "tmp.out");

	//get options
	while ((copt=getopt(argc,argv,"i:j:m:e:o:s:k:l:g:d:z:r:vnuqxy")) != EOF)	{
		switch(copt) {
			case 'i':
			  sscanf(optarg,"%s", read1file);
			  continue;
			case 'j':
			  sscanf(optarg,"%s", read2file);
			  continue;
			case 'm':
			  sscanf(optarg,"%s", samfile);
			  continue;
			case 'e':
			  sscanf(optarg,"%s", edgefile);
			  continue;
			case 'n':
			  stranded = false;
			  eng -> set_stranded_false();
			  continue;
			case 's':
			  sscanf(optarg,"%s", edgeseqfile);
			  continue;
			case 'o':
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'k':
			  sscanf(optarg,"%d", &inp);
			  eng -> set_kmersize(inp);
			  continue;
			case 'l':
			  sscanf(optarg,"%d", &inp);
		          eng -> set_overlaplen(inp);
			  continue;
			case 'z':
			  sscanf(optarg,"%d", &insert_size);
			  continue;
			case 'q':
			  fastq = true;
			  continue;
			case 'u':
			  SOAP2 = false;
			  continue;
			case 'g':
			  sscanf(optarg,"%d", &inp);
		  	  eng -> set_ngap(inp);
			  continue;
			case 'd':
			  sscanf(optarg,"%d", &inp);
		  	  eng -> set_maxmut(inp);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &output_minlen);
			  continue;
			case 'x':
			  sepmult = true;
			  continue;
			case 'y':
			  savegraph = true;
			  continue;
			default:
			  printusage("unknown input");
			}
	
		optind--;
	}

	if(read1file[0] == 0) {
		printusage("read1file not specified");
	}
	if(samfile[0] == 0) {
		printusage("samfile not specified");
	}

	//load assembly graph: inputs, edgefile & edgeseqfile
	eng -> loadsoap(SOAP2, edgefile, edgeseqfile);

	//map reads to assembly graph & save results to outfile, singleoutfile
	//inputs: read1file & read2file (reads file, paired ends if read2file given)
	//inputs: samfile (mapped segments, can be paired)
	eng -> loadsam(samfile, insert_size);

        time_t t1 = time(NULL);
	eng -> read2graph(read1file, read2file, fastq);
        time_t t2 = time(NULL);

	//assembly the transcripts based on the mapping results
	eng -> assembly();

	//write transcript assemblies
	eng -> writefile(outfile, sepmult, savegraph, output_minlen);

	delete eng;

        time_t tn = time(NULL);
	if(int(difftime(t2, t1)/60) == 0) {
	        printf("Total time used: %d sec (%d sec for reads mapping to graph)\n",int(difftime(tn, t0)), int(difftime(t2, t1)));
	}
	else {
	        printf("Total time used: %d min (%d min for reads mapping to graph)\n",int(difftime(tn, t0)/60), int(difftime(t2, t1)/60));
	}
}

void printusage(char *error)
{
	
	cout<<"Error: "<<error<<endl;
	cout<<"tag version 0.91"<<endl;
	cout<<"Usage: tag -i read1file -j read2file -m samfile -e edgefile -s edgeseqfile -o -k kmersize -o outfile\n";
	cout<<"-i read1file: The input file name of 1st reads"<<endl;
	cout<<"-j read2file: The input file name of 2nd reads"<<endl;
	cout<<"-m samfile: sam file"<<endl;
	cout<<"-e edgeFile: The input edge file name"<<endl;
	cout<<"-s edgeSeqFile: The input edge sequence (contig) file name"<<endl;
	cout<<"-o OutFile(base name only): Transcript files"<<endl;
	cout<<"-k kmersize (must be the same as the kmer size used for the assembly of matched metagenome2!!)"<<endl;
	cout<<"-l minimum-transcript-len (default: 100)"<<endl;
	cout<<"\nOther options"<<endl;	
	cout<<"-q (input fastq files when set)"<<endl;
	cout<<"-l minimum overlap length"<<endl;
	cout<<"-g maximum allowed gaps"<<endl;
	cout<<"-d maximum mutations (mismatch + indel)"<<endl;
	cout<<"-z maximum insert size"<<endl;
	cout<<"-n (default: no tag, RNA-seq reads were acquired by using a strand-specific protocol)"<<endl;
	cout<<"-x (default: off, transcripts spanning multiple edges output to the same file as other transcripts)"<<endl;
	cout<<"-y (default: off, no outputs for transcript graph)"<<endl;
	cout<<"-u (SOAP when set; default off for SOAP2)"<<endl;
	exit(-1);
}
