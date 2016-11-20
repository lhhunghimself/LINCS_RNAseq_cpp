#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);  

//takes the first nLines from a large fastq file 
//takes two parameters fastq file and number of lines to extract
  
int main(int argc, char *argv[]){
 gzFile fp1=0;  
 kseq_t *seq;
 char *inputFile=0;
 int opt,maxLines=10000;
 maxLines=atoi(argv[2]);
 inputFile=(argv[1]);
          
 if(!(fp1 = gzopen(inputFile, "r"))){
	 fprintf(stderr,"unable to open R1 file %s\n",inputFile),	 
		exit(0);
	}	
 seq = kseq_init(fp1);
 int nLines=0;
 while ( kseq_read(seq) >= 0 && nLines < maxLines) { // STEP 4: read sequence  
  printf("@%s\n%s\n+\n%s\n",seq->name.s,seq->seq.s,seq->qual.s);
  nLines++;
 }	
 kseq_destroy(seq);
 gzclose(fp1);
 return 0;  
}
