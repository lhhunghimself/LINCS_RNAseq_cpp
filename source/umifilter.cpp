#include <stdio.h>
#include <iostream>
#include <string.h>  
#include "umitools.hpp"

extern "C" {
 #include "optparse.h"  
} 
int main(int argc, char *argv[]){
 vector<string> inputFiles;
 int opt;
 char *arg=0;
 bool verbose;
 int maxCounts=1000;
 unsigned int umicounts[NUMIS];
 string errmsg="-m <maxCounts>\n-c countFile\n<inputFiles>\n";
 //parse flags
	struct optparse options;
 optparse_init(&options, argv);	 
 
 while ((opt = optparse(&options, "h?m:c:")) != -1) {
  switch (opt){
			case 'v':
			 verbose=1;
			break;
			case 'c':
				if(!readCountsFile(options.optarg,umicounts)){
  		 fprintf(stderr,"unable to read %s\n",options.optarg);
	    exit(1);
				}
			break;
   case 'm':
    maxCounts=atoi(options.optarg);
   break;  
   case '?':
    fprintf(stderr, "%s parameters are: %s\n", argv[0], errmsg.c_str());
    exit(0);
   break;
   case 'h':
    fprintf(stderr, "%s parameters are: %s\n", argv[0], errmsg.c_str());
    exit(0);
   break; 
  }
 }
 //parse file arguments
 //these will be R1 followed by R2 arguments
 while ((arg = optparse_arg(&options))){
	 inputFiles.push_back(string(arg));
	}

	for(int i=0;i<inputFiles.size();i++){
		char line[1024]; //max lines so 813
	 FILE *fp=fopen(inputFiles[i].c_str(),"r");
	 if(!fp){
	 	fprintf(stderr,"unable to read %s\n",inputFiles[i].c_str());
		 return 0;
		}
		while(fgets(line,sizeof(line),fp)){
			//remember \n at the end of line 1
			string barcode(line+(strlen(line)-1-UMISIZE),UMISIZE);
			if(umicounts[hashCode4(barcode)] <= maxCounts){
    fputs(line,stdout);
			 fputs(fgets(line,sizeof(line),fp),stdout);
			 fputs(fgets(line,sizeof(line),fp),stdout);
			 fputs(fgets(line,sizeof(line),fp),stdout);
			}
			else{
    fgets(line,sizeof(line),fp);
    fgets(line,sizeof(line),fp);
    fgets(line,sizeof(line),fp);
			}	
		}	
		fclose(fp);
	}	
	return 0;	
}


