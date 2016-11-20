#include <zlib.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>  
#include <omp.h>
#include "kseq.h"
#include "umitools.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

//this follows the python script
extern "C" {
 #include "optparse.h"  
}
#define R1_LENGTH 16 //default UMI length
string findBaseName(string name);
KSEQ_INIT(gzFile, gzread);  

int main(int argc, char *argv[]){

 char verbose=0,barcode=0;  
 int l1,l2,emptySeqs=0,nameMismatch=0,minQual=10,UMILength=R1_LENGTH,nArgs=0;
 unsigned int barcodeSize=0,minDist=4;
 int adjustedQ=minQual+33;
 char *arg=0,*outputFileName=0;
 string barcodeFileName;
 struct optparse options;
 int opt;
 int mismatchTol=1,NTol=2;
 optparse_init(&options, argv);
 int nThreads=1;
 int filter=0;
 string outputDir;
 vector<string> inputFiles;
 int maxSeqs=5000000; //number of sequence per fastq file
 
 //parse flags
	string errmsg="umisplit_orig vh?s:l:o:t:q:\n-v (Verbose mode)\n-s <Maximum number of reads in a fastq file (5000000)>\n-l <Length of UMI (16):\n-t <number of threads(1)>\n-q <minimum quality of UMI base pair before changed to N (10)>\n-o <Output Directory>\n-h -? (display this message)\n<R1file.1> <R2file.1>..<R1file.N> <R2file.N>\n\nExample:\numisplit_orig -o Aligns sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz\n";
	 
 while ((opt = optparse(&options, "vh?s:l:o:t:q:")) != -1) {
  switch (opt){
			case 'v':
			 verbose=1;
			break;	
			case 's':
			 maxSeqs=atoi(options.optarg);
			break;
			case 'l':
    UMILength=atoi(options.optarg);//not the well barcode
   break;
			case 'o':
			 outputDir=string(options.optarg);
			break;
   case 't':
    nThreads=atoi(options.optarg);
   break;         
   case 'q':
    minQual=atoi(options.optarg);
    adjustedQ=33+minQual;
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

 //these will be R1 followed by R2 arguments
 while ((arg = optparse_arg(&options))){
	 inputFiles.push_back(string(arg));
	}
 //extract fileName
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"UMI Length %d\n",UMILength);
		fprintf(stderr,"Threads used %d\n",nThreads);
		fprintf(stderr,"maxnumber of sequences to write out %d\n",maxSeqs);
		fprintf(stderr,"Minimum quality %d\n",minQual);
		int i=0;
		while (i<inputFiles.size()){
		 fprintf(stderr,"R1 file %s\n",inputFiles[i++].c_str());
		 if(i == inputFiles.size()){
				fprintf(stderr,"missing corresponding R2 file\n");
				exit(EXIT_FAILURE);
			}
		fprintf(stderr,"R2 file %s\n",inputFiles[i++].c_str());		
		}		
	}
 
#pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i+=2){
		const int tid=omp_get_thread_num();
		string basename=findBaseName(inputFiles[i]);
		gzFile fp1=0, fp2=0;  
  kseq_t *seq1,*seq2;
  int l1,l2;
		fp1 = gzopen(inputFiles[i].c_str(), "r");
		if(fp1 == Z_NULL){
			fprintf(stderr,"unable to open %s\n",inputFiles[i].c_str());
			exit(EXIT_FAILURE);
		}
		fp2 = gzopen(inputFiles[i+1].c_str(), "r");	
		if(fp2 == Z_NULL){
			fprintf(stderr,"unable to open %s\n",inputFiles[i+1].c_str());
			exit(EXIT_FAILURE);
		}
  seq1 = kseq_init(fp1);
  seq2 = kseq_init(fp2);
  
  string outputFile=outputDir+"/"+basename+"temp." +std::to_string(tid);
  FILE *ofp=fopen(outputFile.c_str(),"w");
  if(!ofp)exit(EXIT_FAILURE);
  fprintf(stderr,"opened %s\n",outputFile.c_str());
	 int nSeqs=0;
  while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >=0) {
	  //check for errors 
		 if(l1==0 || l2 == 0 || strcmp(seq1->name.s,seq2->name.s)){
		 	if(l1==0 || l2 == 0) emptySeqs++;
		 	if(strcmp(seq1->name.s,seq2->name.s)){
		 		nameMismatch++;
		 		if(verbose)fprintf(stderr,"Warning - mismatch of names in R1 and R2\n");
		 	}	
		 	continue;
		 }
			char *cptr=seq1->seq.s;
		 char *qptr=seq1->qual.s;
		 int k=0; 
			while(*cptr && k < UMILength){
		  if(*qptr<adjustedQ)*cptr='N';
		  cptr++;qptr++;k++;
		 }
			fputc('@',ofp);
		 fputs(seq2->name.s,ofp);
		 fputc(':',ofp);				
		 fwrite(seq1->seq.s,UMILength,1,ofp);
			fputc('\n',ofp);	
			fputs(seq2->seq.s,ofp);
		 fputs("\n+\n",ofp);		 
   fputs(seq2->qual.s,ofp);
   fputc('\n',ofp);
   nSeqs++;				 
		
	  if(nSeqs && nSeqs%maxSeqs ==0){
	 	 fclose(ofp);
	 	 string newFile=outputDir+"/"+basename+std::to_string(nSeqs)+".fastq";
	 	 fs::path outputPath{outputFile};
	 	 fs::rename(outputPath,newFile);
				ofp=fopen(outputFile.c_str(),"w");
				if(!ofp)exit(EXIT_FAILURE);
		 }
		}
		if(nSeqs && nSeqs%maxSeqs !=0){
			fclose(ofp);
	 	string newFile=outputDir+"/"+basename+std::to_string(nSeqs)+".fastq";
	 	fs::path outputPath{outputFile};
	 	fs::rename(outputPath,newFile);	
		}
		else{
   //clean up empty temp file
			fclose(ofp);
			fs::remove(outputFile);	
		}							
  kseq_destroy(seq1);
  kseq_destroy(seq2);
  gzclose(fp1);
  gzclose(fp2);		 	 	
	}
 return 0;  
}

string findBaseName(string name){
	string basename,namestem;
	fs::path namePath{name};
	if (namePath.extension().string() == "gz" && (namePath.stem().extension().string() == "fastq" || (namePath.stem().extension().string() == ".fq")))   
  namestem=namePath.stem().stem().string();
 else 
  namestem=namePath.stem().string();
 vector <string> split;
 splitStr(namestem,"_",split);
 if(split.size() < 3){
		fprintf(stderr,"format of file is supposed to be sampleName_lane_R(1/2)_fastq(.gz)\n");
		for(int i=0;i<split.size();i++)
		 cerr<<split[i]<<endl;
		exit (EXIT_FAILURE);
	}
	int i=0;
	for (;i<split.size()-3;i++){
		basename=basename+split[i]+"_";	
	}	
	for(;i<split.size()-1;i++)
	basename=basename+split[i]+".";	
	return basename;
}	
