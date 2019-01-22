#include <zlib.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>  
#include <omp.h>
#include "kseq.h"
#include "umitools.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

extern "C" {
 #include "optparse.h"  
}
bool Ncheck(const char *seq, const int size);

#define R1_LENGTH 16 //default UMI length


KSEQ_INIT(gzFile, gzread);  
string errmsg="umisplit h?vft:m:N:o:b:l:q:\n-h -? (display this message)\n-v (Verbose mode)\n-f (filter and discard reads with ambiguous UMI and barcode (default is to keep))\n-l <Length of UMI (16):\n-t <number of threads(1)>\n-q <minimum quality of UMI base pair before changed to N (10)>-b <barcode file>\n-m <maximum number of mismatches tolerated in barcode (0)>\n-c <the file to output the umicounts to - default is UMIcounts.bin in the output directory>\n-N <maximum number of ambiguous base pairs tolerated in barcode (0)>\n-o <Output Directory>\n<R1file.1> <R2file.1>..<R1file.N> <R2file.N>\n\nExample:\numisplit -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -o Aligns sample1_R1.fastq.gz sample1_R2.fastq.gz sample2_R1.fastq.gz sample2_R2.fastq.gz\n";
  
int main(int argc, char *argv[]){
	bool compressFlag=0;
 char verbose=0,barcode=0;
 uint64_t maxSizeKB=0;  
 int l1,l2,emptySeqs=0,nameMismatch=0,minQual=10,UMILength=R1_LENGTH,nArgs=0;
 unsigned int barcodeSize=0;
 int adjustedQ=minQual+33;
 char *arg=0,*outputFileName=0;
 string barcodeFileName,countsFile="";
 struct optparse options;
 int opt;
 int mismatchTol=0,NTol=0;
 optparse_init(&options, argv);
 int nThreads=1;
 int filter=0;
 string outputDir="";
 vector<string> inputFiles;
 
 //parse flags
 while ((opt = optparse(&options, "h?vft:z:s:m:N:o:b:l:q:c:")) != -1) {
  switch (opt){
			case 'v':
			 verbose=1;
			break;
			case 'f':
    filter=1;
   break;
   case 'z':
    compressFlag=1;
   break;
   case 'c':
    countsFile=string(options.optarg);
   break; 
   case 's':
    maxSizeKB=std::stoull(options.optarg);
   break;     
			case 'o':
			 outputDir=string(options.optarg);
			break;
   case 'm':
    mismatchTol=atoi(options.optarg);
   break;
   case 't':
    nThreads=atoi(options.optarg);
   break;         
   case 'N':
    NTol=atoi(options.optarg);
   break;
   case 'b':
    barcodeFileName=options.optarg;
   break;
   case 'l':
    UMILength=atoi(options.optarg);//not the well barcode
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
	if(outputDir=="" || barcodeFileName==""){
		fprintf(stderr,"must give an output directory and barcodeFilename\n");
		exit(EXIT_FAILURE);
	}	
 //parse file arguments
 //these will be R1 followed by R2 arguments
 while ((arg = optparse_arg(&options))){
	 inputFiles.push_back(string(arg));
	}
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"UMI Length %d\n",UMILength);
		fprintf(stderr,"Barcode fileName %s\n",barcodeFileName.c_str());
		fprintf(stderr,"Output directory is %s\n",outputDir.c_str());
		fprintf(stderr,"Minimum quality %d\n",minQual);		
		fprintf(stderr,"Number of threads %d\n",nThreads);		
		fprintf(stderr,"Maximum number of mismatches in barcode %d\n",mismatchTol);		
		fprintf(stderr,"Maximum number of unknown bases in barcode %d\n",NTol);
		if(filter)fprintf(stderr,"filtering out unmatched barcodes \n");
	 else (stderr,"saving unmatched barcodes to 'X' directory \n");
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
 fs::create_directory(fs::system_complete(outputDir));
  //change this if using 384 wells
 #if NWELLS > 256
 umipanel<uint32_t,uint16_t> **barcodePanel=new umipanel<uint32_t,uint16_tr>*[nThreads];
 for(int i=0;i<nThreads;i++)
	 barcodePanel[i]=new umipanel<uint32_t,uint16_t> (barcodeFileName,mismatchTol,NTol);
	#else
	 umipanel<uint32_t,unsigned char> **barcodePanel=new umipanel<uint32_t,unsigned char>*[nThreads];
 for(int i=0;i<nThreads;i++)
	 barcodePanel[i]=new umipanel<uint32_t,unsigned char> (barcodeFileName,mismatchTol,NTol);
	#endif
	
	//create directories for each well
	for(int i=0;i<NWELLS;i++){
		auto p=fs::path(outputDir+"/"+barcodePanel[0]->wells[i]);
		fs::create_directory(fs::system_complete(p));
	}	
	//create bad directory	
	if(!filter){
	 auto p=fs::path(outputDir+"/X");
	 fs::create_directory(fs::system_complete(p));
	}
	
 //keep track of UMIs
 
  unsigned int **umiBarcodes=new unsigned int* [nThreads];
  umiBarcodes[0]=new unsigned int [nThreads*NUMIS]; //NUMIS is 4**10 i.e. 4 to power of the size of the UMI
  memset(umiBarcodes[0],0,nThreads*NUMIS*sizeof(unsigned int));
  for (int i=1;i<nThreads;i++){
	  umiBarcodes[i]=umiBarcodes[i-1]+NUMIS;
	 }


#pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i+=2){
		const uint64_t maxSize=1024*maxSizeKB;
		const int tid=omp_get_thread_num();
		const int wellSequenceSize=barcodePanel[tid]->barcodeSize;
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
  
  fs::path R1(inputFiles[i]);
  fs::path R2(inputFiles[i+1]);
  string R1stem;
  if (R1.extension().string() == "gz" && (R1.stem().extension().string() == "fastq" || (R1.stem().extension().string() == ".fq")))   
   R1stem=R1.stem().stem().string();
  else 
   R1stem=R1.stem().string();
  FILE *ofp,*ofps[NWELLS+1];
  gzFile ofpgz,ofpsgz[NWELLS+1];
  uint64_t fileSizes[NWELLS+1];
  uint16_t numberofFiles[NWELLS+1];
  memset(fileSizes,0,sizeof(fileSizes));
  memset(numberofFiles,0,sizeof(numberofFiles));
  
  for(int j=0;j<NWELLS+1;j++){
			ofps[j]=0;
			ofpsgz[j]=0;
			string file;
			if(!filter && !j){
				file=outputDir+"/X"+"/"+R1stem+"R2_"+"X_"+std::to_string(numberofFiles[0])+".fq";
				numberofFiles[0]++;
				if (compressFlag){
					file=file+".gz";
					ofpsgz[0]=gzopen(file.c_str(),"wb");			
				}
				else{
	    ofps[0]=fopen(file.c_str(),"w");			
	    if(!ofps[j]){
				  fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
			  }
				}	
			}	
			else if(j){
    file =outputDir+"/"+barcodePanel[tid]->wells[j-1]+"/"+R1stem+"R2_"+barcodePanel[tid]->wells[j-1]+"_"+std::to_string(numberofFiles[j])+".fq";
				if (compressFlag){
					file=file+".gz";
					ofpsgz[j]=gzopen(file.c_str(),"wb");			
				}
				else{
     ofps[j]=fopen(file.c_str(),"w");
     if(!ofps[j]){
				  fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
			  }	
				}
    numberofFiles[j]++;
			}
		}

	 
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
		 //adjust quality and find the well

			char *cptr=seq1->seq.s;
		 char *qptr=seq1->qual.s;
		 int k=0; 
			while(*cptr && k < UMILength){
		  if(*qptr<adjustedQ)*cptr='N';
		  cptr++;qptr++;k++;
		 }
		 //check if there is a N
		 if(filter && Ncheck((seq1->seq.s)+wellSequenceSize,UMILength-wellSequenceSize)) continue;
		 const unsigned int barcodeIndex=barcodePanel[tid]->bestMatch(seq1->seq.s);
		 if(filter && !barcodeIndex)continue;
			
			//skip if ambiguous barcode and get unique barcode index from sequence
			string umiBarcode(seq1->seq.s+6,UMISIZE);
			unsigned int  umiBarcodeIndex=0;
			if (filter && ambigCheck(umiBarcode,umiBarcodeIndex)) continue;
			umiBarcodes[tid][umiBarcodeIndex]++;
			if(compressFlag) ofpgz=ofpsgz[barcodeIndex];
			else ofp=ofps[barcodeIndex];
			
		 if(maxSize){
			 uint16_t lineBytes=strlen(seq2->name.s)+strlen(seq2->seq.s)+UMILength+strlen(seq2->qual.s)+7;
			 if (barcodeIndex && lineBytes+fileSizes[barcodeIndex] > maxSize){
					string file =outputDir+"/"+barcodePanel[tid]->wells[barcodeIndex-1]+"/"+R1stem+"R2_"+barcodePanel[tid]->wells[barcodeIndex-1]+"_"+std::to_string(numberofFiles[barcodeIndex])+".fq";
					if (compressFlag){
					 gzclose(ofpgz);
					 file=file+".gz";
					 ofpsgz[barcodeIndex]=gzopen(file.c_str(),"wb");	
					 ofpgz=ofpsgz[barcodeIndex];
					}
					else{
						fclose(ofp);
					 ofps[barcodeIndex]=fopen(file.c_str(),"w");
					 ofp=ofps[barcodeIndex];
					}
	    numberofFiles[barcodeIndex]++;
					fileSizes[barcodeIndex]=0;
				}
			 fileSizes[barcodeIndex]+=lineBytes;
			}
			if (compressFlag){
				char buffer[1024];
				memset(buffer,0,sizeof(buffer));
				buffer[0]='@';
				char *bufPtr=buffer+1;
				char *dest=seq2->name.s;
			 while(*dest){
					*bufPtr++=*dest++;	
				}
			 *bufPtr++=':';
				memcpy(bufPtr,seq1->seq.s,UMILength);
				bufPtr+=UMILength;
				*bufPtr++='\n';
				dest=seq2->seq.s;
				while(*dest){
					*bufPtr++=*dest++;	
				}
				*bufPtr++='\n';
				*bufPtr++='+';
				*bufPtr++='\n';
				dest=seq2->qual.s;
				while(*dest){
					*bufPtr++=*dest++;	
				}
				*bufPtr++='\n';
				gzwrite(ofpgz,buffer,bufPtr-buffer);
			}
			else{
				fputc('@',ofp);
			 fputs(seq2->name.s,ofp);
			 fputc(':',ofp);				
			 fwrite(seq1->seq.s,UMILength,1,ofp);
				fputc('\n',ofp);	
				fputs(seq2->seq.s,ofp);
			 fputs("\n+\n",ofp);		 
	   fputs(seq2->qual.s,ofp);
	   fputc('\n',ofp);
			}
		}
		if(!filter){
			if (compressFlag){
				gzclose(ofpsgz[0]);
			}
			else fclose(ofps[0]);
		}		   
		for(int j=1;j<NWELLS+1;j++){
			if (compressFlag) gzclose(ofpsgz[j]);
			else fclose(ofps[j]);
		}  
  kseq_destroy(seq1);
  kseq_destroy(seq2);
  gzclose(fp1);
  gzclose(fp2);		 	 	
	}
	//combine the UMI counts
	for(int tid=1;tid<nThreads;tid++){
		for(int j=0;j<NUMIS;j++){
			umiBarcodes[0][j]+=umiBarcodes[tid][j];
		}
	}	
	//print out the UMI counts to a binary
	{
		if (countsFile==""){
	  countsFile=outputDir+"/"+"UMIcounts.bin";
		}
	 FILE *fp=fopen(countsFile.c_str(),"w");
		if(!fp){
		 fprintf(stderr,"unable to open %s\n",countsFile.c_str());
			exit(1);
		}
		if(!fwrite(umiBarcodes[0],sizeof(unsigned int)*NUMIS,1,fp)){
			fprintf(stderr,"error in writing UMI counts\n");
			exit(1);
		}			
		fclose(fp); 
	}
 delete[] umiBarcodes[0];
 delete[] umiBarcodes;
 delete[] barcodePanel;
 return 0;  
}

bool Ncheck(const char *seq, const int size){
	for(int i=0;i<size;i++)
	 if(seq[i] == 'N')return 1;
	return 0;
}
