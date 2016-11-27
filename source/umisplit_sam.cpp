#include <zlib.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>  
#include <omp.h>
#include <glob.h>
#include "umitools.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

extern "C" {
 #include "optparse.h"  
}
  
int main(int argc, char *argv[]){
 char verbose=0;  
 unsigned int barcodeSize=0;
 char *arg=0,*outputFileName=0;
 string barcodeFileName;
 struct optparse options;
 int opt;
 int mismatchTol=0,NTol=0;
 optparse_init(&options, argv);
 int nThreads=1;
 int filter=0;
 string outputDir="output",inputDir="";
 vector<string> inputFiles;
string errmsg="umisplit_sam h?vfi:t:m:N:o:b:\n-h -? (display this message)\n-v (Verbose mode)\n-f (filter and discard reads with ambiguous UMI and barcode (default is to keep))\n-i <input directory - instead of or in conjunction of entering sam files - an input directory can be provided with the sam files>):\n-t <number of threads(1)>\n-q <minimum quality of UMI base pair before changed to N (10)>-b <barcode file>\n-m <maximum number of mismatches tolerated in barcode (0)>\n-N <maximum number of ambiguous base pairs tolerated in barcode (0)>\n-o <Output Directory (output)>\n<samfile1>..<samfileN>\n\nRequired parameters: barcodefile, input directory and/or individual sam files\n\nExample:\numisplit_sam -b References/Broad_UMI/barcodes_trugrade_96_set4.dat -i Aligns -o well_sorted_Aligns\n";
 
 //parse flags
 while ((opt = optparse(&options, "vfi:t:m:N:o:b:")) != -1){
  switch (opt){
			case 'v':
			 verbose=1;
			break;
			case 'i':
			 inputDir=string(options.optarg);
			break;	
			case 'f':
    filter=1;
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
 while ((arg = optparse_arg(&options))){
	 inputFiles.push_back(string(arg));
	}
	if(barcodeFileName == ""){
		fprintf(stderr,"a file with barcodes is required\n");
		exit(EXIT_FAILURE);
	}	
	if(inputFiles.size() && inputDir != ""){
		fprintf(stderr,"can not enter inputFiles and use -d flag at the same time\n");
		exit(EXIT_FAILURE);
	}
	
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		if(inputDir != "") fprintf(stderr,"Input directory is %s\n",inputDir.c_str());
		fprintf(stderr,"Maximum number of mismatched bases in match %d\n",mismatchTol);
		fprintf(stderr,"Maximum number of unknown bases in match %d \n",NTol);
		fprintf(stderr,"Number of threads is %d\n",nThreads);
		fprintf(stderr,"barcodeFile is %s\n",barcodeFileName.c_str());
		if(filter)fprintf(stderr,"filtering out unmatched barcodes \n");
	 else (stderr,"saving unmatched barcodes to 'X' directory \n");
	 if(inputDir != ""){
			fprintf(stderr,"will look for sam files in directory %s\n",inputDir.c_str());
		}	
  for (int i=0;i<inputFiles.size();i++){
		 fprintf(stderr,"will open file %s\n",inputFiles[i].c_str());		
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
	//get sam files from directory if directory is given
	if(inputDir != ""){
		glob_t glob_result;
	 string globString=inputDir+"/*.sam";
		if(verbose) fprintf(stderr,"using query %s to search for files\n",globString.c_str()); 
		glob(globString.c_str(),GLOB_TILDE,NULL,&glob_result); 
	 for(unsigned int i=0; i<glob_result.gl_pathc; ++i) {
			if(verbose) fprintf(stderr,"adding %s to list of input files\n",glob_result.gl_pathv[i]);
		 inputFiles.push_back(string(glob_result.gl_pathv[i]));
		}
	}	
	const int wellSequenceSize=barcodePanel[0]->barcodeSize;
 #pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i+=2){
	 std::ifstream fstream(inputFiles[i]);
		if(!fstream.is_open()){		fprintf(stderr,"unable to open file %s\n", inputFiles[i].c_str());continue;}
		const int tid=omp_get_thread_num();
		fs::path inputPath (inputFiles[i]);
		string inputFilestem=inputPath.stem().string();
		
  FILE *ofp,*ofps[NWELLS+1];	
  for(int j=0;j<NWELLS+1;j++){
			ofps[j]=0;
			string file;
			if(!filter && !j){
				file=outputDir+"/X"+"/"+inputFilestem+".X"+".sam";
	   ofps[0]=fopen(file.c_str(),"w");			
	   if(!ofps[j]){
				 fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
			 }	
			}	
			else if(j){
    file =outputDir+"/"+barcodePanel[tid]->wells[j-1]+"/"+inputFilestem+"."+barcodePanel[tid]->wells[j-1]+".sam";
    ofps[j]=fopen(file.c_str(),"w");
			 if(!ofps[j]){
				 fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
			 }				
			}
		}
		string fullLine;
		while(getline(fstream,fullLine)){
		 if(fullLine[0] != '@'){
				vector <string> items;
				vector <string> tempItems;
				splitStr(fullLine," \t",items);
				splitStr(items[0],":",tempItems);
	 	 string fullBarcode=tempItems[tempItems.size()-1];
	 	 int wellIndex=barcodePanel[tid]->bestMatch(fullBarcode.c_str());
	 	 if(!wellIndex && filter)continue;
				fputs(fullLine.c_str(),ofps[wellIndex]);		 
			}
		}
		fstream.close();
		if(!filter)
		 fclose(ofps[0]);		   
		for(int j=1;j<NWELLS+1;j++)
		 fclose(ofps[j]); 	
	}
	for(int i=0;i<nThreads;i++)
  if(barcodePanel[i])delete barcodePanel[i];
 delete[] barcodePanel;
 return 0;  
}
