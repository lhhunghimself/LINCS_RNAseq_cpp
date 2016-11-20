#include <zlib.h>  
#include <stdio.h>
#include <string.h>
#include <string> 
#include <iostream> 
#include <unordered_map> 
#include <unordered_set> 
#include <set> 
#include <vector> 
#include "kseq.h"
#include <glob.h>
#include "umitools.hpp"

using namespace std;   

int main(int argc, char *argv[]){
	if(argc != 8){
		fprintf(stderr,"%d parms requires params: sample_id, sym2ref, ercc_fasta, barcodes, aligned_dir, dge_dir, loose_barcodes\n",argc);
		for(int i=1;i<argc;i++){
			fprintf(stderr,"%d\t%s\n",i,argv[i]);
		}	
		exit(EXIT_FAILURE);
	}
	int k=1;	
 char *sample_id=argv[k++];
 char *sym2ref=argv[k++];
 char *ercc_fasta=argv[k++];
 char *barcodes=argv[k++];
 char *aligned_dir=argv[k++];
 char *dge_dir=argv[k++];
 char *looseBarcodes=argv[k++];
 
 unordered_map<string,string> refseq_to_gene;
 vector<string>erccList, geneList,wellList;
 unordered_map<string,unsigned int>well_to_index;
 unordered_map<string,unsigned int>ercc_to_index;
 unordered_map<string,unsigned int>gene_to_index;
 vector<string> unknown_list;
 
 string plateID=readWells(barcodes,well_to_index,wellList);
 readERCC(ercc_fasta,erccList);
 readRefseq(sym2ref,refseq_to_gene ,geneList);
 
 for(int i=0;i<geneList.size();i++)
  gene_to_index[geneList[i]]=i;
 for(int i=0;i<erccList.size();i++)
  ercc_to_index[erccList[i]]=i;
 
 string mmId[2]={".unq.",".all."};

 
 //get all the sam files
 glob_t glob_result;
 char glob_str[1024];
 sprintf(glob_str,"%s/*.sam",aligned_dir);
 glob(glob_str,GLOB_TILDE,NULL,&glob_result);

 
 int mismatchTol=0,NTol=0,mismatchSet=0,NTolSet=0;
 
 //check for loose barcodes
 if(looseBarcodes[0] == 't' || looseBarcodes[0] == 'T'){
		if(!mismatchSet) mismatchTol=1;
		if(!NTolSet) NTol=1;
	}	
 #if NWELLS > 256
 umipanel<uint32_t,uint16_t> barcodePanel(string(barcodes),mismatchTol,NTol);  		
	#else
	umipanel<uint32_t,unsigned char> barcodePanel(string(barcodes),mismatchTol,NTol);  		
	#endif
	
	//for hashvalue - gives rise to 96*1024*1024*28000 ~ 42 bits - faster than 20 byte text to hash

	const uint64_t barcodeBase=geneList.size()+erccList.size()+1; //add1 for chrM
	const uint64_t wellBase=(1<<20)*barcodeBase;
	for(int mm=0;mm<2;mm++){
		 //outputfile names
	 string logFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"log.dat";
	 string unknownFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"unknown_list";
	 string totalAlignFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"refseq.total.dat";
	 string umiFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"refseq.umi.dat";
	 string erccFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"spike.total.dat";
	 string erccUMIFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"spike.umi.dat";
	 string wellTotalFile=string(dge_dir)+"/"+string(sample_id)+mmId[mm]+"well_summary.dat";
	 unsigned int 
		 //unq counts
		 count_total_reads = 0,
		 count_assigned_reads = 0,
		 count_assigned_aligned_reads = 0,	
		 count_assigned_mito_reads = 0,
		 count_assigned_mito_umi = 0,
	 	count_assigned_unknown_reads = 0,
		 count_assigned_unknown_umi = 0,
	  
	  //2D arrays for spike ins
	  //arrays for pointers
			*count_spike_total[NWELLS],
			*count_spike_umi[NWELLS],   
			*count_total[NWELLS],
			*count_umi[NWELLS];   
	  
	  //allocate memory for total
			count_spike_total[0]=new unsigned int[erccList.size()*NWELLS];
			count_spike_umi[0]=new unsigned int[erccList.size()*NWELLS];      
			count_total[0]=new unsigned int[geneList.size()*NWELLS];
			count_umi[0]=new unsigned int[geneList.size()*NWELLS];
	
		//set to zero
			
			memset(count_spike_total[0],0,erccList.size()*NWELLS*sizeof(unsigned int));
			memset(count_spike_umi[0],0,erccList.size()*NWELLS*sizeof(unsigned int));      
			memset(count_total[0],0,geneList.size()*NWELLS*sizeof(unsigned int));
			memset(count_umi[0],0,geneList.size()*NWELLS*sizeof(unsigned int));
	
			//initialize pointers
			for(int i=1;i<NWELLS;i++){
			 count_spike_total[i]=count_spike_total[i-1]+erccList.size();
				count_spike_umi[i]=count_spike_umi[i-1]+erccList.size();      
				count_total[i]=count_total[i-1]+geneList.size();
				count_umi[i]=count_umi[i-1]+geneList.size();
			}
					
		//hashes
	 	unordered_set<uint64_t>umi_seen;
	 	unordered_set<string>unknown_set;
	 	unordered_set<string>unknown_umi_seen;
	
	
	 	for(unsigned int i=0; i<glob_result.gl_pathc; ++i) {
		 FILE *fp =fopen(glob_result.gl_pathv[i],"r");
	  if(!fp){
	 		fprintf(stderr,"unable to open file %s\n", glob_result.gl_pathv[i]);
	 		continue;
	 	}

			
		 char line[1024];
	
			while(fgets(line,sizeof(line),fp)){
		  if(line[0] != '@'){
					vector <string> items;
					vector <string> tempItems;
					count_total_reads++;
					splitStr(line," \t",items);
		   //get rid of all the bad counts
	
		   //parse for barcode
	
	 	  splitStr(items[0],":",tempItems);
	 	  string fullBarcode=tempItems[tempItems.size()-1];
	 	  
	 	  //string well=tempItems[tempItems.size()-2];
	 	  int wellIndex=barcodePanel.bestMatch(fullBarcode.c_str());
	 	  if(!wellIndex)continue;
	 	  wellIndex--;
		   string barcode=fullBarcode.substr(6,10); //change if barcode size changes
		   
		   //skip if ambiguous barcode
		   unsigned int barcodeIndex=0;
		   if(ambigCheck(barcode,barcodeIndex))continue;
		   count_assigned_reads++;
		   string aligned_id=items[2];
		   if(aligned_id == "*")continue;
		   string read=items[9];
		   int edit_dist=stoi(splitStrIndex(items[12],":",-1));
		   int best_hits=stoi(splitStrIndex(items[13],":",-1));
		    
	    if(edit_dist > MAX_EDIT_DISTANCE || best_hits > MAX_BEST || polyACheck(read)) continue;
	    
	    vector<string> best_hits_list;
	    if(items.size() > 19){				
						string best_hits_loc = splitStrIndex(items[19],":",-1);
		    vector<string> split_loc;
		    splitStr(best_hits_loc,";",split_loc);
		    for (int k=0;k<split_loc.size()-1;k++){//extra ; at end so need the -1
							best_hits_list.push_back(splitStrIndex(split_loc[k],",",0));			
						} 
					}
					count_assigned_aligned_reads++;
					
					if (aligned_id.substr(0,4) == "ERCC"){
					if(best_hits_list.size() && !mm) continue;
						const int erccIndex=ercc_to_index[aligned_id];
						uint64_t umi=wellIndex*wellBase+barcodeIndex*barcodeBase+erccIndex+geneList.size();
						count_spike_total[wellIndex][erccIndex]++;						
						if(!umi_seen.count(umi)){
							umi_seen.insert(umi);
							count_spike_umi[wellIndex][erccIndex]++;
						}
					}
					else if 	(aligned_id.substr(0,4) == "chrM"){
						if(best_hits_list.size() && !mm) continue;	
						uint64_t umi=wellIndex*wellBase+barcodeIndex*barcodeBase+barcodeBase-2;
						count_assigned_mito_reads++;
						if(!umi_seen.count(umi)){
						 count_assigned_mito_umi++;
						 umi_seen.insert(umi);
						}						
					}
					else if (refseq_to_gene.count(aligned_id)){
						string gene = refseq_to_gene[aligned_id];
						const unsigned int geneIndex =gene_to_index[gene];

						bool multiGene=multiGeneHit(best_hits_list,gene,refseq_to_gene);
						if(!mm && multiGene) continue;
							uint64_t umi=wellIndex*wellBase+barcodeIndex*barcodeBase+geneIndex;
						count_total[wellIndex][geneIndex]++;
					 if(!umi_seen.count(umi)){
						 count_umi[wellIndex][geneIndex]++;
						 umi_seen.insert(umi);
					 }		
					}
					else{
						string umi=barcodePanel.wells[wellIndex] + barcode + aligned_id;
						if(!mm && best_hits_list.size())continue;
						count_assigned_unknown_reads++;
						if(!unknown_set.count(aligned_id)){
						 unknown_set.insert(aligned_id);
						}
					 if(!unknown_umi_seen.count(umi)){
						 count_assigned_unknown_umi++;	
							unknown_umi_seen.insert(umi);
						}			
					}		 
			 }
			}
		 fclose(fp);
		}
	
	 
	 FILE *fp=fopen(logFile.c_str(),"w");
	 if(!fp)exit(EXIT_FAILURE);
	 fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
	 fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	           sample_id,
	           count_total_reads,
												count_assigned_reads,
												count_assigned_aligned_reads,
												sumCounts(count_spike_total,NWELLS,erccList.size()),
												sumCounts(count_spike_umi,NWELLS,erccList.size()),
												count_assigned_mito_reads,
												count_assigned_mito_umi,
												sumCounts(count_total,NWELLS,geneList.size()),
												sumCounts(count_umi,NWELLS,geneList.size()),
												count_assigned_unknown_reads,
												count_assigned_unknown_umi
										);
	 fclose(fp);
	 
	 if(!(fp=fopen(unknownFile.c_str(),"w")))exit(EXIT_FAILURE);
	 for (auto itr = unknown_set.begin(); itr != unknown_set.end(); ++itr) {
	  fprintf(fp,"%s\n",itr->c_str());
	 }
	 fclose(fp);
	  
	 if(!(fp=fopen(totalAlignFile.c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 for(int i=0;i<geneList.size();i++){
	  fprintf(fp,"%s\t",geneList[i].c_str());
	  for(int j=0;j<NWELLS-1;j++){
				fprintf(fp,"%d\t",count_total[j][i]);
			}
			fprintf(fp,"%d\n",count_total[NWELLS-1][i]);
		} 
	 fclose(fp);
	 
	 if(!(fp=fopen(umiFile.c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");  
	 for(int i=0;i<geneList.size();i++){
	  fprintf(fp,"%s\t",geneList[i].c_str());
	  for(int j=0;j<NWELLS-1;j++){
				fprintf(fp,"%d\t",count_umi[j][i]);
			}
			fprintf(fp,"%d\n",count_umi[NWELLS-1][i]);
		}  
	 fclose(fp);
  if(!(fp=fopen(erccFile.c_str(),"w")))exit(EXIT_FAILURE);
  for(int i=0;i<wellList.size();i++){
	  fprintf(fp,"\t%s",wellList[i].c_str());
	 }
	 fprintf(fp,"\n"); 
	 for(int i=0;i<erccList.size();i++){
	  fprintf(fp,"%s",erccList[i].c_str());
	  for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",count_spike_umi[j][i]);
			}
			fprintf(fp,"\n");
		}
	 fclose(fp);   
	 if(!(fp=fopen(erccUMIFile.c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 for(int i=0;i<erccList.size();i++){
	  fprintf(fp,"%s",erccList[i].c_str());
	  for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",count_spike_total[j][i]);
			}
			fprintf(fp,"\n");
		}
	 fclose(fp);
	
	 if(!(fp=fopen(wellTotalFile.c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 
	 fprintf(fp,"Refseq_Total"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_total,i,geneList.size()));
		}
		fprintf(fp,"\n");
		
	 fprintf(fp,"Refseq_UMI"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_umi,i,geneList.size()));
		}
		fprintf(fp,"\n");
		
	 fprintf(fp,"Spike_Total"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_spike_total,i,erccList.size()));
		}
		fprintf(fp,"\n");
		
  fprintf(fp,"Spike_UMI"); 
  for(int i=0;i<NWELLS;i++){
	 	fprintf(fp,"\t%d",sumCountsi(count_spike_umi,i,erccList.size()));
	 }
	 fprintf(fp,"\n");
  fclose(fp);

 	 //clean up
	 delete[]count_spike_total[0];
	 delete[]count_spike_umi[0];   
	 delete[]count_total[0];
	 delete[]count_umi[0];
 }   
 return 1;		
	 
}
