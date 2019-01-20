
#Cpp based tools for barcoded RNAseq

##TL;DR Summary

Set of cpp tools for UMI based RNAseq. These are decendants of python scripts used to split fastq files before alignment and merge files after alignment by bwa. All the binaries are multithreaded unless otherwise stated.

Installation: From the source directory 
	make clean; make  (NWELLS=96/384 - default is 96) 
 
 This will result in executable binaries created in the source/w96 or source/w384 directories
 
For all binaries the -h flag gives documentation about the available flags and an example of how to use the binary: 

umisplit: takes pairs of R1 and R2 fastq files and generates a new fastq file with the UMI from the R2 file merged to the header line of the R1 file. It also separates the reads based on the plate well barcode part of the UMI. Each thread handles a different pair of input files.

umimerge_parallel: takes sam files organized by plate wells and counts the reads mapped to each gene. The original gene level filtering scheme used by the Python script (i.e. reads with the same barcode and same gene mapping are filtered) is supported by the -g option. Otherwise position level filtering is used where reads with the same barcode that map to the same position in the gene are filtered which is the scheme from the original UMI paper. Each thread handles a different sam file.


multibwa.sh: A bash script to search a directory for sam files and run bwa in parallel. The input are files split into directories corresponding to the plate well barcode. 
## How to use:

The executable binaries are set up to use the same parameters as the original python scripts. 
## Summary of implementation changes
The major changes from the Python scripts are:

1. Use of unique numerical translation of UMI sequence and mapping position for hash to determine whether the read is a duplicate. This allows for deduplication of reads that map to the same or similar position rather than reads that map to the same gene.

2. Splitting reads by wells to increase speed, decrease memory requirements and allow for parallel processing of sam files. 

3. Use of a lookup table to detemine whether an ambiguous sequence maps uniquely to a known barcode

4. Slightly faster bwa alignment using UNIX threads to parallelize at the application level rather than relying on bwa's built-in parallelization.

##Building and executing Docker container
The build.sh script creates a minimum container called rna-umi-cpp from scratch. This is done in 2 steps - a full build environment to compile the code and then a minimal environment for the runtime execs

To execute the docker container (for 96 well plates):
	sudo docker run --rm rna-umi-cpp <cmd> <arguments>
For 384 well plates:	
	sudo docker run --rm -e NWELLS=384 rna-umi-cpp <cmd> <arguments>


