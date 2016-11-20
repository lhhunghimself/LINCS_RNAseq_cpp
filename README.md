
#Cpp based tools for barcoded RNAseq

##TL;DR Summary

Set of cpp tools for UMI based RNAseq. These are decendants of python scripts used to split fastq files before alignment and merge files after alignment by bwa. All the binaries are multithreaded unless otherwise stated.

Installation: From the source directory 
	make clean; make  (NWELLS=96/384 - default is 96) 

For all binaries the -h flag gives documentation about the available flags and an example of how to use the binary: 

umisplit_orig	: takes pairs of R1 and R2 fastq files and generates a new fastq file with the UMI from the R2 file merged to the header line of the R1 file

umisplit: takes pairs of R1 and R2 fastq files and generates a new fastq file with the UMI from the R2 file merged to the header line of the R1 file. It also separates the reads based on the plate well barcode part of the UMI

umisplit_sam: takes a set of sam files produced by umisplit_orig or the original split_and_align python script and splits the lines into separate files based on the barcode of the UMI

umimerge: takes a set of sam files produced by umisplit and counts the reads mapped to each gene and plate well - single threaded

umimerge_single_pass: same as umi_merge but attempts to count multiple alignments and unique alignments in a single pass at the cost of more RAM (~20 GB required ) - single threaded

umimerge_parallel: takes sam files organized by plate wells and counts the reads mapped to each gene - multithreaded and has option of using position level filtering by UMIs

umisample: simple utility to create smaller files from subsets of lines of fastq files for testing purposes

multbwa_orig.pl and multibwa.pl: Perl scripts to search a directory for sam files and run bwa in parallel. multbwa_orig is for used with the original style of split sam files. multibwa is for use with files split into directories corresponding to the plate well barcode

##Implementation Changes

The major changes from the python scripts are:

1. Use of unique numerical translation of UMI sequence and mapping position for hash to determine whether the read is a duplicate. This allows for deduplication of reads that map to the same or similar position rather than reads that map to the same gene

2. Splitting by wells to increase speed, decrease memory requirements and allow for parallel processing of sam files.

3. umisplit_sam is provided to split by wells after alignment This allows for reanalysis of existing data without need for realignment. 

4. Use of a lookup table to detemine whether an ambiguous sequence maps uniquely to a known barcode

5. Slightly faster bwa alignment using a perl threads to parallelize at the application level rather than relying on bwa's built-in parallelization.

##Building and executing Docker container
The build.sh script creates a minimum container called rna-umi-cpp from scratch. This is done in 2 steps - a full build environment to compile the code and then a minimal environment for the runtime execs

To execute the docker container (for 96 well plates):
	sudo docker run --rm rna-umi-cpp <cmd> <arguments>
For 384 well plates:	
	sudo docker run --rm -e NWELLS=384 rna-umi-cpp <cmd> <arguments>


