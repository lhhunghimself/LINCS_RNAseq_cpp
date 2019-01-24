
#Cpp based tools for barcoded RNAseq

This a set of cpp binaries and scripts for (Unique Molecular Identifier) UMI based RNAseq. The pre-print is available [here.](https://www.biorxiv.org/content/10.1101/345819v2)  These are decendants of python scripts describe  [here.](https://www.biorxiv.org/content/early/2014/03/05/003236)  (https://www.biorxiv.org/content/early/2014/03/05/003236) used to split fastq files before alignment and merge files after alignment by bwa. All the binaries are multithreaded unless otherwise stated.

##To compile and run
### Clone repo
	git clone https://github.com
### Compile executable
	cd LINCS_RNAseq_cpp/source
	make clean; make all96; make all384
This will result in executable binaries created in the LINCS_RNAseq_cpp/source/w96 and LINCS_RNAseq_cpp/source/w384 directories
### Run example
	cd LINCS_RNAseq_cpp
	./runExample.sh
	
This script download all the necessary support files (reference indices and barcode files) and a small dataset. The transcript counts will be in the LINCS_RNAseq_cpp/data/LINCS/Counts directory

For all binaries the -h flag gives documentation about the available flags and an example of how to use the binary: 

##To run using Docker

As an alternative to compiling and running the executables, we also provide a docker container with the executables installed.

### Building Docker containers locally

As an alternative to pulling the Docker container from our Biodepot repository, the user can also build the containers locally using the provided Docker files. 

	cd  LINCS_RNAseq_cpp
	scripts/build.sh

If the container is not built, it will be automatically downloaded when the demo script is run.

### Running the example using Docker

	cd LINCS_RNAseq_cpp
	./runDockerExample.sh

### How to adapt the software for user data


### Documentation for individual executables in the suite
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



