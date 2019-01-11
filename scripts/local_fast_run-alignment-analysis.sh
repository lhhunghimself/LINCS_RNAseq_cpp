#!/bin/bash

#RUN as root
# This script converts a series of mRNA sequencing data file in FASTQ format
# to a table of UMI read counts of human genes in multiple sample conditions.

#set the correct environment variable - this chooses between the binaries for 96 and 384 wells

NWELLS=96
THREAD_NUMBER=1

#DOCKER COMMANDS - comment out if you want to run with the local binaries
DPREFIX="";
DOCKER_CPP=""
DOCKER_BWA=""


# 1 Parameters

# 1.1 Global

TOP_DIR="${PWD}/data"

# 1.2 Dataset
SERIES="20150409"
SAMPLE_ID="RNAseq_${SERIES}"
LANES=6
DATA_DIR=$TOP_DIR/UMITestShort
#SEQ_DIR=$TOP_DIR/UMI/Seqs
SEQ_DIR=$TOP_DIR/shortSeqs
ALIGN_DIR="${DATA_DIR}/Aligns"
COUNT_DIR="${DATA_DIR}/Counts"

# 1.3 Reference
REF_DIR="$TOP_DIR/Broad_UMI"
SPECIES_DIR="${REF_DIR}/Human_RefSeq"
REF_SEQ_FILE="${SPECIES_DIR}/refMrna_ERCC_polyAstrip.hg19.fa"
SYM2REF_FILE="${SPECIES_DIR}/refGene.hg19.sym2ref.dat"
ERCC_SEQ_FILE="${REF_DIR}/ERCC92.fa"
BARCODE_FILE="${REF_DIR}/barcodes_trugrade_96_set4.dat"

# 1.4 Program
PROG_DIR="$TOP_DIR/Programs/R20160310"
BWA_ALN_SEED_LENGTH=24
BWA_SAM_MAX_ALIGNS_FOR_XA_TAG=20


#ARGUMENTS


# 2 Computation

# 2.1 Alignment

# Align sequence fragments to reference genome library.
let "IDX = 1"
while [ "$IDX" -le "${LANES}" ]; do
        SUBSAMPLE_ID="Lane$IDX"
#        SEQ_FILE_R1="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R1.fastq"
#        SEQ_FILE_R2="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R2.fastq"
        SEQ_FILE_R1="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R1.fastq.gz"
        SEQ_FILE_R2="${SEQ_DIR}/${SAMPLE_ID}_${SUBSAMPLE_ID}_R2.fastq.gz"
        SEQ_FILES="${SEQ_FILES} ${SEQ_FILE_R1} ${SEQ_FILE_R2}"
        let "IDX = $IDX + 1"
done
#echo "${DOCKER_CPP} umisplit  -v -l 16 -m 0 -N 0 -f -o $ALIGN_DIR -t $THREAD_NUMBER -b $BARCODE_FILE $SEQ_FILES"
#${DOCKER_CPP} umisplit  -v -l 16 -m 0 -N 0 -f -o $ALIGN_DIR -t $THREAD_NUMBER -b $BARCODE_FILE $SEQ_FILES
# 2.2 Counting
# Count the number of sequence alignments for reference genes.
#echo "$DOCKER_BWA multibwa.pl $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER"
#$DOCKER_BWA multibwa.pl $TOP_DIR $REF_DIR $SPECIES_DIR $ALIGN_DIR $BWA_ALN_SEED_LENGTH $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $THREAD_NUMBER 


echo "$DOCKER_CPP umimerge_parallel -i $SAMPLE_ID -s $SYM2REF_FILE -e $ERCC_SEQ_FILE -b $BARCODE_FILE -a $ALIGN_DIR -o $COUNT_DIR -t $THREAD_NUMBER -m 1000 -p 0 "
#$DOCKER_CPP  umimerge_parallel -i $SAMPLE_ID -s $SYM2REF_FILE -e $ERCC_SEQ_FILE -b $BARCODE_FILE -a $ALIGN_DIR -o $COUNT_DIR -t $THREAD_NUMBER -m 1000 -p 0

