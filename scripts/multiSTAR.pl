#!/usr/bin/perl
use strict; use warnings;
use threads;


# 1.3 Reference
our ($TOP_DIR,$REF_DIR,$SPECIES_DIR,$ALIGN_DIR,$BWA_ALN_SEED_LENGTH,$BWA_SAM_MAX_ALIGNS_FOR_XA_TAG,$nThreads)=@ARGV;
our $REF_SEQ_DIR="$SPECIES_DIR/STAR_chrM_ERCC_hg19_indices";

#system("mkdir -p /tmp/locks.$$");

my @wells=split(' ',`ls $ALIGN_DIR`);
foreach my $well (@wells){
	doOperation("$ALIGN_DIR/$well");
}	


sub doOperation{
	my($dir)=@_;
	my $start_string="\@SQ";
	printf stderr "working on $dir\n";
	print STDERR "cat $dir/*.fq >  $dir/tempall.fastq\n";
	system("cat $dir/*.fq >  $dir/tempall.fastq");
        print STDERR "STAR --genomeDir $REF_SEQ_DIR --runThreadN 16 --readFilesIn  $dir/tempall.fastq  --genomeLoad LoadAndKeep  --outFileNamePrefix $dir/temp\n";
	my $cmd="STAR --genomeDir $REF_SEQ_DIR --runThreadN 16 --readFilesIn  $dir/tempall.fastq  --genomeLoad LoadAndKeep  --outFileNamePrefix $dir/temp";
	system("$cmd");       
	print STDERR "cat $dir/tempAligned.out.sam | grep -v '$start_string' > $dir/combined.sam\n";
	system("cat $dir/tempAligned.out.sam | grep -v '$start_string' > $dir/combined.sam");
	#print STDERR "cp $dir/tempAligned.out.sam  $dir/combined.sam";
	#system("cp $dir/tempAligned.out.sam  $dir/combined.sam");
	system("rm -r $dir/temp*");
}
