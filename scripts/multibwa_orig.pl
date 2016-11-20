#!/usr/bin/perl
use strict; use warnings;
use threads;


# 1.3 Reference
our ($TOP_DIR,$REF_DIR,$SPECIES_DIR,$ALIGN_DIR,$BWA_ALN_SEED_LENGTH,$BWA_SAM_MAX_ALIGNS_FOR_XA_TAG,$nThreads)=@ARGV;

our $REF_SEQ_FILE="$SPECIES_DIR/refMrna_ERCC_polyAstrip.hg19.fa";

system("mkdir -p /tmp/locks.$$");
my @threads = initThreads($nThreads);
print "Using threads @threads\n";
our @cmds,
our @done;

my @files=split(' ',`find $ALIGN_DIR -type f`);
foreach my $file (@files){
	if(substr($file,-4) eq ".sam"){system("rm $file");}
	else{
  push(@cmds,"$file");
	}	
}
foreach(@threads){
	$_ = threads->create(\&doOperation);
}

foreach(@threads){
 $_->join();
}
system("rm -rf /tmp/locks*");
				

sub doOperation{
	my $id = threads->tid();
	foreach my $i (0..$#cmds){
		my $lockfile="/tmp/locks.$$/lock.$i";
 	my $donefile="/tmp/locks.$$/done.$i";
		if (mkdir($lockfile)){
	  my $file=$cmds[$i];
	  my $cmd="(bwa aln -l $BWA_ALN_SEED_LENGTH -t 1 $REF_SEQ_FILE $file | bwa samse -n $BWA_SAM_MAX_ALIGNS_FOR_XA_TAG $REF_SEQ_FILE - $file > $file.sam) &> /dev/null";
	  printf stderr "thread %d: %s\n",$id,$cmd;
	  system("$cmd");
	  mkdir($donefile);
		}
		else{
			#printf stderr "%d %s $i exists\n",$id,$lockfile;
		}	
	}
	threads->exit();
}
sub initThreads{
	my($numThreads)=@_;
	my @initThreads;
	for(my $i = 1;$i<=$numThreads;$i++){
		push(@initThreads,$i);
	}
	return @initThreads;
}
