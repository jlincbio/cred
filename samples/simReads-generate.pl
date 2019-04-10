#!/usr/bin/env perl
# sample script to generate simulated IonTorrent data

use strict;
use warnings;

sub which {
	my $cmd = shift;
	my $path = `which $cmd`; chomp $path;
	return $path;
}

my $n_miss = 0; my %cmds;
foreach my $cmd (<DATA>) {
	chomp $cmd;
	my $p = which($cmd);
	if ($p eq '') {
		$n_miss++;
		printf STDERR "ERROR: Cannot locate $cmd!\n";
	} else {
		$cmds{$cmd} = $p;
	}
}
die "ERROR: Missing executables found!\n" if ($n_miss > 0);

my $ref_fasta = shift @ARGV or die "ERROR: No reference fasta specified!\n";
my $ref_slop = 'simReads_refsites.bed';
my $ctrl_size = 50000000;
my $pos_size = 3000000;

my $fastq_base  = "simReads_hg19_control";
my $fastq_site  = "simReads_hg19_site";
my $fastq_pdown = "simReads_hg19_treatment";

printf STDERR "Press any key to continue..."; <STDIN>;
printf STDERR "Mandatory cool-off period..."; sleep 1;
printf STDERR "\n";

printf STDERR "Generating reads...\n";
system "$cmds{dwgsim} -N $ctrl_size -c 2 -1 149 -2 0 -e 0 -E 0 -f TACGTACGTCTGAGCATCGATCGATGTACAGC -r 0 $ref_fasta\.fa $fastq_base";
system "$cmds{dwgsim} -N $pos_size -c 2 -1 149 -2 0 -r 0 -e 0 -E 0 -f TACGTACGTCTGAGCATCGATCGATGTACAGC -x $ref_slop $ref_fasta\.fa $fastq_site"; 

printf STDERR "Making pulldown fastq...\n";
system "cat $fastq_base\.bwa.read1.fastq >> $fastq_pdown\.fastq";
system "cat $fastq_site\.bwa.read1.fastq >> $fastq_pdown\.fastq";
system "mv $fastq_base\.bwa.read1.fastq $fastq_base\.fastq";

printf STDERR "Aligning control reads by LAST...\n";
system "$cmds{lastal} -Q1 -D100 $ref_fasta $fastq_base\.fastq | $cmds{last-split} > $fastq_base\.maf";
system "$cmds{maf-convert} sam $fastq_base\.maf > $fastq_base\.sam";
system "$cmds{gzip} $fastq_base\.maf";
system "$cmds{samtools} view -@ 4 -b -t $ref_fasta\.fa.fai -o $fastq_base\.bam $fastq_base\.sam";
system "$cmds{samtools} sort -@ 4 -o $fastq_base\-sorted.bam $fastq_base\.bam";
system "$cmds{samtools} index $fastq_base\-sorted.bam";

printf STDERR "Aligning pulldown reads by LAST...\n";
system "$cmds{lastal} -Q1 -D100 $ref_fasta $fastq_pdown\.fastq | $cmds{last-split} > $fastq_pdown\.maf";
system "$cmds{maf-convert} sam $fastq_pdown\.maf > $fastq_pdown\.sam";
system "$cmds{gzip} $fastq_pdown\.maf";
system "$cmds{samtools} sort -@ 4 -o $fastq_pdown\-sorted.bam $fastq_pdown\.bam";
system "$cmds{samtools} index $fastq_pdown\-sorted.bam";

my @temps = (
	"$fastq_base\.sam",
	"$fastq_base\.bam",
	"$fastq_pdown\.sam",
	"$fastq_pdown\.bam"
);
unlink $_ for @temps;

__DATA__
dwgsim
lastal
last-split
maf-convert
gzip
samtools