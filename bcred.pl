#!/usr/bin/env perl
# BCRED: a companion batch helper to CRED (Apr 8, 2019)

use strict;
use warnings;
use Parallel::ForkManager;
use Scalar::Util 'looks_like_number';
use Getopt::Std;
use Cwd;
use File::Basename;
use File::Temp 'tempdir';

use constant EXIT_SUCCESS => 0;
use constant EXIT_GENERAL => 1;
sub numeric {
	my ($n, $t) = @_;
	if (defined $n) {
		if (lc $t eq 'int') {
			return(looks_like_number(int $n));
		} elsif (lc $t eq 'pdec') {
			$n = eval($n) if ($n =~ m/^\d+\/\d+$/);
			return(($n > 0) && ($n <= 1) && (looks_like_number $n));
		}
		return(looks_like_number($n));
	}
	return 0;
}
sub usage {
	my ($msg, $s) = @_;
	$s = EXIT_GENERAL unless (defined $s);
	my @u = (
	'bcred: a batch assistant for CRED Ver. 0.1 (Apr 2019 Initial Release)',
	'Usage: bcred [options] -t TREATMENT.BAM -c CONTROL.BAM > OUTPUT.BED',
	'Required :',
	'  -t  TREATMENT.BAM  Path to the treatment ("pulldown") track [BAM]',
	'  -c  CONTROL.BAM    Path to the control ("input") Chem-seq track [BAM]',
	'Optional :',
	'  -p  [P-VALUE]      Significance level [default 0.0001]',
	'  -q  [SCORE]        Minimum MAPQ quality for reads to count [default 30]',
	'  -w  [INTEGER]      Size of differential windows [default 1200 bp]',
	'  -n  [INTEGER]      Number of threads to utilize [default 1]',
	'  -k                 Evaluate site significance with Kolmogorov-Smirnov',
	'Reminders:',
	'   1. BAM files must be sorted.',
	'   2. Use a pipe (">") to capture CRED output.');
	printf STDERR "$_\n" for @u;
	printf STDERR "\# $msg\n";
	exit $s;
}
sub which {
	my ($x, $rp) = @_;
	my @path = @{$rp}; my $y;
	while (not defined $y) {
		my $p = shift @path;
		last unless (defined $p);
		my $cmd = $p.'/'.$x;
		$y = $cmd if (-e $cmd);
	}
	return $y;
}
sub proceed {
	# s = 2: output mode
	my ($cmd, $s) = @_;
	my $status; my @output;
	if (defined $s) {
		if ($s == 2) {
			@output = `$cmd`;
			chomp @output;
			$status = \@output;
		} else {
			printf STDERR "$cmd\n";
			$status = -1;
		}
	} else {
		$status = system $cmd;
	}
	return $status;
}

my $pwd = cwd(); my $twd = tempdir(CLEANUP => 1);
my %opts; my @path = ($pwd, split(':', $ENV{PATH}));
getopts('t:c:m:w:p:q:n:x:kd', \%opts);
usage('Treatment BAM not specified!') unless (defined $opts{t});
usage('Control BAM not specified!') unless (defined $opts{c});
$opts{w} = 1200 unless (numeric($opts{w}, 'int'));
$opts{n} = 0 unless (numeric($opts{n}, 'int'));
$opts{p} = 0.0001 unless (numeric($opts{p}, 'pdec'));
$opts{q} = 30 unless (numeric($opts{q}, 'int'));
# -x{num} and -d are for debugging purpose and thus left undocumented
my $samtools = which('samtools', \@path);
my $cred = which('cred', \@path);

die "# [BCRED] Treatment BAM cannot be located!\n" unless (-e $opts{t});
die "# [BCRED] Control BAM cannot be located!\n" unless (-e $opts{c});
die "# [BCRED] SAMtools cannot be located! Check \$PATH!\n" unless (defined $samtools);
die "# [BCRED] CRED cannot be located! Check \$PATH!\n" unless (defined $cred);

# process contig names and prepare for files
my @header = `samtools view -H $opts{t}`; chomp @header;
my (%treatment, %control, %contig);
my $n_stop = 1; my $n_sort = 0;
foreach my $h (@header) {
	if ($h =~ m/^\@SQ/) {
		my @x = split '\t', $h;
		(my $chr = $x[1]) =~ s/^SN\://;
		(my $tchr = basename($opts{t})) =~ s/\.[bB][aA][mM]$/_$chr\.bam/;
		(my $cchr = basename($opts{c})) =~ s/\.[bB][aA][mM]$/_$chr\.bam/;
		$treatment{$chr} = $twd.'/'.$tchr;
		$control{$chr} = $twd.'/'.$cchr;
		$contig{$n_sort} = $chr;
		$n_sort++;
		$n_stop++;
	}
	if (defined $opts{x}) {
		last if ($n_stop > $opts{x});
	}
}

my @cstats; my @output; my @errors;
my $pm = new Parallel::ForkManager($opts{n});
$pm -> run_on_finish (sub {
	my ($pid, $err, $ident, $sig, $dump, $rstats) = @_;
	if (defined($rstats)) {
		printf STDERR "# [BCRED] Putting together output...";
		push @output, @{$rstats};
		printf STDERR "\n";
	} else {
		push @errors, $pid;
	}
});

printf STDERR "# [BCRED] TEMPDIR: $twd\n";
if ($opts{n} > 0) {
	printf STDERR "# [BCRED] Multithread processing with %u cores\n", $opts{n};
}
# main processing loop
foreach my $j (sort {$a <=> $b} (keys %contig)) {
	my $k = $contig{$j};
	$pm->start and next;
	my $ks = ''; $ks = '-k' if (defined $opts{k});
	printf STDERR "# [BCRED] Processing contig $k...\n";
	my %comments = (
		1 => "# [BCRED] Parsing treatment contig $k...\n",
		2 => "# [BCRED] Parsing control contig $k...\n",
		3 => "# [BCRED] Differential calling on contig $k...\n"
	);
	my %cmd = (
		1 => "$samtools view -b -o $treatment{$k} $opts{t} $k",
		2 => "$samtools view -b -o $control{$k} $opts{c} $k",
		3 => "$cred -w $opts{w} -q $opts{q} -p $opts{p} $ks -t $treatment{$k} -c $control{$k}"
	);
	foreach my $i (sort keys %cmd) {
		my $debug = $opts{d};
		if ($i == 2) {
			$debug = 2 unless (defined $opts{d});
		}
		printf STDERR $comments{$i};
		push @cstats, proceed($cmd{$i}, $debug);
	}
	unlink $control{$k};
	unlink $treatment{$k};
	$pm->finish($k, \@cstats);
}
$pm->wait_all_children;

my @combined;
foreach my $r (@output) {
	if (ref($r) eq 'ARRAY') {
		push @combined, @{$r};
	}
}

if (@errors) {
	printf STDERR "# [BCRED] Errors encountered in the following contigs:\n";
	my $chrom_err = join ' ', (sort @errors);
	printf STDERR "# $chrom_err\n";
}
exit (scalar @errors);
