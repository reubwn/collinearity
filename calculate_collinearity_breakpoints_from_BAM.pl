#!/usr/bin/env perl

## author: reubwn June 2017

use strict;
use warnings;
use Getopt::Long;
use List::Util qw /sum/;
use Statistics::Descriptive;
use Sort::Naturally;

my $usage = "
SYNOPSIS
  Calculates some statistics for sliding windows across a BAM file which can inform the presence of collinearity breakpoints

OUTPUT
  A file of stats per window, including total reads, average coverage, \%reads mapping to different contig, average insert,
  \%reads with insert > specified threshold, \%reads split...

OPTIONS
  -a|--bam    [FILE] : BAM file (required)
  -r|--region [STR]  : specify a region(s) to analyse (default: all)
  -b|--bed    [FILE] : windows in BED format (required)
  -g|--genome [FILE] : BEDtools style genome lengths file (required)
  -i|--insert [INT]  : insert size threshold; proportion of pairs with ISIZE >= this value reported (default: 500)
  -d|--dryrun        : do a first pass to estimate 5% & 95% quantiles for ISIZE distribution, to inform -i parameter
  -o|--out           : outfile name (default: STDOUT)
  -h|--help          : print this message

USAGE

\n";

my ($bam,$bed,$genome,$dryrun,$help);
my $insert_threshold = 500;
my $region = "";
my %stats_hash;
my @isize_dryrun;
my $n = 1;

GetOptions (
  'a|bam=s'  => \$bam,
  'r|region:s' => \$region,
  'b|bed=s'  => \$bed,
  'g|genome=s' => \$genome,
  'i|insert:i'  => \$insert,
  'd|dryrun' => \$dryrun,
  'h|help'   => \$help
  );

die $usage if $help;
die $usage unless ($bam && $bed && $genome);

print STDERR "[INFO] BAM file: $bam\n";
print STDERR "[INFO] Genome file: $genome\n";
print STDERR "[INFO] Windows file: $bed\n";
print STDERR "[INFO] ISIZE threshold: $insert\n";
print STDERR "[INFO] Region specified: $region\n" if ($region =~ /.+/);
print STDERR "[INFO] Dry run, will only calculate ISIZE\n" if ($dryrun);
print STDOUT join (
  "\t",
  "CHROM\tSTART\tEND", #window,
  "TOTAL_READS",       #total reads,
  "DEPTH_MEAN",        #average depth (numreads/windowsize)
  "MAPPED_DIFFSCAFF",  #reads mapped to a different scaffold
  "MAPPED_DIFFSCAFFp", #proportion
  "MATE_UNMAPPED",     #reads with unmapped mate $F[1]&8 bitset
  "MATE_UNMAPPEDp",
  "BIGINSERT",         #reads with insert > threshold
  "BIGINSERTp",
  "SPLIT",             #reads with split mapping =~ SA:Z tag
  "SPLITp",
  "INSERT_MEAN",       #mean insert over window
  "INSERT_MEDIAN",     #median insert over window
  "\n"
) unless ($dryrun);

open (my $BED, $bed) or die "[ERROR] Cannot open $bed: $!\n";
while (my $line = <$BED>) {
  next unless $line =~ /^$region/;
  print STDERR "\r[INFO] Working on window \#$n" if ($n % 10 == 0);$| = 1;
  chomp($line);

  my @isize;
  my @window = split (/\s+/, $line); ##split window to get coords
  my ($total,$same,$insert,$insert_avg,$split_reads,$mate_unmapped) = (0,0,0,0,0,0);

  ## select window from SAM
  open(my $SAM, "samtools view -F1536 $bam $window[0]:$window[1]-$window[2] |");
  while (<$SAM>) {
    my @F = split (/\s+/, $_);
    if ($dryrun) {
      if ($F[6] eq "=") { ##mate on same scaffold
        push (@isize_dryrun, abs($F[8]));
      }
    } else {
      if ($F[6] eq "=") { ##mate on same scaffold
        $same++;
        $insert++ if (abs($F[8]) > $insert_threshold); ##just look at abs insert, ie ignore read orientation... good idea?? dunno
        push (@isize, abs($F[8]));
      }
      $total++; ##total reads in window
      $mate_unmapped++ if ($F[1]&8); ##bitset for mate is unmapped
      $split_reads++ if ($_ =~ m/SA:Z/); ## BWA mem flag for split read
    }
  }
  close $SAM;
  unless ($dryrun) {
    # $stats_hash{
    #   'window' => $line,
    #   'total'  => $total,
    #   'diff'   => ($total - $same),
    #   'prop_diff' => (($total - $same)/$total),
    #   'biginsert' => $insert,
    #   'prop_biginsert' => ($insert/$total),
    #   'insert_avg' => (sum(@insert_arr)/scalar(@insert_arr))
    # };
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@isize);
    print STDOUT join (
      "\t",
      $line,                              #window
      $total,                             #total reads
      ($total/($window[2]-$window[1])),   #depth average
      ($total - $same),                   #reads mapped to diff scaff
      (($total - $same)/$total),          #proportion
      $unmapped,
      ($unmapped/$total),
      $insert,                            #reads with insert > threshold
      ($insert/$total),                   #proporiotn reads insert > threshold
      $split_reads,
      ($split_reads/$total),
      ($stat->mean()),
      ($stat->median()),
      "\n"
    );
  }
  $n++;
}
close $BED;
print STDERR "\n";

my $stat_dryrun = Statistics::Descriptive::Full->new();
$stat_dryrun->add_data(@isize_dryrun);
print STDERR "\n[INFO] ISIZE mean: ".$stat_dryrun->mean();
print STDERR "\n[INFO] ISIZE median: ".$stat_dryrun->median();
print STDERR "\n[INFO] ISIZE 5\% and 95\%: ".$stat_dryrun->percentile(5).", ".$stat_dryrun->percentile(95);

print STDERR "\n[INFO] Finished on ".`date`."\n";
