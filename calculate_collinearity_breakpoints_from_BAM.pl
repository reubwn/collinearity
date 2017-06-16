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
my $insert = 500;
my $region = "";
my %stats_hash;
my @inserts_dryrun;
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
print STDERR "[INFO] Region specified: $region\n" if ($region =~ /.+/);
print STDERR "[INFO] Windows file: $bed\n";
print STDERR "[INFO] Genome file: $genome\n";
print STDERR "[INFO] Dry run, will only calculate ISIZE\n" if ($dryrun);
print STDOUT join (
  "\t",
  "CHROM\tSTART\tEND",#$window,
  "TOTAL",#$total,
  "DIFF",#$same,
  "PROP_DIFF",#($same/$total),
  "BIGINSERT",#$insert,
  "PROP_BIGINSERT",#($insert/$total),
  "INSERT_AVG",#(sum(@insert_arr)/scalar(@insert_arr)),
  "\n"
) unless ($dryrun);

my $stat = Statistics::Descriptive::Full->new();
my $stat_dryrun = Statistics::Descriptive::Full->new();

open (my $BED, $bed) or die "[ERROR] Cannot open $bed: $!\n";
while (my $line = <$BED>) {
  next unless $line =~ /^$region/;
  print STDERR "\r[INFO] Working on window \#$n";$| = 1;
  chomp($line);
  my @window = split (/\s+/, $line);

  ## print single bed entry to tmp.bed
  # open (my $TMP, ">tmp.bed");
  # print $TMP "$window\n";
  # close $TMP;

  my @isize;
  my ($total,$same,$insert,$insert_avg) = (0,0,0,0);
  #open(my $SAM, "samtools view -F1536 -b $bam $region | bedtools intersect -sorted -g $genome -a stdin -b tmp.bed | samtools view - |");
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
        $total++;
        $insert++ if ($F[8] > $insert);
        push (@isize, abs($F[8]));
      } else {
        $total++;
      }
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
      $line,
      $total,
      ($total - $same),
      (($total - $same)/$total),
      $insert,
      ($insert/$total),
      ($stat->mean()),
      ($stat->median())
      "\n"
    );
  }
  $n++;
}
close $BED;

$stat_dryrun->add_data(@isize_dryrun);
print STDERR "\n\n[INFO] ISIZE mean: ".$stat->mean();
print STDERR "\n[INFO] ISIZE median: ".$stat->median();
print STDERR "\n[INFO] ISIZE 5\% and 95\%: ".$stat->percentile(5).", ".$stat->percentile(95);

print STDERR "\n[INFO] Finished on ".`date`."\n";
