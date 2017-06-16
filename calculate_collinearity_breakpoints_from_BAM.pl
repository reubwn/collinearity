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

my ($bam,$bed,$genome,$dryrun,$help,$debug);
my $insert = 500;
my $region = "";
my @inserts_dryrun;
my $n = 1;

GetOptions (
  'a|bam=s'  => \$bam,
  'r|region:s' => \$region,
  'b|bed=s'  => \$bed,
  'g|genome=s' => \$genome,
  'i|insert:i'  => \$insert,
  'd|dryrun' => \$dryrun,
  'h|help'   => \$help,
  'd|debug'  => \$debug
);

die $usage if $help;
die $usage unless ($bam && $bed && $genome);

print STDERR "[INFO] BAM file: $bam\n";
print STDERR "[INFO] Region specified: $region\n" if ($region =~ /.+/);
print STDERR "[INFO] Windows file: $bed\n";
print STDERR "[INFO] Genome file: $genome\n";
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
);

open (my $BED, $bed) or die "[ERROR] Cannot open $bed: $!\n";
while (my $window = <$BED>) {
  chomp($window);
  print STDERR "\r[INFO] Working on window \#$n: $window";$| = 1;

  ## print single bed entry to tmp.bed
  open (my $TMP, ">tmp.bed");
  print $TMP "$window\n";
  close $TMP;

  my ($total,$same,$insert,$insert_avg) = (0,0,0,0);
  my @insert_arr;
  open(my $SAM, "samtools view -F1536 -b $bam $region | bedtools intersect -sorted -g $genome -a stdin -b tmp.bed | samtools view - |");#`bedtools intersect -sorted -g $genome -a $bam -b <(printf "$_") | samtools view - | perl -lane 'if($F[6]eq"="){if($F[8]>500){$insert++};$same++;$total++}else{$total++}END{print "$total\t$same\t".($same/$total)."\t$insert\t".($insert/$total)}'`;
  while (<$SAM>) {
    my @F = split (/\s+/, $_);
    if ($dryrun) {
      if ($F[6] eq "=") { ##mate on same scaffold
        push (@inserts_dryrun, abs($F[8]));
      }
    } else {
      if ($F[6] eq "=") { ##mate on same scaffold
        $same++;
        $total++;
        $insert++ if ($F[8] > $insert);
        push (@insert_arr, abs($F[8]));
      } else {
        $total++;
      }
    }
  }
  close $SAM;
  unless ($dryrun) {
    print STDOUT join (
      "\t",
      $window,
      $total,
      ($total - $same),
      (($total - $same)/$total),
      $insert,
      ($insert/$total),
      (sum(@insert_arr)/scalar(@insert_arr)),
      "\n"
    );
  }
}
close $BED;

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@inserts_dryrun);
print STDERR "\n[INFO] ISIZE mean: ".$stat->mean()."\n";
print STDERR "\n[INFO] ISIZE median: ".$stat->mean()."\n";
print STDERR "\n[INFO] ISIZE 5\% and 95\%: ".$stat->percentile(5).", ".$stat->percentile(95)."\n";

print STDERR "\n[INFO] Finished on ".`date`."\n";
