#!/usr/bin/env perl

## author: reubwn June 2017

use strict;
use warnings;
use Getopt::Long;
use List::Util qw /sum/;
use String::Random;
use Sort::Naturally;

my $usage = "
SYNOPSIS

OUTPUT

OPTIONS
  -a|--bam    [FILE] : BAM file (required)
  -r|--region [STR]  : specify a region(s) to analyse (default: all)
  -b|--bed    [FILE] : windows in BED format (required)
  -g|--genome [FILE] : BEDtools style genome lengths file (required)
  -i|--insert [INT]  : insert size threshold; proportion of pairs with ISIZE >= this value reported (default: 500)
  -o|--out           : outfile name (default: STDOUT)
  -h|--help          : print this message

USAGE

\n";

my ($bam,$bed,$genome,$help,$debug);
my $insert = 500;
my $region = "";
my $n = 1;

GetOptions (
  'a|bam=s'  => \$bam,
  'r|region:s' => \$region,
  'b|bed=s'  => \$bed,
  'g|genome=s' => \$genome,
  'i|insert'  => \$insert,
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
  "SAME",#$same,
  "S:T",#($same/$total),
  "INSERT",#$insert,
  "I:T",#($insert/$total),
  "INSERT_AVG",#(sum(@insert_arr)/scalar(@insert_arr)),
  "\n"
);

open (my $BED, $bed) or die "[ERROR] Cannot open $bed: $!\n";
while (my $window = <$BED>) {
  #print STDERR "\r[INFO] Working on window \#$n: $_";$| = 1;
  #`printf "$window" > tmp.bed`;
  my $rand = new String::Random;
  my $tmp = $rand->randregex('\w\w\w\w\w\w\w').".bed";
  open (my $TMP, ">$tmp");
  print $TMP $window;
  close $TMP;
  my ($total,$same,$insert,$insert_avg) = (0,0,0,0);
  my @insert_arr;
  open(my $SAM, "samtools view -b $bam $region | bedtools intersect -sorted -g $genome -a stdin -b $tmp | samtools view - |");#`bedtools intersect -sorted -g $genome -a $bam -b <(printf "$_") | samtools view - | perl -lane 'if($F[6]eq"="){if($F[8]>500){$insert++};$same++;$total++}else{$total++}END{print "$total\t$same\t".($same/$total)."\t$insert\t".($insert/$total)}'`;
  while (<$SAM>) {
    my @F = split (/\s+/, $_);
    if ($F[6] eq "=") { ##mate on same scaffold
      $same++;
      $total++;
      $insert++ if ($F[8] > $insert);
      push (@insert_arr, abs($F[8]));
    } else {
      $total++;
    }
  }
  close $SAM;
  chomp($window);
  print STDOUT join (
    "\t",
    $window,
    $total,
    $same,
    ($same/$total),
    $insert,
    ($insert/$total),
    (sum(@insert_arr)/scalar(@insert_arr)),
    "\n"
  );
}

close $BED;
print STDERR "\n[INFO] Finished on ".`date`."\n";
