#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;
use Getopt::Long;
use List::Util qw /sum/;

my $usage = "
SYNOPSIS
  Calculates 'collinearity' score based on the number of collinear genes divided by the total number of genes within that defined block.
  Takes the collinearity file and the 'gff' file used in MCScanX analyses.

  If ka/ks values are present, eg by running the MCScanX 'add_kaks_to_MCScanX.pl' program first, script will also print average Ka and Ks values per block if -k option is set.

OUTPUT
  Prints to a file 'Xyz.collinearity.score'; prints score for each block plus an average.
  Also prints a file 'Xyz.collinearity.reformatted', which removes some of the formatting issues in the original MCScanX 'Xyz.collinearity' file.

OPTIONS
  -i|--in    [FILE] : collinearity file
  -g|--gff   [FILE] : GFF file
  -k|--kaks         : also calculate average Ka & Ks per block
  -h|--help         : print this message

USAGE
  >> calculate_collinarity_metric.pl -i xyz.collinearity -g xyz.gff
  >> calculate_collinarity_metric.pl -i xyz.collinearity.kaks -g xyz.gff -k
\n";

my ($collinearity, $gff, $kaks, $help, $debug);

GetOptions (
  'i|in=s'  => \$collinearity,
  'g|gff=s' => \$gff,
  'k|kaks'  => \$kaks,
  'h|help'  => \$help,
  'd|debug' => \$debug
);

die $usage if $help;
die $usage unless ($collinearity && $gff);

print STDERR "[INFO] Collinearity file: $collinearity\n";
print STDERR "[INFO] Add average Ka and Ks: TRUE\n" if ($kaks);

my (%blocks,$chrom1,$chrom2,$orientation);
open (my $COL, $collinearity) or die $!;
open (my $REFORMAT, ">$collinearity.refomatted") or die $!;
while (<$COL>) {
  chomp;
  if ($_ =~ m/^#/) {
    print $REFORMAT "$_\n";
    if ($_ =~ m/\s(\w+\d+)\&(\w+\w\d+)\s(plus|minus)$/) { ## get chrom names and strand orientation of block 2
      $chrom1 = $1;
      $chrom2 = $2;
      $orientation = $3;
      next;
    } else {
      next;
    }
  }
  $_ =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces

  my @F = split (m/\s+/, $_);
  my $aln_number;
  if ($F[0]=~m/\d+\-\d+\:/) { ## sometimes columns not formatted properly... :/
    my @a = split (m/\-/, $F[0]);
    push @{ $blocks{$a[0]}{block1} }, $F[1];
    push @{ $blocks{$a[0]}{block2} }, $F[2];
    $aln_number = $a[0];
    ## print to $REFORMAT; this should be easier to use for downstream analyses
    $a[1] =~ s/\://;
    my @N = @F;
    print $REFORMAT join "\t", @a, splice (@N,1), "\n";
  } else {
    $F[0] =~ s/\-//;
    push @{ $blocks{$F[0]}{block1} }, $F[2];
    push @{ $blocks{$F[0]}{block2} }, $F[3];
    $aln_number = $F[0];
    $F[1] =~ s/\://;
    print $REFORMAT join "\t", @F, "\n";
  }

  ## dump genes and plus/minus info into %blocks
  $blocks{$aln_number}{chrom1} = $chrom1;
  $blocks{$aln_number}{chrom2} = $chrom2;
  $blocks{$aln_number}{orientation} = $orientation;
  print STDOUT "$aln_number $chrom1 $chrom2 $orientation\n" if $debug;
  if ($kaks) {
    push @{ $blocks{$aln_number}{ks} }, $F[-1]; ## ks is in final column
    push @{ $blocks{$aln_number}{ka} }, $F[-2]; ## ka is in second to last column
  }
}
close $COL;
close $REFORMAT;

open (my $OUT, ">$collinearity.score") or die $!;
if ($kaks) {
  print $OUT join "\t", "block_num","chrom1","chrom2","collinear_genes","total_genes1","total_genes2","orientation","score_block1","score_block2","score_avg","ka_avg","ks_avg","\n";
} else {
  print $OUT join "\t", "block_num","chrom1","chrom2","collinear_genes","total_genes1","total_genes2","orientation","score_block1","score_block2","score_avg","\n";
}

foreach (sort {$a<=>$b} keys %blocks) {

  ## get orientation of block2
  my $orientation = $blocks{$_}{orientation};

  ## get genes of block1
  my @block1_genes = @{ $blocks{$_}{block1} };
  my $bl1_start = shift @block1_genes;
  my $bl1_end   = pop @block1_genes;
  my $bl1_length = `perl -e 'while (<>){print if (/\t\Q$bl1_start\E\t/../\t\Q$bl1_end\E\t/);}' $gff | wc -l`;
  chomp ($bl1_length);
  my $score_block1 = sprintf("%.5f",(scalar(@block1_genes)/$bl1_length));

  ## get genes of block2
  my @block2_genes = @{ $blocks{$_}{block2} };
  my $bl2_start = shift @block2_genes;
  my $bl2_end   = pop @block2_genes;
  my $bl2_length;
  ## if block 2 is in minus orientation, need to reverse the search!
  if ($orientation eq "plus") {
    $bl2_length = `perl -e 'while (<>){print if (/\t\Q$bl2_start\E\t/../\t\Q$bl2_end\E\t/);}' $gff | wc -l`;
  } elsif ($orientation eq "minus") {
    $bl2_length = `perl -e 'while (<>){print if (/\t\Q$bl2_end\E\t/../\t\Q$bl2_start\E\t/);}' $gff | wc -l`;
  } else {
    die "\nUnknown strand orientation for block 2: $orientation\n\n";
  }
  chomp ($bl2_length);
  print STDERR "\r[INFO] Calculating scores for block: $_";
  my $score_block2 = sprintf("%.5f",(scalar(@block2_genes)/$bl2_length));

  ## get kaks values if present
  if ($kaks) {
    my @ka = grep {$_ >= 0} @{ $blocks{$_}{ka} }; ## exclude negative values from calculation;
    my @ks = grep {$_ >= 0} @{ $blocks{$_}{ks} }; ## these are "-2", output when ka/ks cannot be calulated for some reason

    print $OUT join "\t",
      $_,
      $blocks{$_}{chrom1},
      $blocks{$_}{chrom2},
      scalar(@block1_genes),
      $bl1_length,
      $bl2_length,
      $orientation,
      $score_block1,
      $score_block2,
      sprintf("%.5f",(($score_block1+$score_block2)/2)),
      sprintf("%.5f",(avg(@ka))),
      sprintf("%.5f",(avg(@ks))),
      "\n";
  } else {
    print $OUT join "\t",
      $_,
      scalar(@block1_genes),
      $bl1_length,
      $bl2_length,
      $orientation,
      $score_block1,
      $score_block2,
      sprintf("%.5f",(($score_block1+$score_block2)/2)),
      "\n";
  }
}
close $OUT;
print STDERR "\n[INFO] Finished on ".`date`."\n";

sub avg {
  return sum(@_)/@_;
}

__END__
