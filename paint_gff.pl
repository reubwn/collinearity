#!/usr/bin/env perl

## author: reubwn Aug 2017

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw /pairwise/;
use Sort::Naturally;

my $usage = "
INFO
  Obsoleted by calculate_collinearity_breakpoints_GFF.pl

SYNOPSIS
  Prints block number participation for all genes based on Ks threshold in new column of GFF

OPTIONS
  -i|--in     [FILE]  : collinearity file (REFORMATTED ONLY)
  -g|--gff    [FILE]  : GFF file
  -s|--score  [FILE]  : score file with average Ks per block
  -k|--ks     [FLOAT] : Ks threshold to define homologous block (default <= 0.5)
  -o|--out            : outfile (default=GFF.painted)
  -h|--help           : print this message

USAGE
  paint_gff.pl -i Xy.collinearity.kaks.reformatted -g Xy.gff -s Xy.collinearity.kaks.score
\n";

my ($collinearityfile,$gfffile,$scorefile,$help,$verbose);
my $ks = 0.5;

GetOptions (
  'i|collinearity=s' => \$collinearityfile,
  'g|gff=s'          => \$gfffile,
  's|score=s'        => \$scorefile,
  'k|ks:f'           => \$ks,
  'v|verbose'        => \$verbose,
  'h|help'           => \$help
);

die $usage if $help;
die $usage unless ($collinearityfile && $gfffile && $scorefile);
print STDERR "[INFO] Ks threshold set to $ks\n";

my (%collinearity_hash, %homologous_blocks_hash, %score_hash, %score_file_hash, %seen);

open (my $SCORE, $scorefile) or die $!;
while (<$SCORE>) {
  chomp;
  my @F = split (/\s+/, $_);
  $score_hash{$F[0]}{SCORE} = $F[9];
  $score_hash{$F[0]}{Ka} = $F[10];
  $score_hash{$F[0]}{Ks} = $F[11];
  $score_file_hash{$F[0]} = $_; ##file hash
}
close $SCORE;

open (my $COLL, $collinearityfile) or die $!;
while (<$COLL>) {
  chomp;
  if ($_ =~ /^#/) {
    next;
  } else {
    my @F = split (/\s+/, $_);
    push ( @{$collinearity_hash{$F[2]}}, $F[0] ); ##key= gene name; val= @{all blocks that gene is a member of}
    push ( @{$collinearity_hash{$F[3]}}, $F[0] );
  }
}
close $COLL;
print STDERR "[INFO] Parsed ".scalar(keys %collinearity_hash)." genes from $collinearityfile\n";

## print blocks per gene, maybe useful:
open (my $GENES, ">$gfffile.blocks_per_gene") or die $!;
foreach my $gene (nsort keys %collinearity_hash) {
  print $GENES "$gene\t@{$collinearity_hash{$gene}}\n";
  foreach my $block (@{$collinearity_hash{$gene}}) {
    push ( @{$homologous_blocks_hash{$gene}}, $block ) if $score_hash{$block}{Ks} <= $ks; ## %homologous_blocks_hash contains ONLY homologous blocks with Ks <= threshold
  }
}
close $GENES;


open (my $OUT, ">$gfffile.painted") or die $!;
open (my $ORDER, ">$scorefile.reordered") or die $!;
print $ORDER "block_num\tchrom1\tchrom2\tcollinear_genes\ttotal_genes1\ttotal_genes2\torientation\tscore_block1\tscore_block2\tscore_avg\tka_avg\tks_avg\n";
open (my $GFF, $gfffile) or die $!;

while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  if (exists($homologous_blocks_hash{$F[1]})) {
    my @blocks = @{$homologous_blocks_hash{$F[1]}}; ## all homologous blocks for that gene; should be 1 but sometimes more
    if (scalar(@blocks)>1) {
      print STDERR "[WARN] Gene $F[1] has >1 homologous block: @blocks\n" if $verbose;
    } else {
      print $OUT join ("\t", @F, @blocks, "\n");
      print $ORDER "$score_file_hash{$blocks[0]}\n" unless exists($seen{$blocks[0]}); ##print once?
      $seen{$blocks[0]}++;
    }
  } else {
    print $OUT join ("\t", @F, "-", "\n");
  }
}
close $GFF;
close $OUT;
close $ORDER;

print STDERR "[INFO] Finished on ".`date`."\n";
