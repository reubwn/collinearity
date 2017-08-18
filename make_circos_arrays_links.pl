#!/usr/bin/env perl

## author: reubwn Aug 2017

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw /pairwise/;
use Sort::Naturally;

my $usage = "
SYNOPSIS

OPTIONS
  -i|--in     [FILE]  : arrays file
  -g|--gff    [FILE]  : GFF file
  -o|--out            : outfile (default=INFILE.circos.links)
  -h|--help           : print this message

USAGE
  make_circos_links.pl -i Xyz.collinearity.kaks.arrays -g Xzy.gff
\n";

my ($arraysfile,$gfffile,$help);
my $ks = 0.5;

GetOptions (
  'i|collinearity=s' => \$arraysfile,
  'g|gff=s'          => \$gfffile,
  'h|help'           => \$help
);

die $usage if $help;
die $usage unless ($arraysfile && $gfffile);

my (%gff_hash);

open (my $GFF, $gfffile) or die $!;
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  $gff_hash{$F[1]}{CHROM} = $F[0];
  $gff_hash{$F[1]}{START} = $F[2];
  $gff_hash{$F[1]}{END} = $F[3];
}
close $GFF;
print STDERR "[INFO] Parsed ".scalar(keys %gff_hash)." genes from $gfffile\n";

open (my $ARR, $arraysfile) or die $!;
open (my $OUT, ">$arraysfile.circos.links") or die $!;
while (<$ARR>) {
  chomp;
  my @F = split (/\s+/, $_);
  my $color;
  if ($F[6] eq "T") { ##tandem
    $color = "252,187,161";
  } else {
    $color = "165,15,21";
  }
  print $OUT join (
    " ",
    $gff_hash{$F[2]}{CHROM},
    $gff_hash{$F[2]}{START},
    $gff_hash{$F[2]}{END},
    $gff_hash{$F[4]}{CHROM},
    $gff_hash{$F[4]}{START},
    $gff_hash{$F[4]}{END},
    "color=$color,thickness=5p",
    "\n"
  );
}
close $ARR;
close $OUT;

print STDERR "[INFO] Finished on ".`date`."\n";
