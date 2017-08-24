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
  -i|--in     [FILE]  : collinearity file (REFORMATTED ONLY)
  -g|--gff    [FILE]  : GFF file
  -b|--breaks [FILE]  : breaks file
  -r|--bcolor [STRING]: colour to use for detected breaks links
  -s|--score  [FILE]  : score file with average Ks per block
  -k|--ks     [FLOAT] : Ks threshold; don't print links with Ks > this value
  -c|--chrom  [STRING]: colour all links from chrom with this name...
  -l|--color  [STRING]: ... this colour (see etc/colors.conf for legit colnames)
  -o|--out            : outfile (default=INFILE.circos.links)
  -h|--help           : print this message

USAGE
  make_circos_links.pl -i Xyz.collinearity -g Xzy.gff
\n";

my ($collinearityfile,$gfffile,$breaksfile,$scorefile,$chrom,$help);
my $bcolor = "vdblue";
my $color = "lblue";
my $ks = 0.5;

GetOptions (
  'i|collinearity=s' => \$collinearityfile,
  'g|gff=s'          => \$gfffile,
  'b|breaks:s'       => \$breaksfile,
  'r|bcolor:s'       => \$bcolor,
  's|score=s'        => \$scorefile,
  'k|ks:f'           => \$ks,
  'c|chrom:s'        => \$chrom,
  'l|color:s'        => \$color,
  'h|help'           => \$help
);

die $usage if $help;
die $usage unless ($collinearityfile && $gfffile && $scorefile);
print STDERR "[INFO] Ks threshold from $scorefile set to $ks\n";

my (%gff_hash, %breaks_hash, %score_hash);

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

if ($breaksfile) {
  open (my $BREAKS, $breaksfile) or die $!;
  while (<$BREAKS>) {
    chomp;
    my @F = split (/\s+/, $_);
    if ($F[-1] =~ "break") {
      $breaks_hash{$F[2]}{CHROM1} = $F[0];
      $breaks_hash{$F[2]}{CHROM2} = $F[4];
    }
  }
}
print STDERR "[INFO] Found ".scalar(keys %breaks_hash)." breaks from $breaksfile\n";

open (my $SCORE, $scorefile) or die $!;
while (<$SCORE>) {
  chomp;
  my @F = split (/\s+/, $_);
  $score_hash{$F[0]}{SCORE} = $F[9];
  $score_hash{$F[0]}{Ka} = $F[10];
  $score_hash{$F[0]}{Ks} = $F[11];
}
close $SCORE;

open (my $COLL, $collinearityfile) or die $!;
open (my $OUT, ">$collinearityfile.circos.links") or die $!;
while (<$COLL>) {
  chomp;
  if ($_ =~ /^#/) {
    next;
  } else {
    my @F = split (/\s+/, $_);
    if ($score_hash{$F[0]}{Ks} <= $ks) {
      if ($breaksfile) {
        if ( $breaks_hash{$F[0]} ) {
          print $OUT join (
            " ",
            $gff_hash{$F[2]}{CHROM},
            $gff_hash{$F[2]}{START},
            $gff_hash{$F[2]}{END},
            $gff_hash{$F[3]}{CHROM},
            $gff_hash{$F[3]}{START},
            $gff_hash{$F[3]}{END},
            "color=$bcolor,thickness=2p",
            "\n"
          );
        } elsif ($chrom) {
          if ($gff_hash{$F[2]}{CHROM} eq $chrom) {
            print $OUT join (
              " ",
              $gff_hash{$F[2]}{CHROM},
              $gff_hash{$F[2]}{START},
              $gff_hash{$F[2]}{END},
              $gff_hash{$F[3]}{CHROM},
              $gff_hash{$F[3]}{START},
              $gff_hash{$F[3]}{END},
              "color=$color,thickness=2p", ##158,202,225
              "\n"
            );
          } else {
            print $OUT join (
              " ",
              $gff_hash{$F[2]}{CHROM},
              $gff_hash{$F[2]}{START},
              $gff_hash{$F[2]}{END},
              $gff_hash{$F[3]}{CHROM},
              $gff_hash{$F[3]}{START},
              $gff_hash{$F[3]}{END},
              "color=vlgrey,thickness=1p", ##158,202,225
              "\n"
            );
          }
        } else {
          print $OUT join (
            " ",
            $gff_hash{$F[2]}{CHROM},
            $gff_hash{$F[2]}{START},
            $gff_hash{$F[2]}{END},
            $gff_hash{$F[3]}{CHROM},
            $gff_hash{$F[3]}{START},
            $gff_hash{$F[3]}{END},
            "color=vlgrey,thickness=1p", ##158,202,225
            "\n"
          );
        }
      } else {
        print $OUT join (
          " ",
          $gff_hash{$F[2]}{CHROM},
          $gff_hash{$F[2]}{START},
          $gff_hash{$F[2]}{END},
          $gff_hash{$F[3]}{CHROM},
          $gff_hash{$F[3]}{START},
          $gff_hash{$F[3]}{END},
          "color=vlgrey,thickness=1p",
          "\n"
        );
      }
    }
  }
}
close $COLL;
close $OUT;

print STDERR "[INFO] Finished on ".`date`."\n";
