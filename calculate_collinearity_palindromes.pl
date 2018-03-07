#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw /pairwise/;
use Sort::Naturally;

my $usage = "
SYNOPSIS
  Searched for cases of homologous regions physically linked on the same scaffolds.
  Decomposes any such cases into tandem arrays or palindromes.

OUTPUT
  1. arrays file       : all results
  2. arrays.tandem     : tandem arrays
  3. arrays.palindrome : palindromic arrays

OPTIONS
  -i|--in   [FILE] : collinearity file from MCScanX (annotated with Ks)
  -g|--gff  [FILE] : GFF file
  -k|--ks          : only examine blocks with Ks <= this threshold
  -h|--help        : print this message
\n";

my ($collinearityfile,$gfffile,$help,$debug);
my $ks = 0.5;

GetOptions (
  'i|in=s'  => \$collinearityfile,
  'g|gff=s' => \$gfffile,
  'k|ks:f'  => \$ks,
  'h|help'  => \$help,
  'd|debug' => \$debug
);

die $usage if $help;
die $usage unless ($collinearityfile && $gfffile);

print STDERR "[INFO] Collinearity file: $collinearityfile\n";
print STDERR "[INFO] GFF file: $gfffile\n";
print STDERR "[INFO] Ks threshold: $ks\n";

my (%gff,%blocks,%palindrome,%tandem,%mixed);
my ($chrom1,$chrom2,$orientation);

## parse GFF:
open (my $GFF, $gfffile) or die $!;
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/);
  $gff{$F[1]} = {
    'chrom' => $F[0],
    'start' => $F[2],
    'end'   => $F[3]
  };
}
close $GFF;

## parse collinearity:
open (my $IN, $collinearityfile) or die $!;
while (<$IN>) {
  chomp;
  if ($_ =~ m/^#/) {
    if ($_ =~ m/\s(\w+\d+)\&(\w+\w\d+)\s(plus|minus)$/) { ## get chrom names and strand orientation of block 2
      $chrom1 = $1;
      $chrom2 = $2;
      $orientation = $3;
      next;
    }
  } else {
    $_ =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces
    my @F = split (m/\s+/, $_); ##split
    (my $block = $F[0]) =~ s/\-$//;
    if ( ($chrom1 eq $chrom2) && ($F[-1] <= $ks) ) {

      #print "$block $chrom1 $chrom2 $F[-5] $F[-4] $gff{$F[-5]}{start} $gff{$F[-4]}{start} $distance\n";
      $blocks{$block}{chrom} = $chrom1;
      push ( @{$blocks{$block}{genes1}}, $F[-5] );
      push ( @{$blocks{$block}{genes2}}, $F[-4] );
      push ( @{$blocks{$block}{starts1}}, $gff{$F[-5]}{start} );
      push ( @{$blocks{$block}{starts2}}, $gff{$F[-4]}{start} );

    }
  }
}

open (my $OUT, ">$collinearityfile.arrays") or die $!;
print $OUT join ("\t",
  "BLOCK",
  "CHROM",
  "GENE1",
  "START1",
  "GENE2",
  "START2",
  "TYPE",
  "\n"
);

foreach my $block (nsort keys %blocks) {
  my @types;
  for my $i (0..($#{$blocks{$block}{genes1}} -1)) {
    my $type;
    if ( ${$blocks{$block}{starts2}}[($i+1)] < ${$blocks{$block}{starts2}}[$i] ) {
      push (@types, "P");
      $type = "P";
    } else {
      push (@types, "T");
      $type = "T";
    }
    print $OUT join ("\t",
      $block,
      $blocks{$block}{chrom},
      ${$blocks{$block}{genes1}}[$i],
      ${$blocks{$block}{starts1}}[$i],
      ${$blocks{$block}{genes2}}[$i],
      ${$blocks{$block}{starts2}}[$i],
      $type,
      "\n"
    );
  }
  if ( ("@types" =~ m/P/) && ("@types" !~ m/T/) ) {
    $palindrome{$block} = $blocks{$block}{chrom};
  } elsif ( ("@types" =~ m/T/) && ("@types" !~ m/P/) ) {
    $tandem{$block} = $blocks{$block}{chrom};
  } else {
    $mixed{$block} = $blocks{$block}{chrom};
  }
}
close $OUT;

if (scalar(keys %palindrome) > 0) {
  open (my $PAL, ">$collinearityfile.arrays.palindrome") or die $!;
  print $PAL join ("\t", "BLOCK","CHROM","TYPE","\n");
  foreach (nsort keys %palindrome) {
    print $PAL join ("\t", $_,$palindrome{$_},"palindrome","\n");
  }
  close $PAL;
}
if (scalar(keys %tandem) > 0) {
  open (my $TAN, ">$collinearityfile.arrays.tandem") or die $!;
  print $TAN join ("\t", "BLOCK","CHROM","TYPE","\n");
  foreach (nsort keys %tandem) {
    print $TAN join ("\t", $_,$tandem{$_},"tandem","\n");
  }
  close $TAN;
}

print STDERR "[INFO] Number of blocks on same chrom: ".scalar(keys %blocks)."\n";
print STDERR "[INFO] Number of blocks in tandem array: ".scalar(keys %tandem)." (".percentage(scalar(keys %tandem),scalar(keys %blocks)).")\n";
print STDERR "[INFO] Number of blocks in palindromic array: ".scalar(keys %palindrome)." (".percentage(scalar(keys %palindrome),scalar(keys %blocks)).")\n";
print STDERR "[INFO] Finished on ".`date`."\n";

############################# subs

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub percentage {
    my $numerator = $_[0];
    my $denominator = $_[1];
    my $result;
    if ($denominator > 0) {
      my $places = "\%.2f"; ## default is two decimal places
      if (exists $_[2]){$places = "\%.".$_[2]."f";};
      my $float = (($numerator / $denominator)*100);
      my $result = sprintf("$places",$float);
    } else {
      $result = 0;
    }
    return "$result\%";
}
 __END__
