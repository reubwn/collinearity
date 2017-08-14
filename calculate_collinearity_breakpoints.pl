#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw /pairwise/;
use Sort::Naturally;

my $usage = "
SYNOPSIS

OUTPUT

OPTIONS
  -i|--in [FILE] : collinearity score file (annotated with Ks)
  -k|--ks        : only examine blocks with Ks <= this threshold
  -o|--out       : outfile (default=INFILE.breaks)
  -h|--help      : print this message

USAGE
  calculate_collinearity_breakpoints.pl -i Xyz.collinearity.kaks.score -k 0.5

\n";

my ($scorefile,$help,$debug);
my $ks = 0.5;

GetOptions (
  'collinearity|i=s' => \$scorefile,
  'ks|k:f'           => \$ks,
  'help|h'           => \$help,
  'debug|d'          => \$debug
);

die $usage if $help;
die $usage unless ($scorefile);
die "[ERROR] Hmm, suspicious lack of \*.score on $scorefile... is this the right file?\n\n";

print STDERR "[INFO] Collinearity file: $scorefile\n";
print STDERR "[INFO] Ks threshold: $ks\n";

my (%blocks1,%blocks2,%chroms1,%chroms2);

open (my $IN, $scorefile) or die $!;
while (<$IN>) {
  next if $. == 1; ##skip header on line1
  my @F = split (/\s+/, $_);
  if ($F[11] <= $ks) { ##only select homologs with Ks < some threshold
    $chroms1{$F[0]} = $F[1]; ##key= blockname; val=chrom where block exists
    push ( @{$blocks1{$F[1]}}, $F[0] ); ##key=chrom1 name; val= @{ blocks on chrom1 name }
    $chroms2{$F[0]} = $F[2]; ##key= blockname; val=chrom where block exists
    push ( @{$blocks2{$F[2]}}, $F[0] ); ##key=chrom2 name; val= @{ blocks on chrom2 name }
  }
}
close $IN;

open (my $OUT, ">$scorefile.breaks") or die $!;
print $OUT join ("\t",
                 "FOCAL",
                 "A.CHROM",
                 "A.BLOCKS",
                 "B.CHROM",
                 "B.BLOCKS",
                 "DESCRIPTION",
                 "RESULT",
                 "\n"
                );

my ($collinear_blocks, $noncollinear_blocks, $total_blocks) = (0,0,0);
foreach my $block (nsort keys %chroms1) {

  my @blocksOnSameChrom = @{ $blocks1{$chroms1{$block}} }; ##all blocks residing on the same chrom
  my @blocksOnHomologousChrom = @{ $blocks2{$chroms2{$block}} }; ##all other blocks residing on (what should be) the HOMOLOGOUS chrom
  my( $index1 ) = grep { $blocksOnSameChrom[$_] == $block } 0..$#blocksOnSameChrom; ##get index of block in series of blocks on same chrom
  my( $index2 ) = grep { $blocksOnHomologousChrom[$_] == $block } 0..$#blocksOnHomologousChrom; ##get index of HOMOLOGOUS block on HOMOLOGOUS chrom
  my $description = "NULL";
  my $result = "collinear";

  if ( ("@blocksOnSameChrom" =~ "@blocksOnHomologousChrom") || ("@blocksOnSameChrom" =~ join(" ",reverse(@blocksOnHomologousChrom))) ) {
    $description = "B subset of A";
    $collinear_blocks++;
  } elsif ( ("@blocksOnHomologousChrom" =~ "@blocksOnSameChrom") || (join(" ",reverse(@blocksOnHomologousChrom)) =~ "@blocksOnSameChrom") ) {
    $description = "A subset of B";
    $collinear_blocks++;
  } else {
    if ( ($index1 == 0 || $index1 == $#blocksOnSameChrom) && ($index2 == 0 || $index2 == $#blocksOnHomologousChrom) ) {
      $description = "both terminal";
      $collinear_blocks++;
    } elsif ( ($index1 != 0 || $index1 != $#blocksOnSameChrom) && ($index2 == 0 || $index2 == $#blocksOnHomologousChrom) ) {
      $description = "B#$block is terminal";
      my ($i,$j)=($index1,$index2);

      ## traverse blocks forwards:
      if ($index2==0) {##subject block is 3' terminal
        unless ($blocksOnSameChrom[($i+1)] == $blocksOnHomologousChrom[($j+1)]) {
          $description .= " (3')";
          $result = "break";
          $noncollinear_blocks++;
        } else {$collinear_blocks++;}
      ## traverse blocks backwards:
    } elsif ($index2==$#blocksOnHomologousChrom) { ##subject block is 5' terminal
        unless ($blocksOnSameChrom[($i-1)] == $blocksOnHomologousChrom[($j-1)]) {
          $description .= " (5')";
          $result = "break";
          $noncollinear_blocks++;
        } else {$collinear_blocks++;}
      }

    } elsif ( ($index1 == 0 || $index1 == $#blocksOnSameChrom) && ($index2 != 0 || $index2 != $#blocksOnHomologousChrom) ) {
      $description = "A#$block is terminal";
      my ($i,$j)=($index1,$index2);

      ## traverse blocks forwards:
      if ($index1==0) {##A block is 3' terminal
        unless ($blocksOnSameChrom[($i+1)] == $blocksOnHomologousChrom[($j+1)]) {
          $description .= " (3')";
          $result = "break";
          $noncollinear_blocks++;
        } else {$collinear_blocks++;}
      ## traverse blocks backwards:
    } elsif ($index1==$#blocksOnSameChrom) { ##focal block is 5' terminal
        unless ($blocksOnSameChrom[($i-1)] == $blocksOnHomologousChrom[($j-1)]) {
          $description .= " (5')";
          $result = "break";
          $noncollinear_blocks++;
        } else {$collinear_blocks++;}
      }

    } else {
      $description = "neither terminal";
      $result = "break";
      $noncollinear_blocks++;
    }
  }

  ## print to OUT
  print $OUT join ("\t",
              $block,
              $chroms1{$block},
              join("|",@blocksOnSameChrom),
              $chroms2{$block},
              join("|",@blocksOnHomologousChrom),
              "$description",
              "$result",
              "\n"
              );

  $total_blocks++;
}
close $OUT;

print STDERR "[INFO] Total blocks: ".commify($total_blocks)."\n";
print STDERR "[INFO] Number of collinear blocks: ".commify($collinear_blocks)." (".percentage($collinear_blocks,$total_blocks).")\n";
print STDERR "[INFO] Number of collinearity breaks: ".commify($noncollinear_blocks)." (".percentage($noncollinear_blocks,$total_blocks).")\n";
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
    my $places = "\%.2f"; ## default is two decimal places
    if (exists $_[2]){$places = "\%.".$_[2]."f";};
    my $float = (($numerator / $denominator)*100);
    my $rounded = sprintf("$places",$float);
    return "$rounded\%";
}
 __END__

# ## BREAKPOINT DEFINITION:
# ## a collinearity break is introduced if there is a mistmatch in the identity of the upstream
# ## or downstream homologous blocks from one collinear region to another. Thus, if a[i] is the identity of the focal block on chrom a,
# ## homologous with block b[j] on chrom b, collinearity is broken if a[i+/-1] != b[j+/-1]. A caveat is when the focal block is located
# ## at the terminus of a scaffold - here, collinearity is not broken if the homologous scaffold can be orientated 'away' from the focal scaffold
#
# ## some examples:
# ## SUBJECT is SINGLETON
# ## a----1--2-----3-4----5---
# ## b     --2---
# ## a and b are inferred collinear, since b contains only 1 block
#
# ## FOCAL and SUBJECT are TERMINAL
# ## a----1--2-----3-4----5---
# ## c                 ---5---6--7--8--
# ## a and c are collinear, since orientation of scaffolds WRT each other is unknown and may align as shown
#
# ## FOCAL or SUBJECT are TERMINAL
# ## a----1--2-----3-4----5---
# ## d              -4--10-11----12---
# ## a and d show a break in collinearity, since i=4, a[4+1]=5, but b[j+1]=10, and similarly a[i-1]=2 but still cannot align with b if orientated in the other direction
#
# ## neither FOCAL or SUBJECT are TERMINAL
# ## a----1--2-----3-4----5---
# ## e  18---2--19-20----21-----
# ## a and e are clearly not collinear
#
# ## (1) collinearity cannot be broken by scaffolds containing only 1 collinear block. These scaffolds/blocks are essentially ignored.
# ## (2) collinearity cannot be broken when focal and subject blocks are BOTH terminal on their respective scaffolds
# ## (3)
