#!/usr/bin/env perl

## author: reubwn August 2017

use strict;
use warnings;
use Getopt::Long;
use Array::Utils qw(:all);
use List::MoreUtils qw /pairwise/;
use Sort::Naturally;

my $usage = "
SYNOPSIS
  Searches for breaks in collinearity defined as occurrences where homologous
  blocks cannot be aligned along scaffolds without some rearrangement

  Uses GNU 'sort' to presort GFF file.

OPTIONS
  -i|--in     [FILE]  : collinearity file (use .reformatted!)
  -g|--gff    [FILE]  : GFF file
  -s|--score  [FILE]  : score file with average Ks per block
  -k|--ks     [FLOAT] : Ks threshold to define homologous block (default <= 0.5)
  -o|--out            : outfile (default=INFILE.breaks)
  -b|--blocks         : also print blocks_per_gene file (default=no)
  -h|--help           : print help
  -m|--morehelp       : print more help

USAGE
  calculate_collinearity_breakpoints.pl -i Xy.collinearity.kaks.reformatted -g Xy.gff -s Xy.collinearity.kaks.score
\n";

my $definitions = "~~~
BREAKPOINT DEFINITION
  A collinearity break is introduced if there is a mistmatch in the identity of the upstream
  or downstream homologous blocks from one collinear region to another. Thus, if a[i] is the identity of the focal block on chrom a,
  homologous with block b[j] on chrom b, collinearity is broken if a[i+/-1] != b[j+/-1]. A caveat is when the focal block is located
  at the terminus of a scaffold - here, collinearity is not broken if the homologous scaffold can be orientated 'away' from the focal scaffold

EXAMPLES
  SUBJECT is SINGLETON
    a----1--2-----3-4----5---
    b     --2---
  a and b are inferred collinear, since b contains only 1 block

  FOCAL and SUBJECT are TERMINAL
    a----1--2-----3-4----5---
    c                 ---5---6--7--8--
  a and c are collinear, since orientation of scaffolds WRT each other is unknown and may align as shown

  FOCAL or SUBJECT are TERMINAL
    a----1--2-----3-4----5---
    d              -4--10-11----12---
  a and d show a break in collinearity, since i=4, a[4+1]=5, but b[j+1]=10, and similarly a[i-1]=2 but still cannot align with b if orientated in the other direction

  neither FOCAL or SUBJECT are TERMINAL
    a----1--2-----3-4----5---
    e  18---2--19-20----21-----
  a and e are clearly not collinear

RULES
  (1) collinearity cannot be broken by scaffolds containing only 1 collinear block.
  (2) collinearity cannot be broken when focal and subject blocks are BOTH terminal on their respective scaffolds
\n";

my ($collinearityfile,$gfffile,$scorefile,$blockspergenefile,$help,$morehelp,$verbose);
my $ks = 0.5;

GetOptions (
  'i|collinearity=s' => \$collinearityfile,
  'g|gff=s'          => \$gfffile,
  's|score=s'        => \$scorefile,
  'b|blocks'         => \$blockspergenefile,
  'k|ks:f'           => \$ks,
  'v|verbose'        => \$verbose,
  'h|help'           => \$help,
  'm|morehelp'       => \$morehelp
);

die $usage if $help;
die $usage.$definitions if $morehelp;
die $usage unless ($collinearityfile && $gfffile && $scorefile);

print STDERR "[INFO] Collinearity file: $collinearityfile\n";
print STDERR "[INFO] GFF file: $gfffile\n";
print STDERR "[INFO] Score file: $scorefile\n";
print STDERR "[INFO] Ks threshold: $ks\n";

## die unless refomatted collinearity file:
unless ($collinearityfile =~ m/refomatted$/) {
  die "[ERROR] Collinearity file: $collinearityfile is not reformatted\n";
}

## things we'll need:
my (%collinearity_hash, %homologous_blocks_hash, %score_hash, %score_file_hash, %gff_hash, %blocks_hash, %seen);
my ($collinear_blocks, $noncollinear_blocks, $total_blocks) = (0,0,0);

## sort GFF file:
if (system("sort -k1,1 -k3,3 -V $gfffile > $gfffile.sorted") != 0) {
  die "[ERROR] Something up with GNU sort\n";
} else {
  print STDERR "[INFO] Sorted GFF file: $gfffile.sorted\n";
}

## get blockwise average Ks and scores:
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

## get genewise block participation:
open (my $COLL, $collinearityfile) or die $!;
while (<$COLL>) {
  chomp;
  if ($_ =~ /^#/) {
    next;
  } else {
    my @F = split (/\s+/, $_);
    ## key= gene name; val= @{all blocks that gene is a member of}
    push ( @{$collinearity_hash{$F[2]}}, $F[0] );
    push ( @{$collinearity_hash{$F[3]}}, $F[0] );
  }
}
close $COLL;
print STDERR "[INFO] Parsed ".scalar(keys %collinearity_hash)." genes from $collinearityfile\n";

my $GENES;
if ($blockspergenefile){
  open ($GENES, ">$gfffile.blocks_per_gene") or die $!;
}

foreach my $gene (nsort keys %collinearity_hash) {
  print $GENES "$gene\t@{$collinearity_hash{$gene}}\n" if ($blockspergenefile);
  foreach my $block (@{$collinearity_hash{$gene}}) {
    push ( @{$homologous_blocks_hash{$gene}}, $block ) if $score_hash{$block}{Ks} <= $ks; ## %homologous_blocks_hash contains ONLY homologous blocks with Ks <= threshold
  }
}
close $GENES if ($blockspergenefile);

## add block participation to sorted GFF:
open (my $OUT1, ">$gfffile.sorted.painted") or die $!;
open (my $GFF, "$gfffile.sorted") or die $!;
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  if (exists($homologous_blocks_hash{$F[1]})) {
    my @blocks = @{$homologous_blocks_hash{$F[1]}}; ## all homologous blocks for that gene; should be 1 but sometimes more
    if (scalar(@blocks)>1) {
      print STDERR "[WARN] Gene $F[1] has >1 homologous block: @blocks\n" if $verbose;
    } else {
      print $OUT1 join ("\t", @F, @blocks, "\n");
      $seen{$blocks[0]}++;
    }
  } else {
    print $OUT1 join ("\t", @F, "-", "\n");
  }
}
close $GFF;
close $OUT1;

## reopen PAINTED GFF file:
open (my $PAINTED, "$gfffile.sorted.painted") or die $!;
while (<$PAINTED>) {
  chomp;
  my @F = split (/\s+/, $_);
  if ($F[4] =~ /\-/) {
    next;
  } else {
    push ( @{ $blocks_hash{$F[4]} }, $F[0]) unless $seen{$F[0]}{$F[4]}; ##key= block, val= @[chroms involved in block]
    push ( @{ $gff_hash{$F[0]} }, $F[4] ) unless $seen{$F[0]}{$F[4]}; ##key= chrom, val= @[order of blocks along chrom]
    $seen{$F[0]}{$F[4]}++;
  }
}
close $PAINTED;

open (my $OUT2, ">$gfffile.sorted.painted.breaks") or die $!;
print $OUT2 join ("\t",
                 "A.CHROM",
                 "A.BLOCKS",
                 "FOCAL",
                 "HOMCHROMS",
                 "B.CHROM",
                 "B.BLOCKS",
                 "DESCRIPTION",
                 "RESULT",
                 "\n"
                );

## main code for detecting breakpoints:
foreach my $chrom (nsort keys %gff_hash) {
  my @blocks1 = @{$gff_hash{$chrom}};

  foreach my $focal_block (@blocks1) {
    my $description = "NULL";
    my $result = "collinear";

    my $i = 0;
    my @arr = @{$blocks_hash{$focal_block}};
    $i++ until $arr[$i] eq $chrom;
    splice(@arr, $i, 1);

    my @blocks2 = @{$gff_hash{$arr[0]}} if (defined($arr[0]));
    my( $index1 ) = grep { $blocks1[$_] == $focal_block } 0..$#blocks1; ##get index of block in series of blocks on same chrom
    my( $index2 ) = grep { $blocks2[$_] == $focal_block } 0..$#blocks2; ##get index of HOMOLOGOUS block on HOMOLOGOUS chrom
    print STDERR "$chrom, $index1\n" unless ($index1);

    if ( ("@blocks1" =~ "@blocks2") || ("@blocks1" =~ join (" ", reverse(@blocks2))) ) {
      $description = "B subset of A";
      $collinear_blocks++;
    } elsif ( ("@blocks2" =~ "@blocks1") || ("@blocks2" =~ join (" ", reverse(@blocks1))) ) {
      $description = "A subset of B";
      $collinear_blocks++;
    } else {
      if ( ($index1 == 0 || $index1 == $#blocks1) && ($index2 == 0 || $index2 == $#blocks2) ) {
            $description = "both terminal";
            $collinear_blocks++;
          } elsif ( ($index1 != 0 || $index1 != $#blocks1) && ($index2 == 0 || $index2 == $#blocks2) ) {
            $description = "B#$focal_block is terminal";
            my ($i,$j)=($index1,$index2);

            ## traverse blocks forwards:
            if ($index2==0) {##subject block is 3' terminal
              unless ($blocks1[($i+1)] == $blocks2[($j+1)]) {
                $description .= " (3')";
                $result = "break";
                $noncollinear_blocks++;
              } else {$collinear_blocks++;}
            ## traverse blocks backwards:
          } elsif ($index2==$#blocks2) { ##subject block is 5' terminal
              unless ($blocks1[($i-1)] == $blocks2[($j-1)]) {
                $description .= " (5')";
                $result = "break";
                $noncollinear_blocks++;
              } else {$collinear_blocks++;}
            }

          } elsif ( ($index1 == 0 || $index1 == $#blocks1) && ($index2 != 0 || $index2 != $#blocks2) ) {
            $description = "A#$focal_block is terminal";
            my ($i,$j)=($index1,$index2);

            ## traverse blocks forwards:
            if ($index1==0) {##A block is 3' terminal
              unless ($blocks1[($i+1)] == $blocks2[($j+1)]) {
                $description .= " (3')";
                $result = "break";
                $noncollinear_blocks++;
              } else {$collinear_blocks++;}
            ## traverse blocks backwards:
          } elsif ($index1==$#blocks1) { ##focal block is 5' terminal
              unless ($blocks1[($i-1)] == $blocks2[($j-1)]) {
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

    print $OUT2 "$chrom\t";
    print $OUT2 join ("|", @blocks1);
    print $OUT2 "\t$focal_block\t";
    print $OUT2 join ("|", @{$blocks_hash{$focal_block}});
    print $OUT2 "\t@arr\t";
    print $OUT2 join ("|", @blocks2) if (defined($arr[0]));
    print $OUT2 "\t$description\t$result\n";
    $total_blocks++;
  }
}
close $OUT2;

print STDERR "[INFO] Number of collinear blocks: ".commify($collinear_blocks)." (".percentage($collinear_blocks,$total_blocks).")\n";
print STDERR "[INFO] Number of noncollinear blocks: ".commify($noncollinear_blocks)." (".percentage($noncollinear_blocks,$total_blocks).")\n";
print STDERR "[INFO] Number of collinearity breaks: ".($noncollinear_blocks/2)."\n";
print STDERR "[INFO] Results written to: $gfffile.sorted.painted.breaks\n";
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
# ## (1) collinearity cannot be broken by scaffolds containing only 1 collinear block.
# ## (2) collinearity cannot be broken when focal and subject blocks are BOTH terminal on their respective scaffolds
