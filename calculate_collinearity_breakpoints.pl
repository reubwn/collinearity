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
  -i|--in     [FILE]  : collinearity file
  -g|--gff    [FILE]  : GFF file
  -s|--score  [FILE]  : score file with average Ks per block
  -k|--ks     [FLOAT] : Ks threshold to define homologous block (default <= 0.5)
  -o|--out            : outfile (default=INFILE.breaks)
  -b|--blocks         : also print blocks_per_gene file (default=no)
  -h|--help           : print help
  -m|--morehelp       : print more help

USAGE
  calculate_collinearity_breakpoints.pl -i Xyz.collinearity.kaks -g Xyz.gff -s Xyz.collinearity.kaks.score
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

NOTES
  (1) Cases where homologous blocks are found on the same scaffold are also accounted for in the script.
  (2) Use the script 'calculate_collinearity_palindromes.pl' to search in more detail for these 'linked' cases.
  (3) Cases where a focal region shows >1 homologous region (e.g., recent duplications of A, B and B') are discounted
      as determining orthology from paralogy is not possible with these data (these sites are simply skipped).
\n";

my ($collinearityfile,$gfffile,$scorefile,$blockspergenefile,$help,$morehelp,$verbose);
my $ks = 0.5;

GetOptions (
  'i|in=s'     => \$collinearityfile,
  'g|gff=s'    => \$gfffile,
  's|score=s'  => \$scorefile,
  'b|blocks'   => \$blockspergenefile,
  'k|ks:f'     => \$ks,
  'v|verbose'  => \$verbose,
  'h|help'     => \$help,
  'm|morehelp' => \$morehelp
);

die $usage if $help;
die $usage.$definitions if $morehelp;
die $usage unless ($collinearityfile && $gfffile && $scorefile);

print STDERR "[INFO] Collinearity file: $collinearityfile\n";
print STDERR "[INFO] GFF file: $gfffile\n";
print STDERR "[INFO] Score file: $scorefile\n";

## things we'll need:
my (%collinearity_hash, %homologous_blocks_hash, %score_hash, %score_file_hash, %gff_hash, %blocks_hash, %seen);
my (%total_blocks, %collinear_blocks, %noncollinear_blocks, %noncollinear_blocks_linked_to_same_scaffold, %blocks_linked_to_multiple_homol_regions);

## sort GFF file:
if (system("sort -k1,1 -k3,3 -V $gfffile > $gfffile.sorted") != 0) {
  die "[ERROR] Something up with GNU sort\n";
} else {
  print STDERR "[INFO] Sorted GFF file: $gfffile.sorted\n";
}

## get genewise block participation:
open (my $COLL, $collinearityfile) or die $!;
while (<$COLL>) {
  chomp;
  if ($_ =~ /^#/) {
    next;
  } else {
    $_ =~ s/^\s+|\s+$//g; ## remove leading and trailing whitespaces
    my @F = split (m/\s+/, $_); ## split
    my $block_number;
    if ($F[0]=~m/\d+\-\d+\:/) { ## sometimes columns not formatted properly... :/
      my @a = split (m/\-/, $F[0]);
      $block_number = $a[0];
    } else {
      $F[0] =~ s/\-//;
      $block_number = $F[0]; ## this should also work if the .reformatted file is used
    }
    ## key= gene name; val= @[all blocks that gene is a member of]
    push ( @{$collinearity_hash{$F[-4]}}, $block_number );
    push ( @{$collinearity_hash{$F[-5]}}, $block_number );
  }
}
close $COLL;
print STDERR "[INFO] Parsed ".commify(scalar(keys %collinearity_hash))." genes from $collinearityfile\n";

## get blockwise average Ks and scores:
open (my $SCORE, $scorefile) or die $!;
while (<$SCORE>) {
  chomp;
  my @F = split (/\s+/, $_);
  $score_hash{$F[0]}{SCORE} = $F[9]; ## per block score
  $score_hash{$F[0]}{Ka} = $F[10]; ## per block Ka
  $score_hash{$F[0]}{Ks} = $F[11]; ## per block Ks
}
close $SCORE;
print STDERR "[INFO] Parsed ".commify(scalar(keys %score_hash))." blocks from $scorefile\n";
# print STDERR "[INFO] Ks threshold set to: $ks\n";

my $GENES;
if ($blockspergenefile){
  open ($GENES, ">$gfffile.blocks_per_gene") or die $!;
}

my %pass;
foreach my $gene (nsort keys %collinearity_hash) {
  print $GENES "$gene\t@{$collinearity_hash{$gene}}\n" if ($blockspergenefile);
  foreach my $block (@{$collinearity_hash{$gene}}) {
    push ( @{$homologous_blocks_hash{$gene}}, $block ) if $score_hash{$block}{Ks} <= $ks; ## %homologous_blocks_hash contains ONLY homologous blocks with Ks <= threshold
    $pass{$block}++ if $score_hash{$block}{Ks} <= $ks;
  }
}
close $GENES if ($blockspergenefile);
print STDERR "[INFO] Number of blocks with Ks <= $ks: ".commify(scalar(keys %pass))."\n";

## add block participation to sorted GFF:
open (my $OUT1, ">$gfffile.sorted.painted") or die $!;
open (my $GFF, "$gfffile.sorted") or die $!;
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  if (exists($homologous_blocks_hash{$F[1]})) {
    my @blocks = @{$homologous_blocks_hash{$F[1]}}; ## all homologous blocks for that gene; should be 1 but sometimes more
    print $OUT1 join ("\t", @F, join("|",@blocks), "\n");
  } else {
    print $OUT1 join ("\t", @F, "-", "\n");
  }
}
close $GFF;
close $OUT1;

## reopen PAINTED GFF file:
open (my $PAINTED, "$gfffile.sorted.painted") or die $!;
GENE: while (<$PAINTED>) {
  chomp;
  my @F = split (/\s+/, $_);
  if ($F[4] =~ m/\-/) { ## block does not pass Ks threshold
    next GENE;
  } elsif ($F[4] =~ m/\|/) { ## gene is involved in >1 block
    ## if all genes in given block are involved in >1 block, that block will be omitted
    next GENE;
  } else { ## only analyse those blocks that link exactly two homol regions
    $blocks_hash{$F[4]}{$F[0]}++; ## key= block ID; val= %{chrom names associated with that block}
    push ( @{ $gff_hash{$F[0]} }, $F[4] ) unless $seen{$F[0]}{$F[4]}; ##key= chrom, val= @[order of blocks along chrom]
    $seen{$F[0]}{$F[4]}++;
  }
}
close $PAINTED;
# print STDERR "[INFO] Number of blocks linking exactly 2 homologous regions: ".commify(scalar(keys %blocks_hash))."\n";

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
CHROM: foreach my $focal_chrom (nsort keys %gff_hash) {
  my @blocks1 = @{$gff_hash{$focal_chrom}}; ## array of blocks along chrom

  BLOCK: foreach my $focal_block (@blocks1) {
    my $description = "NULL";
    my $result = "collinear";
    $total_blocks{$focal_block}++;

    ## get identity of $homol_chrom (ie, non-focal chrom linked to focal chrom via block $focal_block)
    my %chroms_linked_to_block = %{ $blocks_hash{$focal_block} }; ## get all chroms linked by $focal_block
    if (scalar(keys(%chroms_linked_to_block))==1) { ## indication that focal and homol chrom are the same!
      # print STDERR "[INFO]+ Block $focal_block is linked to same scaffold (".(join (" ", keys %chroms_linked_to_block)).")\n";
      print $OUT2 join ("\t", $focal_chrom, (join("|",@blocks1)), $focal_block, "NA", "NA", "NA", "Block links two regions on same scaffold", "break", "\n");
      $noncollinear_blocks_linked_to_same_scaffold{$focal_block}++;
      next BLOCK; ## go straight to next block without evaluating code below
    } elsif (scalar(keys(%chroms_linked_to_block))>2) { ## indication that block may link to more than 2 chroms, eg could be caused by recent duplications leading to >2 homologous regions
      print STDERR "[INFO]+ Block $focal_block is linked to multiple scaffolds (".(join (" ", keys %chroms_linked_to_block)).")\n";
      next BLOCK; ## also skip these cases
    }
    my ( $homol_chrom ) = grep { $_ ne $focal_chrom } keys %chroms_linked_to_block;
    # print STDERR "Focal chrom: $focal_chrom\n";
    # print STDERR join ("\t", $focal_block.":", (keys %chroms_linked_to_block), "\n");
    # print STDERR "Homol chrom: $homol_chrom\n\n";

    my @blocks2 = @{$gff_hash{$homol_chrom}}; ## get all blocks on $homol_chrom
    my( $index1 ) = grep { $blocks1[$_] == $focal_block } 0..$#blocks1; ##get index of block in series of blocks on focal chrom
    my( $index2 ) = grep { $blocks2[$_] == $focal_block } 0..$#blocks2; ##get index of HOMOLOGOUS block on HOMOLOGOUS chrom2

    if ( ("@blocks1" =~ "@blocks2") || ("@blocks1" =~ join (" ", reverse(@blocks2))) ) {
      $description = "B subset of A";
      $collinear_blocks{$focal_block}++;
    } elsif ( ("@blocks2" =~ "@blocks1") || ("@blocks2" =~ join (" ", reverse(@blocks1))) ) {
      $description = "A subset of B";
      $collinear_blocks{$focal_block}++;
    } else {
      if ( ($index1 == 0 || $index1 == $#blocks1) && ($index2 == 0 || $index2 == $#blocks2) ) {
            $description = "both terminal";
            $collinear_blocks{$focal_block}++;
          } elsif ( ($index1 != 0 || $index1 != $#blocks1) && ($index2 == 0 || $index2 == $#blocks2) ) {
            $description = "B#$focal_block is terminal";
            my ($i,$j)=($index1,$index2);

            ## traverse blocks forwards:
            if ($index2==0) {##subject block is 3' terminal
              unless ($blocks1[($i+1)] == $blocks2[($j+1)]) {
                $description .= " (3')";
                $result = "break";
                $noncollinear_blocks{$focal_block}++;
              } else {$collinear_blocks{$focal_block}++;}
            ## traverse blocks backwards:
          } elsif ($index2==$#blocks2) { ##subject block is 5' terminal
              unless ($blocks1[($i-1)] == $blocks2[($j-1)]) {
                $description .= " (5')";
                $result = "break";
                $noncollinear_blocks{$focal_block}++;
              } else {$collinear_blocks{$focal_block}++;}
            }

          } elsif ( ($index1 == 0 || $index1 == $#blocks1) && ($index2 != 0 || $index2 != $#blocks2) ) {
            $description = "A#$focal_block is terminal";
            my ($i,$j)=($index1,$index2);

            ## traverse blocks forwards:
            if ($index1==0) {##A block is 3' terminal
              unless ($blocks1[($i+1)] == $blocks2[($j+1)]) {
                $description .= " (3')";
                $result = "break";
                $noncollinear_blocks{$focal_block}++;
              } else {$collinear_blocks{$focal_block}++;}
            ## traverse blocks backwards:
          } elsif ($index1==$#blocks1) { ##focal block is 5' terminal
              unless ($blocks1[($i-1)] == $blocks2[($j-1)]) {
                $description .= " (5')";
                $result = "break";
                $noncollinear_blocks{$focal_block}++;
              } else {$collinear_blocks{$focal_block}++;}
            }

          } else {
            $description = "neither terminal";
            $result = "break";
            $noncollinear_blocks{$focal_block}++;
          }
    }

    print $OUT2 "$focal_chrom\t";
    print $OUT2 join ("|", @blocks1);
    print $OUT2 "\t$focal_block\t";
    print $OUT2 join ("|", keys %{ $blocks_hash{$focal_block} });
    print $OUT2 "\t$homol_chrom\t";
    print $OUT2 join ("|", @blocks2);
    print $OUT2 "\t$description\t$result\n";
  }
}
close $OUT2;

## get some numbers
my ($c,$n,$nl,$T) = (0,0,0,0);
$c  = scalar(keys %collinear_blocks) if (%collinear_blocks);
$n  = scalar(keys %noncollinear_blocks) if (%noncollinear_blocks);
$nl = scalar(keys %noncollinear_blocks_linked_to_same_scaffold) if (%noncollinear_blocks_linked_to_same_scaffold);
$T  = scalar(keys %total_blocks) if (%total_blocks);

print STDERR "[INFO] ~~~\n";
print STDERR "[INFO] Number of blocks analysed: ".commify($T)."\n";
print STDERR "[INFO] Number of collinear blocks: ".commify($c)." (".percentage($c,$T).")\n";
print STDERR "[INFO] Number of noncollinear blocks (different scaffolds): ".commify($n)." (".percentage($n,$T).")\n";
print STDERR "[INFO] Number of noncollinear blocks (linked on same scaffold): ".commify($nl)." (".percentage($nl,$T).")\n";
print STDERR "[INFO] Total number of collinearity breaks: ".commify(($n+$nl))." (".percentage(($n+$nl),$T).")\n"; ## /2 because each break is counted twice, from the perspective of both involved chroms
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
