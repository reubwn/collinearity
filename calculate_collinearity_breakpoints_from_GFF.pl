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
  Given a sorted GFF annotated with block participation, will search for breaks
  in collinearity defined as occurrences where homologous blocks cannot be
  aligned along scaffolds without some rearrangement

OPTIONS
  -g|--gff [FILE] : painted GFF file annotated with block membership
  -o|--out        : outfile (default=INFILE.breaks)
  -h|--help       : print this message

USAGE
  calculate_collinearity_breakpoints.pl -g Xyz.gff.painted
\n";

my ($scorefile,$gfffile,$help);
my $ks = 0.5;

GetOptions (
  'gff|g=s'          => \$gfffile,
  'help|h'           => \$help
);

die $usage if $help;
die $usage unless ($gfffile);

#print STDERR "[INFO] Collinearity file: $scorefile\n";
print STDERR "[INFO] Painted GFF file: $gfffile\n";
#print STDERR "[INFO] Ks threshold: $ks\n";

my (%gff_hash, %blocks_hash, %seen);
my ($collinear_blocks, $noncollinear_blocks, $total_blocks) = (0,0,0);

open (my $GFF, $gfffile) or die $!;
while (<$GFF>) {
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
close $GFF;
#print STDERR "[INFO] Total chroms with blocks: ".scalar(keys %gff_hash)."\n";

open (my $OUT, ">$gfffile.breaks") or die $!;
print $OUT join ("\t",
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

    print $OUT "$chrom\t";
    print $OUT join ("|", @blocks1);
    print $OUT "\t$focal_block\t";
    print $OUT join ("|", @{$blocks_hash{$focal_block}});
    print $OUT "\t@arr\t";
    print $OUT join ("|", @blocks2) if (defined($arr[0]));
    print $OUT "\t$description\t$result\n";
    $total_blocks++;
  }
}
close $OUT;

print STDERR "[INFO] Total blocks: ".commify($total_blocks)."\n";
print STDERR "[INFO] Number of collinear blocks: ".commify($collinear_blocks)." (".percentage($collinear_blocks,$total_blocks).")\n";
print STDERR "[INFO] Number of noncollinear blocks: ".commify($noncollinear_blocks)." (".percentage($noncollinear_blocks,$total_blocks).")\n";
print STDERR "[INFO] Number of collinearity breaks: ".($noncollinear_blocks/2)."\n";
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
