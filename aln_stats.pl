#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::AlignIO;
use File::Basename;
use Sort::Naturally;
use Bio::SimpleAlign;

my $usage = "
SYNOPSIS
  Calculates average percent identity (or divergence) for nucleotide alignments.

OPTIONS
  -p|--path   [STRING] : path to files [required]
  -s|--suffix [STRING] : suffix of alignment files to glob from path (default 'fasta')
  -m|--format [STRING] : alignment format (default: 'fasta')
  -o|--out    [FILE]   : output file name (default: 'aln_stats.txt')
  -d|--div             : print divergence rather than identity
  -h|--help            : this message
\n";

my ($path, $divergence, $help);
my $suffix = "fasta";
my $format = "fasta";
my $outfile = "aln_stats.txt";

GetOptions (
  'p|path=s'    => \$path,
  's|suffix:s'  => \$suffix,
  'm|format:s'  => \$format,
  'o|out:s'     => \$outfile,
  'd|div'       => \$divergence,
  'h|help'      => \$help
);

die $usage if $help;
die $usage unless ($path);

my @files = glob ("$path/*.$suffix");
print STDERR "[INFO] Globbed ".scalar(@files)." '*.$suffix' files from $path\n";

open (my $OUT_fh, ">$outfile") or die "$!\n";
print STDERR "[INFO] Printing results to '$outfile'\n";
print STDERR "[INFO] Printing divergence rather than identity\n" if ( $divergence );

foreach my $file (nsort @files) {
  print $OUT_fh "$file\t";
  my $in = Bio::AlignIO->new( -file => $file, -format => $format );
  while ( my $aln_obj = $in->next_aln() ) {
    print $OUT_fh $aln_obj->num_sequences . "\t";
    if ( $divergence ) {
      print $OUT_fh (100 - ($aln_obj->percentage_identity)) . "\n";
    } else {
      print $OUT_fh $aln_obj->percentage_identity . "\n";
    }
  }
}
close $OUT_fh;
print STDERR "[INFO] Finished on ".`date`."\n";

__END__
