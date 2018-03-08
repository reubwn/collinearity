#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Sort::Naturally;
use File::Path 'rmtree';
use Bio::Align::Utilities qw(aa_to_dna_aln);

my $usage = "
SYNOPSIS
  Generates alignments for pairs of genes with a given minimum Ks from an MCScanX collinearity file (annotated with Ka/Ks).
  If protein and CDS files (-p & -c) are provided, will do protein alignment followed by backtranslation; if gff and genome files (-g & -f) are provided, will do alignment of whole gene region.
  Aligner (clustalo for proteins, mafft for nucleotides) must be executable and in \$PATH.
  Assumes Ks values are in the final column of the collinearity file (default).

OPTIONS
  -i|--in      [FILE]  : MCScanX collinearity file (annotated with Ka/Ks)
  -p|--prot    [FILE]  : fasta file of protein sequences
  -c|--cds     [FILE]  : fasta file of corresponding CDS (nucleotide)
  -f|--fasta   [FILE]  : fasta file of genome
  -g|--gff     [FILE]  : gff file of gene regions
  -o|--out     [STR]   : outfile prefix (default = 'ALN'); alignments will be written to ALN.<NUM>.fasta
  -d|--outdir  [DIR]   : dirname to save alignments (default = 'alignments')
  -k|--minks   [FLOAT] : minimum Ks between genes (default >= 0.5)
  -t|--threads [INT]   : number of aligner threads (default = 1)
  -h|--help            : this message

USAGE
  (1) To generate codon alignments for CDS using Clustal-Omega:
  >> get_alignments_from_MCScanX.pl -i <Xy.collinearity> -p </path/to/proteins/> -c </path/to/cds/> -d <dirname> -k <minKs> -t <threads>
  (2) To generate whole gene (ie exonts + introns) alignments using MAFFT:
  >> get_alignments_from_MCScanX.pl -i <Xy.collinearity> -f </path/to/fasta/> -c </path/to/gff/> -d <dirname> -k <minKs> -t <threads>

\n";

my ($infile, $proteinfile, $cdsfile, $fastafile, $gfffile, $help);
my $minKs = 0.5;
my $out = "ALN";
my $outdir = "alignments";
my $threads = 1;

GetOptions (
  'in|i=s'      => \$infile,
  'prot|p:s'    => \$proteinfile,
  'cds|c:s'     => \$cdsfile,
  'fasta|f:s'   => \$fastafile,
  'gff|g:s'     => \$gfffile,
  'out|o:s'     => \$out,
  'outdir|d:s'  => \$outdir,
  'minks|k:f'   => \$minKs,
  'threads|t:i' => \$threads,
  'help|h'      => \$help
);

die $usage if $help;
die $usage unless ($infile);
die $usage unless (($proteinfile && $cdsfile) || ($fastafile && $gfffile));

print STDERR "[INFO] Collinearity file: $infile\n";
print STDERR "[INFO] Minimum Ks: $minKs\n";
print STDERR "[INFO] Number of threads set to: $threads\n";
print STDERR "[INFO] Parsing collinearity file...\n";

## parse collinearity file:
open (my $IN, $infile) or die "Cannot open $infile: $!\n\n";
my %pairs;
while (<$IN>) {
  chomp;
  next if /^\#/;
  $_ =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces
  my @F = split (/\s+/,$_);
  ## work from end of array as columns in collinearity file not consistently formatted
  $pairs{$F[-5]}=$F[-4] if $F[-1] >= $minKs;
}
close $IN;
print STDERR "[INFO] Number of pairs with Ks >= $minKs: ".commify(scalar(keys(%pairs)))."\n";

##################################################
## codon based alignment of CDS regions only
##################################################

if ($proteinfile && $cdsfile) {
  print STDERR "[INFO] Will run ClustalO alignment on CDS regions & backtranslate to nucleotides\n";
  print STDERR "[INFO] Parsing protein and CDS files...\n";
  ## parse proteins and CDS:
  my (%protein_hash, %cds_hash);
  my $in_p = Bio::SeqIO->new( -file => $proteinfile, -format => 'fasta' );
  while (my $seq = $in_p->next_seq() ) {
    $protein_hash{$seq->display_id()} = $seq->seq();
  }
  print "[INFO] Fetched ".commify(scalar(keys %protein_hash))." proteins from $proteinfile\n";
  my $in_c = Bio::SeqIO->new( -file => $cdsfile, -format => 'fasta' );
  while (my $seq = $in_c->next_seq() ) {
    $cds_hash{$seq->display_id()} = $seq->seq();
  }
  print "[INFO] Fetched ".commify(scalar(keys %cds_hash))." CDS from $cdsfile\n";
  die "[ERROR] No sequences found in $proteinfile or $cdsfile!\n" if ((scalar(keys %protein_hash) == 0) || (scalar(keys %cds_hash) == 0));

  ## generate alignments for all pairs in %pairs:
  my $n = 1;
  open (my $OUT, ">$infile.map") or die "[ERROR] Cannot open $infile.map: $!\n\n";
  print $OUT "NUM\tSEQ1\tSEQ2\n";
  foreach (nsort keys %pairs) {
    print STDERR "\r[INFO] Working on pair \#$n: $_, $pairs{$_}";$| = 1;
    print $OUT "$n\t$_\t$pairs{$_}\n";

    ## fetch proteins and print to temp file
    open (my $PRO, ">temp.faa") or die $!;
    if ((exists($protein_hash{$_})) && (exists($protein_hash{$pairs{$_}}))) {
      print $PRO ">$_\n$protein_hash{$_}\n>$pairs{$_}\n$protein_hash{$pairs{$_}}";
      close $PRO;
    } else {
      die "[ERROR] Protein ID '$_' or '$pairs{$_}' not found in file $proteinfile!\n";
    }

    ## make CDS hash of nucleotides
    my %cds_seqs;
    if ((exists($cds_hash{$_})) && (exists($cds_hash{$pairs{$_}}))) {
      $cds_seqs{$_} = Bio::Seq->new( -display_id => "$_", -seq => $cds_hash{$_} );
      $cds_seqs{$pairs{$_}} = Bio::Seq->new( -display_id => "$pairs{$_}", -seq => $cds_hash{$pairs{$_}} );
    } else {
      die "[ERROR] CDS ID '$_' or '$pairs{$_}' not found in file $cdsfile!\n";
    }

    ## run alignment
    if (system ("clustalo --infile=temp.faa --outfile=temp.aln --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }

    ## fetch alignment, backtranslate to nucleotides & write
    my $get_prot_aln = Bio::AlignIO -> new(-file=>"temp.aln", -format=>"fasta");
    my $prot_aln = $get_prot_aln -> next_aln();
    my $dna_aln = aa_to_dna_aln($prot_aln, \%cds_seqs);
    my $write_dna_aln = Bio::AlignIO -> new(-file=>">$out.$n.fna", -format=>"fasta");
    $write_dna_aln->write_aln($dna_aln);
    $n++;
  }
  close $OUT;
}

##################################################
## alignment of whole gene region (intron+exons) parsed from gff
##################################################

elsif ($fastafile && $gfffile) {
  print STDERR "[INFO] Will run MAFFT alignment on gene regions\n";
  print STDERR "[INFO] Parsing genome fasta and GFF file...\n";
  ## parse genome fasta and GFF:
  my (%genome_hash, %gff_hash);
  my $in_f = Bio::SeqIO->new( -file => $fastafile, -format => 'fasta' );
  while (my $seq_obj = $in_f->next_seq() ) {
    $genome_hash{$seq_obj->display_id()} = $seq_obj; ##key=seqid; val=seq_obj
  }
  die "[ERROR] No sequences found in $fastafile!\n\n" if scalar(keys %genome_hash) == 0;
  print STDERR "[INFO] Fetched ".commify(scalar(keys %genome_hash))." sequences from $fastafile\n";
  open (my $GFF, $gfffile) or die "[ERROR] Cannot open $gfffile: $!\n\n";
  while (<$GFF>) {
    chomp;
    my @F = split (/\s+/, $_);
    if ($F[2] eq "gene") {
      if ($F[8]=~/ID\=(.+)\;/) {
        ## NOTE .t1 suffix added to hash key!!!
        $gff_hash{"$1.t1"} = { 'scaffold'=>$F[0], 'strand'=>$F[6], 'start'=>$F[3], 'end'=>$F[4] } ##HoH key=geneid; val={scaffold,strand,start,end}
      }
    }
  }
  close $GFF;
  die "[ERROR] No genes found in $gfffile!\n\n" if scalar(keys %gff_hash) == 0;
  print STDERR "[INFO] Fetched ".commify(scalar(keys %gff_hash))." genes from $gfffile\n";

  ## generate alignments for all pairs in %pairs:
  my $n = 1;
  open (my $OUT, ">$infile.map") or die "[ERROR] Cannot open $infile.map: $!\n\n";
  print $OUT "NUM\tSEQ1\tSEQ2\n";
  foreach (nsort keys %pairs) {
    print STDERR "\r[INFO] Working on pair \#$n: $_, $pairs{$_}";$| = 1;
    print $OUT "$n\t$_\t$pairs{$_}\n";

    ## get sequence hashes:
    my %h1 = %{$gff_hash{$_}};
    my %h2 = %{$gff_hash{$pairs{$_}}};
    my $seqobj1 = Bio::Seq->new(-display_id=>$_, -seq=>$genome_hash{$h1{scaffold}}->subseq($h1{start},$h1{end})); ##returns a seqobj
    my $seqobj2 = Bio::Seq->new(-display_id=>$pairs{$_}, -seq=>$genome_hash{$h2{scaffold}}->subseq($h2{start},$h2{end}));
    $seqobj1 = $seqobj1->revcom() if $h1{strand} eq "-"; ##revcom if on the minus strand
    $seqobj2 = $seqobj2->revcom() if $h2{strand} eq "-"; ##revcom if on the minus strand
    die "[ERROR] Empty sequence! Something's up...\n\n" if (($seqobj1->length()==0) || ($seqobj2->length()==0));
    ## print to temporary fna file:
    open (my $NUC, ">temp.fna") or die "[ERROR] Cannot write to temp fasta: $!\n\n";
    print $NUC ">$_ [$h1{scaffold}|$h1{strand}|$h1{start}|$h1{end}]\n".$seqobj1->seq()."\n";
    print $NUC ">$pairs{$_} [$h2{scaffold}|$h2{strand}|$h2{start}|$h2{end}]\n".$seqobj2->seq()."\n";
    close $NUC;

    ## run alignment:
    if ( system("mafft --auto --quiet --thread $threads temp.fna > $out.$n.fna") != 0 ) { die "[ERROR] Problem running MAFFT!\n\n"; }
    $n++;
  }
  close $OUT;
}

## file cleanup:
if (-d $outdir) {
  rmtree([ "$outdir" ]);
  mkdir $outdir;
  system("mv $out* $outdir");
} else {
  mkdir $outdir;
  system("mv $out* $outdir");
}
system ("rm temp.*");
print "\n[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
