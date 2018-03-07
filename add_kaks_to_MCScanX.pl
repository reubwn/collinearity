#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Seq::EncodedSeq;
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);

local $SIG{__WARN__} = sub { warn $_[0] unless $_[0] =~ /WARNING/}; ##suppress WARNINGs from BioPerl

my $usage = "
SYNOPSIS
  Calculates Ka and Ks for pairs of genes from an MCScanX collinearity file.
  Will detect CDS which are non-multiples of 3 and trim from 3' accordingly (script
  which ships with MCScanX does not do this and can lead to errors..)

OPTIONS
  -i|--in      [FILE] : MCScanX collinearity file
  -p|--prot    [FILE] : fasta file of protein sequences
  -c|--cds     [FILE] : fasta file of corresponding CDS (nucleotide)
  -t|--threads [INT]  : number of aligner threads (default = 1)
  -h|--help           : this message

OUTPUTS
  Annotated MCScanX collinearity file, with Ka and Ks as final 2 columns for every gene pair.
  Will write '-2' for cases where Ka/Ks cannot be calculated.
\n";

my ($infile, $proteinfile, $cdsfile, $help);
my $threads = 1;

GetOptions (
  'i|in=s'      => \$infile,
  'p|prot=s'    => \$proteinfile,
  'c|cds=s'     => \$cdsfile,
  't|threads:i' => \$threads,
  'h|help'      => \$help
);

die $usage if $help;
die $usage unless ($infile && $proteinfile && $cdsfile );
my $outfile = "$infile.kaks";
my $warningsfile = "$infile.kaks.warnings";

print STDERR "[INFO] Collinearity file: $infile\n";
print STDERR "[INFO] Proteins file: $proteinfile\n";
print STDERR "[INFO] CDS file: $cdsfile\n";
print STDERR "[INFO] Parsing protein and CDS files...\n";

## parse proteins and CDS:
my (%protein_hash, %cds_hash);
my $in_prot = Bio::SeqIO->new( -file => $proteinfile, -format => 'fasta' );
while (my $seq = $in_prot->next_seq() ) {
  $protein_hash{$seq->display_id()} = $seq->seq();
}
die "[ERROR] No sequences found in $proteinfile!\n" if (scalar(keys %protein_hash) == 0);
print "[INFO] Fetched ".commify(scalar(keys %protein_hash))." proteins from $proteinfile\n";
my $in_cds = Bio::SeqIO->new( -file => $cdsfile, -format => 'fasta' );
while (my $seq = $in_cds->next_seq() ) {
  $cds_hash{$seq->display_id()} = $seq->seq();
}
die "[ERROR] No sequences found in $cdsfile!\n" if (scalar(keys %cds_hash) == 0);
print "[INFO] Fetched ".commify(scalar(keys %cds_hash))." CDS from $cdsfile\n";

print STDERR "[INFO] Parsing collinearity file...\n";

## parse collinearity file:
open (my $OUT, ">$outfile") or die "Cannot open $outfile: $!\n\n";
open (my $WARN, ">$warningsfile") or die "Cannot open $warningsfile: $!\n\n";
open (my $IN, $infile) or die "Cannot open $infile: $!\n\n";
my ($ntrimmed,$na,$nskipped) = (0,0,0);
my $n = 1;
PAIR: while (<$IN>) {
  if ($_ =~ m/^\#/) {
    print $OUT $_;
  } else {
    chomp;
    my $LINE = $_;
    $LINE =~ s/^\s+|\s+$//g; ##remove leading and trailing whitespaces
    my @F = split (/\s+/, $LINE);
    #print STDERR "$LINE\t$F[-3]\t$F[-2]\n";
    my $gene1 = $F[-3]; ##work from end of array as columns in collinearity file not consistently formatted
    my $gene2 = $F[-2];

    print STDERR "\r[INFO] Working on pair \#$n: $gene1, $gene2";$| = 1;

    ## fetch proteins and print to temp file:
    open (my $PRO, ">temp.faa") or die $!;
    if ((exists($protein_hash{$gene1})) && (exists($protein_hash{$gene2}))) {
      print $PRO ">$gene1\n$protein_hash{$gene1}\n>$gene2\n$protein_hash{$gene2}";
      close $PRO;
    } else {
      die "[ERROR] Protein ID '$gene1' or '$gene2' not found in file $proteinfile!\n";
    }

    ## run clustalo alignment:
    if (system ("clustalo --infile=temp.faa --outfile=temp.aln --force --threads=$threads") != 0) { die "[ERROR] Problem with clustalo!\n"; }

    ## make CDS hash of nucleotides:
    my %cds_seqs;
    if ((exists($cds_hash{$gene1})) && (exists($cds_hash{$gene2}))) {
      $cds_seqs{"$gene1"} = Bio::Seq->new( -display_id => "$gene1", -seq => $cds_hash{$gene1} );
      $cds_seqs{"$gene2"} = Bio::Seq->new( -display_id => "$gene2", -seq => $cds_hash{$gene2} );
    } else {
      die "[ERROR] CDS '$gene1' or '$gene2' not found in file $cdsfile!\n";
    }
    ## check if all CDS seqs are multiple of 3:
    foreach (keys %cds_seqs) {
      if ($cds_seqs{$_}->length() % 3 != 0) {
        print $WARN "[WARN] Seq ".$cds_seqs{$_}->display_id()." is not a multiple of 3\n";
        print $WARN "[INFO] Length: ".$cds_seqs{$_}->length()."; Excess: ".($cds_seqs{$_}->length() % 3)."\n";
        print $WARN "[INFO] Translation in current frame: ".($cds_seqs{$_}->translate(-frame=>0, -complete=>1)->seq()."\n");
        ## determine frame of sequence:
        my ($frame, $trimmed);
        if ( $cds_seqs{$_}->translate(-frame=>0, -complete=>1)->seq() !~ m/\*/g ) { ##look for INTERNAL stop codons; -complete should trim the terminal codon from all alignments, if present
          ## seq is in frame, therefore errant codon must be at the 3' end:
          $trimmed = $cds_seqs{$_}->subseq(1,(($cds_seqs{$_}->length()) - ($cds_seqs{$_}->length() % 3))); ##trims remainder off 3' end; returns a STRING, annoyingly
          $frame = 0;
        } elsif ( $cds_seqs{$_}->translate(-frame=>1, -complete=>1)->seq() !~ m/\*/g ) {
          ## seq is OUT of frame (offset 1); therefore trim 1 base from 5' end:
          $trimmed = $cds_seqs{$_}->subseq(2,(($cds_seqs{$_}->length()))); ##NB subseq is coord inclusive
          $frame = 1;
        } elsif ( $cds_seqs{$_}->translate(-frame=>2, -complete=>1)->seq() !~ m/\*/g ) {
          ## seq is OUT of frame (offset 2); therefore trim 2 base from 5' end:
          $trimmed = $cds_seqs{$_}->subseq(3,(($cds_seqs{$_}->length())));
          $frame = 2;
        } else {
          ## there are stops in all 3 frames (?); some kind of error/pseudogene; just skip it
          print $WARN "[ERROR] Frame 0 Translation: ".$cds_seqs{$_}->translate(-frame=>0, -complete=>1)->seq()."\n";
          print $WARN "[ERROR] Frame 1 Translation: ".$cds_seqs{$_}->translate(-frame=>1, -complete=>1)->seq()."\n";
          print $WARN "[ERROR] Frame 2 Translation: ".$cds_seqs{$_}->translate(-frame=>2, -complete=>1)->seq()."\n";
          print $WARN "[ERROR] Stop codons detected in all 3 frames! Skipping this pair...\n";
          $nskipped++;
          next PAIR; ##skip to next pair
        }
        ## replace old seq with trimmed seq:
        $cds_seqs{$_} = Bio::Seq->new( -display_id => $_, -seq => $trimmed );
        $ntrimmed++;

        ## print some info to WARN:
        print $WARN "[INFO] Correct frame: $frame\n";
        print $WARN "[INFO] Translation in correct frame: ".($cds_seqs{$_}->translate(-frame=>0, -complete=>1)->seq()."\n");
      }
    }

    ## fetch alignment, backtranslate to nucleotides:
    my $get_prot_aln = Bio::AlignIO -> new(-file=>"temp.aln", -format=>"fasta");
    my $prot_aln = $get_prot_aln -> next_aln();
    my $dna_aln = aa_to_dna_aln($prot_aln, \%cds_seqs);

    ## get Ka (Dn), Ks (Ds) values for the pair of genes:
    my ($Ka,$Ks) = (-2,-2); ## default values
    eval {
      my $stats = Bio::Align::DNAStatistics->new();
      my @result = @{ $stats->calc_KaKs_pair($dna_aln,$gene1,$gene2) };
      if (exists $result[0]{'D_n'}) {
        $Ka = $result[0]{'D_n'};
      }
      if (exists $result[0]{'D_s'}) {
        $Ks = $result[0]{'D_s'};
      }
    };
    $na++ if ( ($Ka == -2) || ($Ks == -2) );
    print $OUT "$_\t$Ka\t$Ks\n";


    # eval {
    #   my $stats = Bio::Align::DNAStatistics->new();
    #   my $result = $stats->calc_all_KaKs_pairs($dna_aln);
    #   my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
    #   for my $an (@$result) {
    #     for (sort keys %$an ) {
    #       next if /Seq/;
    #       if($_ eq "D_n"){$Dn = $an->{$_}};
    #       if($_ eq "D_s"){$Ds = $an->{$_}};
    #       if($_ eq "S_d"){$S_d = $an->{$_};}
    #       if($_ eq "N_d"){$N_d = $an->{$_};}
    #       if($_ eq "S"){$S = $an->{$_};}
    #       if($_ eq "N"){$N = $an->{$_};}
    #     }
    #   }
    #   $Dn = -2 unless ($Dn); ##default values
    #   $Ds = -2 unless ($Ds);
    #   $na++ if ( ($Dn == -2) || ($Ds == -2) );
    #
    #   ## annotate input file
    #   print $OUT "$_\t$Dn\t$Ds\n";
    # };
    $n++;
  }
}
close $IN;
close $OUT;
close $WARN;
system ("rm temp.*"); ##remove last temp files.

print STDERR "\n";
print STDERR "[INFO] Total number of pairs: $n\n";
print STDERR "[INFO] Number of seqs trimmed because length % 3 != 0: $ntrimmed\n";
print STDERR "[INFO] Number of pairs skipped because stop codons detected in all 3 +frames: $nskipped\n";
print STDERR "[INFO] Number of pairs for which Ka or Ks was not calculated: $na\n";
print STDERR "[INFO] Finished on ".`date`."\n";

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
