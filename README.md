# collinearity
Scripts to parse and analyse MCScanX collinearity output.

### doi

[![DOI](https://zenodo.org/badge/92963110.svg)](https://zenodo.org/badge/latestdoi/92963110)

## MCScanX pipeline
A typical MCScanX analysis following gene finding using Augustus or Braker might consist of:

1. Run Diamond:
   ```
   >> diamond makedb --in augustus.aa -d augustus.aa
   >> diamond blastp -e 1e-5 -p 8 -q augustus.aa -d augustus.aa -a augustus.aa.vs.self
   >> diamond view -a augustus.aa.vs.self.daa -o Xyz.blast
   ```
2. Generate CDS (required for Ka/Ks calculation):
   ```
   >> /path/to/augustus-3.2.1/scripts/getAnnoFasta.pl --seqfile=mygenome.fasta augustus.gff
   ```
3. Generate GFF:
   ```
   >> perl -lane 'print join("\t",$F[0],$F[8],$F[3],$F[4]) if ($F[2]eq"transcript")' augustus.gff > Xyz.gff
   ```
4. Run MCScanX:
   ```
   >> mkdir results
   >> mv Xyz* results/
   >> /path/to/MCScanX/MCScanX results/Xyz
   ```

---

## add_kaks_to_MCScanX.pl
This script annotates the Xyz.collinearity output file from MCScanX analysis with Ka and Ks values for all pairs of genes.

### Prerequisites
Requires [ClustalO](http://www.clustal.org/omega/) to be in `$PATH` and some BioPerl libraries.

### Options
Type `add_kaks_to_MCScanX.pl -h` to see the options:
```
OPTIONS
  -i|--in      [FILE] : MCScanX collinearity file
  -p|--prot    [FILE] : fasta file of protein sequences
  -c|--cds     [FILE] : fasta file of corresponding CDS (nucleotide)
  -t|--threads [INT]  : number of aligner threads (default = 1)
  -h|--help           : this message
```

### Outputs
MCScanX collinearity file annotated with Ka and Ks per pair of genes.

## calculate_collinarity_metric.pl
Calculates collinearity score based on the number of collinear genes divided by the total number of genes within that defined block (from [Flot et al. 2013](http://dx.doi.org/10.1038/nature12326)).
Also calculates average Ka and Ks per block, if run on annotated MCScanX file.

### Options
```
  -i|--in   [FILE] : collinearity file from MCScanX
  -g|--gff  [FILE] : modified gff file from MCScanX
  -k|--kaks        : parse collinearity file to get average ka & ks per block
  -h|--help        : print this message
```

### Outputs
A score file with per-block collinearity scores, average Ka and Ks etc.
Also writes an MCScanX collinearity file reformatted to remove some of the strange leading whitespaces etc.

## calculate_collinarity_breakpoints.pl
Searches for breaks in collinearity defined as occurrences where homologous blocks cannot be aligned along scaffolds without some rearrangement.
Also counts cases where homologous regions are found physically linked on the same scaffold.

Type `calculate_collinarity_breakpoints.pl -m` to see full help and definitions.

### Options
```
  -i|--in     [FILE]  : collinearity file (use .reformatted!)
  -g|--gff    [FILE]  : GFF file
  -s|--score  [FILE]  : score file with average Ks per block
  -k|--ks     [FLOAT] : Ks threshold to define homologous block (default <= 0.5)
  -o|--out            : outfile (default=INFILE.breaks)
  -b|--blocks         : also print blocks_per_gene file (default=no)
  -h|--help           : print help
  -m|--morehelp       : print more help
```

### Definition of a breakpoint
A collinearity break is introduced if there is a mistmatch in the identity of the upstream or downstream homologous blocks from one collinear region to another. Thus, if `a[i]` is the identity of the focal block on chrom `a`, homologous with block `b[j]` on chrom `b`, collinearity is broken if `a[i+/-1] != b[j+/-1]`. A caveat is when the focal block is located at the terminus of a scaffold; here, collinearity is not broken if the homologous scaffold can be orientated 'away' from the focal scaffold.

1. Collinearity cannot be broken by scaffolds containing only 1 collinear block.
2. Collinearity cannot be broken when focal and subject blocks are BOTH terminal on their respective scaffolds.

#### Examples
```
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
```

## calculate_collinarity_palindromes.pl
Searched for cases of homologous regions physically linked on the same scaffolds.
Decomposes any such cases into tandem arrays or palindromes.

### Options
```
-i|--in   [FILE] : collinearity file from MCScanX (annotated with Ks)
-g|--gff  [FILE] : GFF file
-k|--ks          : only examine blocks with Ks <= this threshold
-h|--help        : print this message
```

### Output
1. arrays file: all results
2. arrays.tandem: tandem arrays
3. arrays.palindrome: palindromic arrays

## Other Scripts
There are a bunch of other helper scripts for Circos plotting... type `-h` for options.
