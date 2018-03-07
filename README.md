# collinearity
Scripts to parse and analyse MCScanX collinearity output.

### doi

[![DOI](https://zenodo.org/badge/92963110.svg)](https://zenodo.org/badge/latestdoi/92963110)

## add_kaks_to_MCScanX.pl
This script annotates the Xyz.collinearity output file from MCScanX analysis with Ka and Ks values for all pairs of genes.

### Prerequisites
1. Requires ClustalO to be in `$PATH`
2. `BioPerl` libraries

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
1. Annotated MCScanX collinearity file
2. Annotated MCScanX collinearity file also reformatted to remove some of the strange leading whitespaces etc.

## calculate_collinarity_metric.pl
Calculates collinearity score based on the number of collinear genes divided by the total number of genes within that defined block (from [Flot et al. 2013](http://dx.doi.org/10.1038/nature12326)).
Also calculates average Ka and Ks per block, if run on annotated MCScanX file.

## Options
```
OPTIONS
  -i|--collinearity [FILE] : collinearity file from MCScanX
  -g|--gff          [FILE] : modified gff file from MCScanX
  -k|--kaks                : parse collinearity file to get average ka & ks per block
  -h|--help                : print this message
```

##
