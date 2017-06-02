#!/usr/bin/env perl

use strict;
use warnings;
use Sort::Naturally;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $usage = "
SYNOPSIS:
  Fetches and prints all connections for an initial input scaffold name.
  Use --all flag to also print all downstream subconnections.

OPTIONS:
  -i|--collinearity [FILE]   : collinearity file from MCScanX
  -s|--search       [STRING] : chromosome / scaffold name to get connections for
  -a|--all                   : fetch all downstream connections as well
  -h|--help                  : print this message

USAGE:
  >> get_all_connections_from_collinarity.pl -i xyz.collinearity -s Ar1 --all
\n";

my ($collinearity, $search, $all, $help);

GetOptions (
  'collinearity|i=s' => \$collinearity,
  'search|s=s'       => \$search,
  'all|a'            => \$all,
  'help|h'           => \$help,
);

die $usage if $help;
die $usage unless ($collinearity && $search);

my %connections;
open (my $IN, $collinearity) or die $!;

while (<$IN>) {
  if ($_ =~ m/(\w+\d+)\&(\w+\d+)\s/) {
    push @{ $connections{$1} }, $2;
    push @{ $connections{$2} }, $1;
  }
}

chomp $search;
my @result = @{$connections{$search}};
push @result, $search; ## include the focal chromosome!

## fetch all downstream subconnections; ie any scaffold that shares a collinear region with a scaffold in original @result
if ($all) {
  foreach ( nsort @{$connections{$search}} ) {
    push @result, @{$connections{$_}};
  }
}

my @uniq_result = uniq(@result);
print join (",", (nsort @uniq_result)), "\n";

__END__
