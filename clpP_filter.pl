#!/usr/bin/perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: $0 fasta_file accessionList_file\n\n";

my $fasta_file = shift or die ($usage);
my $accession_file = shift or die ($usage);

my $FH = open_file ($accession_file);

my %accs;

while (<$FH>){
	chomp $_;
	$accs{$_} = 1;
}

my %fasta = fasta2hash($fasta_file);

foreach (sort keys %fasta){
	my @sk = split (/\t/, $_);
	exists ($accs{$sk[1]}) or next;
	my $prot = dna2peptide ($fasta{$_});
	substr($prot, -1) eq '*' or print STDERR "Does not end in stop codon: $_\n";
	substr($prot, 0, -1) =~ /\*/ and print STDERR "Internal stop codon(s): $_\n";
	print ">$_\n$fasta{$_}\n";
}

exit;