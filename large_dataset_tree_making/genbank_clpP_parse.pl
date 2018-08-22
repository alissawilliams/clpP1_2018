#!/usr/bin/perl

use warnings;
use strict;
use sloan;
use Bio::SeqIO;


my $usage = "\nUSAGE: $0 (multi)GenBank_file output\n\n";

my $file = shift or die ($usage);
my $output = shift or die ($usage);
my $fasta = open_output ("$output\.fas");
my $summary = open_output ("$output\.summary.txt");

print $summary "Species\tAccession\tExtracted_Seqs\n";

my $seqio_object = Bio::SeqIO->new(-file => $file);

while (my $seq_object = $seqio_object->next_seq){
	my $species = $seq_object->species->genus . " " . $seq_object->species->species;
	my $accession = $seq_object->accession_number();
	my $count = 0;
	for my $feat_object ($seq_object->get_SeqFeatures){
		if ($feat_object->primary_tag eq "CDS") {
			if ($feat_object->has_tag('gene')) {
				for my $val ($feat_object->get_tag_values('gene')){
					if ($val =~ /clp/i){
	      				++$count;
	      				my $seq = $feat_object->spliced_seq->seq;
	      				my $exonCount = 0;
						if ( $feat_object->location->isa('Bio::Location::SplitLocationI'))  {
							for ( $feat_object->location->sub_Location ) {
								++$exonCount	
							}
						}else{
							$exonCount = 1;
						}
	      				print $fasta ">$species\t$accession\t$val\t$count\texons=$exonCount\n$seq\n";
	      			}
	      		}
	      	}
	   	}
	}
	print $summary "$species\t$accession\t$count\n";
}

exit;