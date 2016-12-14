#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
# usage
my $usage = "$0 -i|input fasta_input_files, -m|map print map of mask locations (Default = FALSE) \n";
# global values
my $input_file;
my $mapCheck;
# read user options
GetOptions(
	"i|input=s" => \$input_file,
	"m|map=s" => \$mapCheck
);
# check for user parameters, set defaults
if( !$input_file ){
	die $usage;
}
# open input fasta file
unless( -e $input_file ){
	die "error: cannot find fasta input file $input_file\n";
}
my $input = Bio::SeqIO->new (-file => "<$input_file", '-format' => 'Fasta')
	or die "error: failure opening fasta $input_file for reading: $!\n";
my $prefix = \$input_file;
my $mask_file = "$prefix", ".sizemaskcount.txt";
my $mask_map_file = "$prefix", ".map.txt";
open (my $mask_summary, '>', $mask_file) or die "could not open $mask_file $!";
#print header
print $mask_summary "ID\tLength\tMask_count\n";
# step through sequences in input fasta file
while (my $seqObject = $input->next_seq){
	
	# get sequence information
	my $id  = $seqObject->id;
	my $seq = $seqObject->seq;
	my @seqarray = split('', $seq);
	my $length = $seqObject->length;
 	my $maskcount=0;
 	my $easycount=0;

		for (@seqarray) {   #look through  sequence
			++$easycount;
						if (($_ eq 'X')||($_ eq 'x')||($_ eq 'N')||($_ eq 'n')) {  #if the nucleotide is an n or an x, add it to the mask count
							$maskcount++;						
						}
		}
		print $mask_summary "$id\t$length\t$maskcount\n";
}
close $mask_summary;
