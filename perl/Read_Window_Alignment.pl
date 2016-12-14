#
####
#written by Bradford Condon, Dec 13, 2016
#Goal of script:
#Read in report file on number of nucleotides aligned for each scaffold
####
####
#Input format
#Query: B71_SZL_masked	Subject: 87-120_masked
#Identifying repeat regions
#B71_SZL_masked.87-120_masked.BLAST
#chr1	0	0
#above are line types 1, 2, 3, and 4.  want line type 4.



use strict;
use warnings;
use Data::Dumper qw(Dumper);

sub ReadWindowAlignment{ 
my $filename  = shift;
my %outputArray;

open(my $fh, "<", $filename) or die "Could not open file '$filename' $!";

my $ref;
my $quer;
while (my $line = <$fh>) {
	chomp $line;
	if ($line =~ /Query/) { #check if line type 1
		my @split = split(' ', $line);
		$ref = $split[1];
		$quer = $split[3];
		$quer =~ s/' '//;
		$ref =~ s/' '//;
		next;
	}
	if ($line =~ /Identifying/) {  #check if line type 2, skip
		next;
	}

	#split line
	my @split = split("\t", $line);

	if (!exists $split[1]){  		#skip line type 3
		next;
	}
	if (not defined $ref) {
		die "error: please check input format and header.\nFirst line should resemble:\nQuery: B71_SZL_masked	Subject: 87-120_masked\n";
	}

	#trim chr
	$split[0]=~s/chr//;

if ($split[0] =~ /ctg/) {
	$split[0]=~s/ctg//;
	$split[0] = $split[0]+ 7; 	#add 7 to it
	}

	#format is scaffold, start, numberOfNucsAligned
	$outputArray{$quer}{$split[0]}{$split[1]} = $split[2];
}
return %outputArray;
}



1;