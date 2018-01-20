#! /usr/bin/perl

#goal of script:
#create specialized clade list for each set of SNP reports against each reference.


use strict;
use warnings;
use Data::Dumper;  #for debugging
use List::Util qw(min max sum);  
use Getopt::Long;



	my $directory;
	my $outfile = '_clade_list.txt';
	my $cladeListReference = 'masterCladeList.txt';
	my @fileList;
	my %cladetracker;

GetOptions( "d|directory=s" => \$directory);

if (!$directory) {die "\nno directory  set with -d option.\n";  }




#list files in directory
    opendir(DIR, $directory) or die $!;

    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);#ignore hidden files
        push @fileList, $file;
    }
    closedir(DIR);

#read in master clade list

open(my $fh, "<", $cladeListReference)
or die "couldnt open '$cladeListReference' $!";		
	while (<$fh>){
	chomp;
	my @split = split(/\t/);
	my $name = $split[1];
	my $clade = $split[7];
	#$name = $directory.$name;
	$cladetracker{$name} = $clade;
		}
close $fh;

#create outfile handle
my $sample= $fileList[1];
$sample =~s/_.*//;

#Open outfile
 open (my $oh, '>', $sample.$outfile) or die "could not open $outfile $!";	
#match 
#result will be output file -> clade

foreach my $thisFile (@fileList) { 
	#remove ref tag
	my $thisFileFull = $thisFile;
	$thisFile =~s/.*_v_//;#remove first half
	$thisFile =~s/_.*//;#remove second half
	#search for matching key
	my $thisClade = $cladetracker{$thisFile};
	#print out original name with clade match
	if (defined $thisClade){
		print $oh "$thisFileFull\t$thisClade\n";
		}else { 
		 print $oh "$thisFileFull\tUndefined\n"; } 
	} 

close($oh);


