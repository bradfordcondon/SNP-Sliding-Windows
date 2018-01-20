 #! /usr/bin/perl
	#this script goes through SNP output.  It counts the number of times a scaffold has SNPs compared to another genome.
	#updated 11-18-16
	#Main goals of this script:
	#For every given window, determine which genome has the fewest SNPs (from an outclade)
	#to do this, build a hash with number of SNPs for each window for each genome.
	#then, determine which clade the least SNPs are from (excluding whatever is set as the in taxon)

	use strict;
	use warnings;
	use Data::Dumper;  #for debugging
	use List::Util qw(min max sum);  
	use Getopt::Long;



	my $usage = "$0 -d|directory directory -p|conserved TRUE/FALSE generate conserved region file \n 
	-c|cladelist 3 column clade list -r TRUE/FALSE generate report file -i mummer (default = FALSE) \n";

#
#Option variable

	my $directory = '/Users/chet/uky/SNP_density_windows/12-21-16-test/reports/';
	my $cladelist;
	my $windowsize = '1000';
	my $stepsizemaster = '1000';
	my $minsnpsite = '0';
	my $mummer = 'FALSE';
	my $maxPrimerSnp = '400';
	my $lengthDataFile = '/Users/chet/uky/SNP_density_windows/12-27-16_full/B71_SZL_lengthlist_named.txt';
	my $outfile = 'outCladeWindows.txt';
	my $outfile2 = 'averageWindowSNPsbyClade.txt';
	my $windowCoverageHashDirectory;# '/Users/chet/uky/SNP_density_windows/12-26-16-test/B71_SZL_windows/'
	my $cladeListReference= '/Users/chet/uky/SNP_density_windows/12-26-16-test/masterCladeList.txt';


#/Users/chet/uky/SNP_density_windows/12-21-16-test
#
#Internal variables
#

	my %reporthash;	

	
	# read user options
	GetOptions(
		"d|directory=s" => \$directory,
		"p|coverageDirectory=s" => \$windowCoverageHashDirectory,
		"w|windowsize=s" =>\$windowsize,
		"s|stepsize=s" =>\$stepsizemaster,
		"m|minsnpsite=s" =>\$minsnpsite,
		"i|mummer=s" =>\$mummer,
		"x|maxPrimerSNP=s" =>\$maxPrimerSnp,
		"e|lengthFile=s" =>\$lengthDataFile,
		
			);

		# check for user parameters	

if (!$lengthDataFile) {die "no length data file set -e.  \n";  }

	my %masterhash;  #declare global variables. 
	my %allhash;
	my %allnucs;

	my @files = <$directory*>;  #separate reference files
	
######
#Load in file of window alignment report
######
print "reading in Report file for window alignments\n";
#my %windowCoverageHash = ReadWindowAlignment($windowCoverageHashFile);


my @reportFiles = <$windowCoverageHashDirectory*>;

my %windowCoverageHash;

foreach my $file (@reportFiles) {   ##open each SNP file
	print "working on $file \n!";

ReadWindowAlignmentSeperateFile($file);

			}	



my %inverseCladeTracker;
my %cladetracker;
my %masterCladeTracker;
#########
#Read in master clade list.  Then, search directory to assign
#
#########

open(my $fh, "<", $cladeListReference)

or die "couldnt open '$cladeListReference' $!";		
	while (<$fh>){
	chomp;
	my @split = split(/\t/);
	my $name = $split[1];
	my $clade = $split[7];
	#$name = $directory.$name;
	$masterCladeTracker{$name} = $clade;
		}
close $fh;


my @fileList;


#list files in directory
    opendir(DIR, $directory) or die $!;
    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);#ignore hidden files
	my $thisFileFull = $file;

	$file =~s/.*_v_//;#remove first half
	$file =~s/_.*//;#remove second half
	#search for matching key
	my $thisClade = $masterCladeTracker{$file};
	#print out original name with clade match
	if (defined $thisClade){
		print  "$file\t$thisClade\n";
		$cladetracker{$file} = $thisClade;
	$inverseCladeTracker{$thisClade}{$file} = '1';					
	}else { 
		 print  "$thisFileFull\tUndefined\n";  
		$cladetracker{$file} = "Undefined";
	$inverseCladeTracker{"Undefined"}{$file} = '1';					

	} 

}
   
   closedir(DIR);





 ###################################################                	   
	####	Go through each SNP report and build a hash
 ###################################################                	   	
	foreach my $file (@files) {   ##open each SNP file
	print "working on $file \n!";
######################
	#### Build hash of each nucleotide change at the SNP
######################			
		#%allnucs = build_SNP_NUC_hash($file, \%allnucs);
	   	#%allhash = build_SNP_hash($file, \%allhash);	 	 				 		
		%allhash = build_Clade_Hash($file, \%allhash); #builds a hash structured by scaffold, nuc number, FILE
			}	
##Check that allhash is filled.
#If not, user likely set the wrong input option

 ###################################################                	   
 ###################################################                	   


###################################################                	   
	#### Load in list of scaffold lengths.
###################################################    

open(my $length_file, "<", $lengthDataFile)
or die "couldn't open length data file\n";

my %lengthinfo;
while (<$length_file>){
	chomp;
	my @split = split(/\s+/);
	my $key = $split[0];
	my $value = $split[1];
	$lengthinfo{$key} = $value;
}

close $length_file;

#test that $lengthinfo has more than 1 key

 if (keys %lengthinfo <= 1){
 	die "$length_file only had a single entry.  Do non-UNIX characters need to be stripped?\n";
 }

		
	print "Calculating SNPs in windows\n";

#step size then windowsize
	#my %SNP_rich_windows= Diversity_sliding_window(\%allhash, \%allhash, $stepsizemaster, $windowsize);	  
	my ($allWindowsByClade, $allWindowsAltArrange) = Clade_sliding_window(\%allhash, $stepsizemaster, $windowsize); #hash, step size, windowsize
####
my %allWindowsByClade = %$allWindowsByClade;
my %allWindowsAltArrange = %$allWindowsAltArrange;


####
#Print out ALL snp values for every taxon.
#Currently commented out in favor of only printing average and minimum for the clade
####
	
# open (my $oh, '>', $outfile) or die "could not open $outfile $!";	
# 	print "outputting SNP counts for each window by taxon\n";
# 	print $oh "Scaffold\tstart\tend\tTaxon\tClade\tnumSNPs\n";
# 	foreach my $thisScaffold (sort keys %allWindowsByClade) {
#   		foreach my $thisTaxon( sort keys $allWindowsByClade{$thisScaffold}) { 
#   				my $thisCleanTaxon = CleanTaxonName($thisTaxon);
#   			foreach my $start (sort keys $allWindowsByClade{$thisScaffold}{$thisTaxon}){
# 					my $snpCount =$allWindowsByClade{$thisScaffold}{$thisTaxon}{$start}{"count"};
#   					my $end = 	$allWindowsByClade{$thisScaffold}{$thisTaxon}{$start}{"end"};
# 					print $oh "$thisScaffold\t$start\t$end\t$thisCleanTaxon\t$snpCount\n"; ###print out keys at end  	
#   			}
# 		}
# 	}			
#close $oh;
	




#determine the best taxon for each
my $AvgWindowVals = DetermineBestClade(\%allWindowsAltArrange); #hash, step size, windowsize
my %AvgWindowVals = %$AvgWindowVals;


		
open (my $oh2, '>', $outfile2) or die "could not open $outfile2 $!";	
	print $oh2 "Scaffold\tstart\tClade\tavgnumSNPs\tminSNPs\n";
	foreach my $thisScaffold (sort keys %AvgWindowVals) {
  		foreach my $thisStart(sort  keys $AvgWindowVals{$thisScaffold} ) { 
  			foreach my $thisClade (keys $AvgWindowVals{$thisScaffold}{$thisStart}){#error here
  					my $minimum;
					my $average = $AvgWindowVals{$thisScaffold}{$thisStart}{$thisClade}{"avg"};
  					if (defined  $AvgWindowVals{$thisScaffold}{$thisStart}{$thisClade}{"min"}){
  					 	$minimum = 	$AvgWindowVals{$thisScaffold}{$thisStart}{$thisClade}{"min"};  						
  					}else{$minimum = 0;}
					print $oh2 "$thisScaffold\t$thisStart\t$thisClade\t$average\t$minimum\n"; ###print out keys at end  	
  			}
  			
		}
	}			
close $oh2;
	



	###Subroutines######
			
			
	########
	########
	#Sliding window subroutine.  
	#This Subroutine takes two 2-layer hash.  It sorts them by scaffolds, and builds an array of nucleotides for each scaffold.
	#It  then slides along at a defined stepsize and window size.  It scrolls through the first (%hash) and checks against teh second %hashcheck
	#If there is a value in the input hash found as it slides along, it counts it.
	#It will output a 2 level hash: scaffold, range,  and count.
			
			##Input is hash reference, step size, windowsize.
			
			
				sub sliding_window {  
				
   				 my %hash = %{shift()};
   				 my %hashcheck = %{shift()};
   				 my $stepsize = shift;
   				 my $windowsize = shift; 
				 my %final;
				
	foreach my $k (sort { $hash{$a} <=> $hash{$b} } keys %hash){   ##sort the input hash																					
	
					my @locsort;	   #define the sorted array for each scaffold
					foreach my $k2 (keys ($hash{$k})) {	     
						push (my @loc, $k2);  #build array of nucleotides
						@locsort = sort {$b<=>$a} @loc;   #sort the array of nucleotides from biggest to smallest.  Backwards because we want to know what the max value is below.
														}
	
				my $p = 1;  #p is my scaffold bookmark.  We will always start at 1 for each scaffold
				my $max = $lengthinfo{$k};  #end of the scaffold. 
	
				for ($p; $p< $max; $p= $p+$stepsize) {   #this moves our window
	
						my $count = 0;				#we're counting SNPs- restart counter at 0 for each window
						my $i = $p;					#start window at scaffold bookmark	
						my $winsize = $windowsize+$i;			#our window range is our start to start + windowsize
	
					for ($i; $i <= $winsize; $i++) {	#check each nucleotide, 1 by 1, in our range, for a hash key		
			$count++	if exists $hashcheck{$k}{$i};  #if we found a hash key, count it.  
					}
								
			my $end = $i;
			if ($end > $max){
				$end = $max;
			}  
			my $start = $i-$windowsize;
			  
			$final{$k}{$start}{"count"} = $count;  #assing range to count value	
			$final{$k}{$start}{"end"}	= $end;							
													}  #done looping through  windows	
	  }	#done looping through scaffolds	
				return (%final);				
				}
				


#######			
#######Sub filter window.  Improves upon above sliding window script.  
#######

#note: currently generating error 8-30-16 with MUMMer data

	sub Diversity_sliding_window {  
	
	 my %hash = %{shift()};
	 my %hashcheck = %{shift()};
	 my $stepsize = shift;
	 my $windowsize = shift; 
	 my %final;
				
	foreach my $k (sort { $hash{$a} <=> $hash{$b} } keys %hash){   ##sort the input hash																					
	
					my @locsort;	   #define the sorted array for each scaffold
					foreach my $k2 (keys ($hash{$k})) {	     
						push (my @loc, $k2);  #build array of nucleotides
						@locsort = sort {$b<=>$a} @loc;   #sort the array of nucleotides from biggest to smallest.  Backwards because we want to know what the max value is below.
														}
	
				my $p = 1;  #p is my scaffold bookmark.  We will always start at 1 for each scaffold
				my $max = $lengthinfo{$k};  #end of the scaffold. 
	
				for ($p; $p< $max; $p= $p+$stepsize) {   #this moves our window
	
						my $count = 0;	#we're counting SNPs- restart counter at 0 for each window
						my $primerCount = 0;
						my $windivscore = 0;   #set the window's diversity score to 0.
						my $i = $p;					#start window at scaffold bookmark	
						my $winsize = $windowsize+$i;			#our window range is our start to start + windowsize
						my $maskcount = 0;
						my $partscore = 0;
					for ($i; $i <= $winsize; $i++) {	#check each nucleotide, 1 by 1, in our range, for a hash key		

			$count++	if exists $hashcheck{$k}{$i};  #if we found a hash key in the hash we are counting things, count it.  
						if ($winsize - $i <= 20  || $winsize - $i > $windowsize - 20){#Checking here for exact conserved start/end for loci.
														#this will ensure that primers are likely to amplify.
											$primerCount++	if exists $hashcheck{$k}{$i}; 			
						}

						if (exists $reporthash{$k}{$i}) {

					my $toadd =  $reporthash{$k}{$i}{"Numnucs"};  #score is simply the number of unique nucleotides- 1-3 (A, T, G, C- reference is assumed as 4th nucleotdie).
					my $participants =  $reporthash{$k}{$i}{"NumSNPs"};  #Number of strains participating in each site
					$windivscore = $windivscore+$toadd;
					$partscore = $partscore+$participants;
						}				
											
					}
								
			my $end = $i;  #build keystring as name range
			my $start = $i-$windowsize;

			#Don't want windows to be longer than the scaffold- fix that here.
			if ($end > $max){
				$end = $max;
			}  
						
			$final{$k}{$start}{"primerCount"} = $primerCount;		
			$final{$k}{$start}{"participants"} = $partscore;
			$final{$k}{$start}{"count"} = $count;  #assign range to count value	
			$final{$k}{$start}{"end"}	= $end;		
			$final{$k}{$start}{"diversity"}	= $windivscore;			#assign diversity score
			$final{$k}{$start}{"mask"} = $maskcount;
													}  #done looping through  windows	
	  }	#done looping through scaffolds	
				return (%final);				
				}
				

	sub build_SNP_hash { 
	
				 my $file = shift;
				 my %outhash = %{shift()};
				my $linecount = 0;
	open (my $fh, "<", $file)
	  	  or die "couldn't open '$file' $!"; 		

	  	  if ($mummer eq 'TRUE' ) {
			 while (<$fh>) {
			 if ($linecount > 5) {#remove first 5 lines
				 chomp;
				 $linecount++;
				 }
				else {
	        	chomp ;
	        	next unless length;	
			    my @split = split(/\s+/); 
	        	next unless $split[11];	  	  	
				    my $queryloc= $split[11];
			  	    my $hitname= $file;
			  	    my $quernuc = $split[1];
			  	    my $hitbase = $split[3];
			  	    my $querbase = $split[2];
			  	 unless ($hitbase eq 'N' ||$hitbase eq '.' || $querbase eq 'N' || $querbase eq '.'){
		  	     	++$outhash{$queryloc}{$quernuc};  ##keep track of the number of SNPs at each basepair, for each scaffold
		  	    	}
		  	    	}
		  	}
			return (%outhash); 	
			}
				
	  	  else {
	   		 while (<$fh>) {
			    chomp ;
			    next unless length;
			    my @split = split(/\s+/);   
			    my $queryloc = $split[0];
			    my $hitname = $split[1];
			    my $quernuc = $split[2];
				++$outhash{$queryloc}{$quernuc};  ##keep track of the number of SNPs at each basepair, for each scaffold
			      } 						
	   		return (%outhash);
	   		 } 
 }


####OPEN SNP FILE AND BUILD HASH OF NUCLEOTIDE OCCURRENCES



 
sub build_SNP_NUC_hash {

		 my $file = shift;
		 my %outhash = %{shift()};
			open (my $fh, "<", $file)  or die "couldn't open '$file' $!"; 
			my $linecount=1;	
				#separate parsing for mummer or non-mummer
			if ($mummer eq 'TRUE' ) {
				 while (<$fh>) {
				 if ($linecount > 5) {#remove first 5 lines
				 chomp;
				 $linecount++;
				 }
				 else{		 
		        	chomp ;
		        	next unless length;	
				    my @split = split(/\s+/); 
		        	next unless $split[11];	  	  	
				    my $queryloc= $split[11];
			  	    my $hitname= $file;
			  	    my $quernuc = $split[1];
			  	    my $hitbase = $split[3];
			  	    my $querbase = $split[2];
			  	    unless ($hitbase eq 'N' ||$hitbase eq '.' || $querbase eq 'N' || $querbase eq '.'){
	  	    		 $outhash{$queryloc}{$quernuc} .="$hitbase";  ##keep track of the number of SNPs at each basepair, for each scaffold
	  	   			}
	  	   		}
	  	    	}
				return (%outhash);
  	 		}	
		 else {   					
			while (<$fh>) {
		    chomp ;
		    next unless length;
		    my @split = split(/\s+/);   
		    my $queryloc = $split[0];
		    my $hitname = $split[1];
		    my $quernuc = $split[2];
		    my $hitnuc = $split[3];
			my $hitbase = $split[5];
			$outhash{$queryloc}{$quernuc} .= "$hitbase";
				}
	  	return (%outhash);
		}
}
				

##
#Last updated 10-21-16
#Permutation on SNP tracker.
#Main distinction is to keep track of the SNPs while keeping file info.
##

sub build_Clade_Hash {
	#cladetracker{$name} = $clade;	
		 my $file = shift;
		 my %outhash = %{shift()};
			open (my $fh, "<", $file)  or die "couldn't open '$file' $!"; 
			my $linecount=1;	
				#separate parsing for mummer or non-mummer
			if ($mummer eq 'TRUE' ) {
				 while (<$fh>) {
				 if ($linecount > 5) {#remove first 5 lines
				 chomp;
				 $linecount++;
				 }
				 else{		 
		        	chomp ;
		        	next unless length;	
				    my @split = split(/\s+/); 
		        	next unless $split[11];	  	  	
				    my $queryloc= $split[11];
			  	    my $hitname= $file;
			  	    my $quernuc = $split[1];
			  	    my $hitbase = $split[3];
			  	    my $querbase = $split[2];
			  	    unless ($hitbase eq 'N' ||$hitbase eq '.' || $querbase eq 'N' || $querbase eq '.'){
	  	    		 $outhash{$queryloc}{$file}{$quernuc}++;  ##keep track of the number of SNPs at each basepair, for each scaffold
	  	   			}
	  	   		}
	  	    	}
				return (%outhash);
  	 		}	
		 else {   					
			while (<$fh>) {
		    chomp ;
		    next unless length;
		    my @split = split(/\s+/);   
		    my $queryloc = $split[0];
		    my $hitname = $split[1];
		    my $quernuc = $split[2];
		    my $hitnuc = $split[3];
			my $hitbase = $split[5];
			$outhash{$queryloc}{$file}{$quernuc}++;
				}
	  	return (%outhash);
		}
}
	

##
#Sliding windows script to check by clade. 
#Input is result of build_Clade_Hash
#
#Output is Scaffold -> Nuc number -> Taxon  Number of SNPs in window
#
#Last updated 10-21-16
##


sub Clade_sliding_window {  
			
	 my %hash = %{shift()};
	 my $stepsize = shift;
	 my $windowsize = shift; 
	 my %final;
	 my %finalAlt;
	  open (my $log, '>', "log.txt") or die "could not open log.txt $!";	

	foreach my $thisScaffold (sort { $hash{$a} <=> $hash{$b} } keys %hash){   ##sort the input hash by scaffold number																				
			my $p = 0;  #p is my nuc number in scaffold.  Always starts at 0
			my $max = $lengthinfo{$thisScaffold};  #end of the scaffold. 

			for ($p; $p<= $max; $p= $p+($stepsize)) {   #this moves our window
				foreach my $thisTaxon (keys $hash{$thisScaffold}){
					my $count = 0;	#we're counting SNPs- restart counter at 0 for each window
					my $i = $p;					#start window at scaffold bookmark	
					my $winsize = $windowsize+$i;			#our window range is our start to start + windowsize
				for ($i; $i <= $winsize; $i++) {	#check each nucleotide, 1 by 1, in our range, for a hash key	
		$count++	if exists $hash{$thisScaffold}{$thisTaxon}{$i};  #if we found a hash key in the hash we are counting things, count it.  					
				}
							
		my $end = $i;  #build keystring as name range
		my $start = $p;

		#Don't want windows to be longer than the scaffold- fix that here.
		if ($end > $max){
			$end = $max;
		}  	

		###Look up this window in refHash to see how many nucs aligned
					my $thisTaxonClip = removeRefandExtra($thisTaxon);
	my $nucsAlignedforWindow;
		if (exists $windowCoverageHash{$thisTaxonClip}{$thisScaffold}{$start}){ 
		$nucsAlignedforWindow = $windowCoverageHash{$thisTaxonClip}{$thisScaffold}{$start};
		}else{
			#print $log "$thisTaxonClip\t$thisScaffold\t$start\n";
			$nucsAlignedforWindow = 0;  #assume that there is no entry because there was no alignment.
		}

		my $adjustedSNP = $count + ($windowsize - $nucsAlignedforWindow);

		if ($adjustedSNP > $windowsize) {
			print $log "$thisTaxonClip\t$thisScaffold\t$start\t$count\t$nucsAlignedforWindow\n";
			#die "Error: more SNPs than window size.  Likely error in crosschecking hashes\n";
		}

		#print "Counted $count SNPs\n.  Now, Adjusting by adding $windowsize - $nucsAlignedforWindow\n";

		$final{$thisScaffold}{$thisTaxonClip}{$start}{"count"} = $count;  #assign range to count value	
		$final{$thisScaffold}{$thisTaxonClip}{$start}{"end"}	= $end;	

		$finalAlt{$thisScaffold}{$start}{$thisTaxonClip}{"count"}	= $count;	
		#$finalAlt{$thisScaffold}{$start}{$thisTaxon}{"end"}	= $end;	
		$finalAlt{$thisScaffold}{$start}{$thisTaxonClip}{"adjustedCount"}	= $adjustedSNP;	



				}#done with this taxon
			}  #done looping through  windows for this scaffold
			print "Sliding windows: finished with $thisScaffold\n";
	} #done with all scaffolds
return( \%final, \%finalAlt);

		}



#DetermineBestClade
#This subroutine goes through all windows, and determines the most likely source of that DNA.
#To do this, take the result of Clade_Sliding_Window and look at the minimum and average number of SNPs across all taxa for each clade
#

sub DetermineBestClade {  
	 my %windowHash = %{shift()};
	 my %outTracker;

		#$finalAlt{$thisScaffold}{$start}{$thisTaxon}{"count"}	= $count;	

	foreach my $thisScaffold (keys %windowHash){  
		foreach my $start (keys $windowHash{$thisScaffold}){
			foreach my $thisClade (keys %inverseCladeTracker) {
				my @taxonValues =();
				my $thisAverage =0;
				my $thisMin =0;
				foreach my $thisTaxon (keys $inverseCladeTracker{$thisClade}){
					$thisTaxon=removeRefandExtra($thisTaxon);
					my $thisValue = $windowHash{$thisScaffold}{$start}{$thisTaxon}{"adjustedCount"};
					
					push @taxonValues, $thisValue;
					}
				#now, calculate average and minimum
				if ( scalar @taxonValues > 1 ){
				$thisAverage = avg(@taxonValues);
				$thisMin =  min(@taxonValues);  #undefined if empty
				} elsif(scalar @taxonValues ==1 ) {
					$thisMin =  $taxonValues[0];
					$thisAverage = 0;
				 }
				else {$thisAverage = 0; $thisMin= 0;}

				$outTracker{$thisScaffold}{$start}{$thisClade}{'min'} = $thisMin;
				
				$outTracker{$thisScaffold}{$start}{$thisClade}{'avg'} = $thisAverage;
			}
		}
	}
	return(\%outTracker);
}

#calculate average of array.
#from Stack Overflow
#http://stackoverflow.com/questions/15696633/taking-average-of-array
sub avg { sum(@_)/@_ }


#sub clean_taxon_name
#input is a string
#output is the same string with:
#everything up to the last '/' removed (ie the path)
#everything after and including .fasta (ie the extension) 
#11-21-16



sub CleanTaxonName {

my $stringToClean = shift;
$stringToClean =~s/.*\///;#remove path prefix
$stringToClean =~s/.fasta.*//;#remove extensions

return $stringToClean;
}



#written by Bradford Condon, Dec 13, 2016
#Goal of subroutine:
#Read in report file on number of nucleotides aligned for each scaffold
####
####
#Input format
#Query: B71_SZL_masked	Subject: 87-120_masked
#Identifying repeat regions
#B71_SZL_masked.87-120_masked.BLAST
#chr1	0	0
#above are line types 1, 2, 3, and 4.  want line type 4.



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

	

	#format is scaffold, start, numberOfNucsAligned
	$outputArray{$quer}{$split[0]}{$split[1]} = $split[2];
}
return %outputArray;
}

##
#Created 12-20-16
#Another subroutine to read in reports on 1kb windows to determine how many nucs were aligned in the BLAST.
#Unlike hte above script, this reads in a file where each file is just a tab delimited file, 
#instead of above where the report had to be parsed and each pairwise comparison seperated.
##

sub ReadWindowAlignmentSeperateFile {
	my $filename  = shift;

open(my $fh, "<", $filename) or die "Could not open file '$filename' $!";

while (my $line = <$fh>) {
	chomp $line;
	if ($line =~ /versus/) { 
	next;
	}#check if line type 1
		my @split = split("\t", $line);

	#format is scaffold, start, numberOfNucsAligned
					
					my $filenameadjust = removeRefandExtra($filename);


	$windowCoverageHash{$filenameadjust}{$split[0]}{$split[1]} = $split[2];

}


}

sub removePathAndExtension{
my $string = shift;
$string =~ s{.*/}{};      # removes path  
$string =~ s{\.[^.]+$}{}; # removes extension
return $string;
}


sub removeRefandExtra{
my $string = shift;

$string =~s{.*_v_}{};
$string =~s{_.*}{};
return $string;
}
