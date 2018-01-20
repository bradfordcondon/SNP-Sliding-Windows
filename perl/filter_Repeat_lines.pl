

my $directory = "/Users/chet/uky/Fusarium_1-2-17/inputData/";



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
	my @split = split(/\t/);



}
   
   closedir(DIR);
