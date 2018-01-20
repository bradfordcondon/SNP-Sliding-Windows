

my $inputFile = "CML3070_masked_v_PH1_masked_out.txt";


my $outputFile = "$inputFile.swapped";

open(my $fh, "<", $inputFile)
or die "couldnt open '$inputFile' $!";	

open(my $Out, "", $outputFile)	
	while (<$fh>){
	chomp;
	my @split = split(/\t/);
	 my $c0 = $split[0];
	my $c1 = $split[1];
	my $c2 = $split[2];
	my $c3 = $split[3];
	my $c4 = $split[4];
	my $c5 = $split[5];
	my $c6 = $split[6];
	my $comment = $split[7];

	if ($comment) {
		print "$comment\n";
		}
	}
close $fh;

