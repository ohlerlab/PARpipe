#!/usr/bin/env perl
( $iniFile, $inputFile, $bitfile ) = @ARGV;

open(USRFL, "<", $iniFile) or die "Can't open $iniFile\n";
while(<USRFL>) {
	if( (not /^"BOWTIE_FILE"/ && not /^"OUTPUT"/) && length($_) > 0 ) {
		print $_;
	}

}
close USRFL;


if ($inputFile =~ /\.aligned.sam/) {$library = $`};

print "GENOME_2BIT_FILE=$bitfile";
print "\n";
print "SAM_FILE=$inputFile\=COLLAPSED\n\n";
print "OUTPUT_DISTRIBUTIONS_FILE="."$library\.distribution\n";
print "OUTPUT_GROUPS_FILE="."$library\.groups\n";
print "OUTPUT_CLUSTERS_FILE="."$library\.clusters\n";
print "OUTPUT_READS_FILE="."$library\_PARalyzer_Utilized.sam\n\n";

