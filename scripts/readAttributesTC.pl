#!/usr/bin/env perl
#use strict;
#use warnings;

(my $infile) = @ARGV;

open(USRFL, "<", $infile) or die "can't open $infile $!\n";
print "Chr,Strand,Start,End,mismatch,alignment,readID,count,copy,length,seq,end_base\n";
my @data = ();
while(<USRFL>) {

	chomp;
	chop($_) if ($_ =~ m/\r$/);
	my @line = split(/\t/, $_);
	if ($line[1] == 4){
		next;
	}
	my $readID = $line[0];
	if ($line[1] == 16){
		$strand = "-";
	}
	else{
		$strand = "+";
	}
	my $chrom = $line[2];
	my $start = $line[3];
	my $read = ($line[9]);
	my $conv = $line[12];
	if( $strand eq "-" ) {
                $read =~ tr/ACGT/TGCA/;
                $conv =~ tr/ACGT/TGCA/;
        }
	my $tag = $conv;
        my $conversionMismatches = 0;
        my $totalMismatches = 0;
        while( $tag =~ /\D*(\d*)([A,C,G,T])(.*)/ ){
                my $misbase = substr($read, $1, 1);
                if ($2 eq "T" and $misbase eq "C"){
                        $conversionMismatches++;
                }
                else{
			$totalMismatches++;
                }
		$tag = $3;
        }
	my $end = $start + length($read);
	if (scalar(@data) > 0){
		if ($readID ne $data[6]){
			while (scalar(@data) > 0){
				print join(',',splice(@data, 0, 12))."\n";
			}
		}
	}

if (length($read) > 19) {

	push(@data, ($chrom,$strand,$start,$end,));

	if ($conversionMismatches + $totalMismatches <= 2){
		if ($conversionMismatches == 0 and $totalMismatches == 0){
			push(@data, "None");
		}
		elsif ($conversionMismatches == 0){
			push(@data, "Other_$totalMismatches");
		}
		elsif ($totalMismatches == 0){
			push(@data, "T2C_$conversionMismatches");
		}
		else{
			push(@data, "T2CandOther_" . ($conversionMismatches + $totalMismatches));
		}
	}
	
	if (scalar(@data) > 12) {
		push (@data, "Multimapper");
		$data[5] = "Multimapper";
	}
	else  {
		push (@data, "Single");
	}

	push (@data, $readID);

	if ($readID =~ /(\d+)-/) {
		push (@data, $');
		if ($' > 1) {
			push (@data, "Redundant");
		}
		else {
			push (@data, "Unique");
		}
	}
	push (@data, length($read));

	if ($strand == "+") {
		my $revread = reverse($read);
		push (@data, $revread);
	}	
	else {
		push (@data, $read);
	}

	if (substr($read, -1, 1) eq "G"){
		push(@data, "G_End");
	}
	else{
		push(@data, "Non_G_End");
	}
		
}
}
while (scalar(@data) > 0){
	print join(',',splice(@data, 0, 12))."\n";
}
close USRFL;
