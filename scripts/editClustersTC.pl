#!/usr/bin/env perl
#use strict;
#use warnings;

sub log10 {  my $n = shift;
	return log($n)/log(10);
}

my %geneSummary;
my $which = 1;
my $cutoff = $ARGV[0];

while(<STDIN>) {
	chomp;
	chop($_) if ($_ =~ m/\r$/);
	my @line = split(/\,/, $_);
#	if($line[13] =~ /[35]\'utr|coding|intron/i) {
#	my($gene, $reads, $conversions, $location, $motif) = ($line[15], $line[6], $line[10], $line[13], $line[16]);
#    my (Chr,Strand,Start,End,alignedTo,mismatch,alignment,ID,count,copy,readlength,seq,clusterID) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[12]);
    if (scalar(@line) == 18){
        ($location,$readCount,$endBase,$Chr,$Start,$End,$alignedTo,$name,$gene,$seq,$reads,$ModeLocation,$ModeScore,$conLocCount,$conEvCount,$nonConEvCount,$annotation,$Strand) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14],$line[15],$line[16],$line[17]);
    }
    else{
	($location,$readCount,$endBase,$Chr,$Start,$End,$alignedTo,$name,$gene,$seq,$reads,$conLocCount,$conEvCount,$annotation,$Strand) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10],$line[11],$line[12],$line[13],$line[14]);
	$which = 0;
    }

	$geneSummary{$gene}{$location}{'reads'} += $readCount;
	if ($readCount > 1){
		$geneSummary{$gene}{'fredundant'}{'red'}++;
		$geneSummary{$gene}{'credundant'}{'red'} += $readCount;
		$geneSummary{$gene}{'credundant'}{'all'} += $readCount;
	}
	else{
		$geneSummary{$gene}{'fredundant'}{'uni'}++;
		$geneSummary{$gene}{'credundant'}{'all'} += $readCount;
	}
	$geneSummary{$gene}{$endBase}{'reads'} += $readCount;
        if (scalar(@line) == 18){
		$geneSummary{$gene}{'summary'} = $Chr . "," . $Strand . "," . $Start . "," . $End . "," . $alignedTo . "," . $name . "," . $annotation . "," . $gene . "," . $seq . "," . $reads . "," . $ModeLocation . "," . $ModeScore . "," . $conLocCount . "," . $conEvCount . "," . $nonConEvCount;
	}
	else{
		$geneSummary{$gene}{'summary'} = $Chr . "," . $Strand . "," . $Start . "," . $End . "," . $alignedTo . "," . $name . "," . $annotation . "," . $gene . "," . $seq . "," . $reads . "," . $conLocCount . "," . $conEvCount;
	}
	$geneSummary{$gene}{'Link'} = $Chr . ":" . $Start . "-" . $End . "\n";
#	$geneSummary{$gene}{'conversions'} += $conversions;
#	$geneSummary{$gene}{$location}++;
#	$geneSummary{$gene}{$motif}++;	
}

if ($which){
	print "Chr,Strand,Start,End,Aligned to,GeneName,AnnotationSource,ClusterID,ClusterSequence,ReadCount,ModeLocation,ModeScore,ConversionLocationCount,ConversionEventCount,NonConversionEventCount,T2Cfraction,ConversionSpecificity,endG_fraction,RedundantSeqFraction,RedundantCopyFraction,UniqueReads,None,Other_1,T2C_1,Link\n";
}
else{
	print "Chr,Strand,Start,End,Aligned to,GeneName,AnnotationSource,GroupID,GroupSequence,ReadCount,ConversionLocationCount,ConversionEventCount,T2Cfraction,ConversionSpecificity,endG_fraction,RedundantSeqFraction,RedundantCopyFraction,UniqueReads,None,Other_1,T2C_1,Link\n";
}
#my @motifs = ("TATTTTTAT", "TATTTTAT", "TATTTAT", "TATTAT", "NA");
my @locations = ("None", "Other_1", "Other_2", "T2C_1", "T2C_2", "T2CandOther_2");
my @storedData = ("reads");
my %clusters;
foreach my $gene (keys(%geneSummary)) {
	my @conversions = ();
	foreach my $dataType ( @storedData ) {
                foreach my $location( @locations ) {
                        if( exists($geneSummary{$gene}{$location}) ) {
                                push(@conversions, $geneSummary{$gene}{$location}{$dataType});
                        }
                        else {
                              	push(@conversions, 0);
                        }
                }
        }
	if (($conversions[3] + $conversions[4] == 0 or log10(($conversions[3] + $conversions[4])/($conversions[1] + $conversions[2] + $conversions[5] + 1)) <= $cutoff) and $which){
		next;
        }
	$clusters{$gene} = $geneSummary{$gene}{'summary'};
	my $T2C_fraction = (($conversions[3] + $conversions[4])/($conversions[0] + $conversions[1] + $conversions[2] + $conversions[3] + $conversions[4] + $conversions[5]));
	$clusters{$gene} .= ",$T2C_fraction";
	if ($conversions[3] + $conversions[4] > 0){
		$convSpec = log10(($conversions[3] + $conversions[4])/($conversions[1] + $conversions[2] + $conversions[5] + 1));
	}
	else{
		$convSpec = "NA";
	}
	$clusters{$gene} .= ",$convSpec";
	if( exists($geneSummary{$gene}{'G_End'}{'reads'}) and exists($geneSummary{$gene}{'Non_G_End'}{'reads'})){
                my $frac = $geneSummary{$gene}{'G_End'}{'reads'}/($geneSummary{$gene}{'G_End'}{'reads'}+$geneSummary{$gene}{'Non_G_End'}{'reads'});
                $clusters{$gene} .= ",$frac";
        }
        elsif( exists($geneSummary{$gene}{'G_End'}{'reads'})){
                $clusters{$gene} .= ",1";
        }
        else{
                $clusters{$gene} .= ",0";
        }
	if( exists($geneSummary{$gene}{'fredundant'}{'uni'}) and exists($geneSummary{$gene}{'fredundant'}{'red'}) ){
		my $frac = $geneSummary{$gene}{'fredundant'}{'red'}/($geneSummary{$gene}{'fredundant'}{'uni'} + $geneSummary{$gene}{'fredundant'}{'red'});
		$clusters{$gene} .= ",$frac";
	}
	elsif(exists($geneSummary{$gene}{'fredundant'}{'red'})){
		$clusters{$gene} .= ",1";
	}
	else{
		$clusters{$gene} .= ",0";
	}
	if( exists($geneSummary{$gene}{'credundant'}{'red'}) ){
                my $frac = $geneSummary{$gene}{'credundant'}{'red'}/$geneSummary{$gene}{'credundant'}{'all'};
		$clusters{$gene} .= ",$frac";
        }
	else{
		$clusters{$gene} .= ",0";
	}
        if( exists($geneSummary{$gene}{'fredundant'}{'uni'}) and exists($geneSummary{$gene}{'fredundant'}{'red'}) ){
                my $frac = $geneSummary{$gene}{'fredundant'}{'uni'} + $geneSummary{$gene}{'fredundant'}{'red'};
                $clusters{$gene} .= ",$frac";
        }
	elsif (exists($geneSummary{$gene}{'fredundant'}{'uni'})){
                my $frac = $geneSummary{$gene}{'fredundant'}{'uni'};
                $clusters{$gene} .= ",$frac";
        }
	elsif (exists($geneSummary{$gene}{'fredundant'}{'red'})){
                my $frac = $geneSummary{$gene}{'fredundant'}{'red'};
                $clusters{$gene} .= ",$frac";
        }
	$clusters{$gene} .= "," . $conversions[0] . "," . $conversions[1] . "," . $conversions[3] . "," . $geneSummary{$gene}{'Link'};
}

foreach $gene(sort(keys(%clusters))){
	print($clusters{$gene});
}
