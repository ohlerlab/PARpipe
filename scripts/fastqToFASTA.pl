#!/usr/bin/env perl
# fastqToFASTA.pl
sub version;

my $option = $ARGV[0];
version if $option eq "-v" or $option eq "--version";

my $lines = 0;

while(<STDIN>) {
	chomp;
	chop($_) if ($_ =~ m/\r$/);
	$lines++;
	if( /^\@/ ) { $ac = substr($_, 1); $x = 1;}
	elsif( /\w/ && $x == 1 ) {
		if( not /N/ ) { print ">$ac\n$_\n"; }
		else { $nCount++; }
		$x = 0;
	}
}
#close USRFL;
my $reads = $lines/4;
my $frac = 100*$nCount/$reads;
print STDERR "Total Reads: $reads\nReads filtered containing N: $nCount\nPercent Reads filtered containing N: $frac\n";


sub version {
print "version 2.0\n"; exit;
}


