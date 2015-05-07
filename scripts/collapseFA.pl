#!/usr/bin/env perl
# collapseFA.pl
sub descSort;
sub version;

my $option = $ARGV[0];
version if $option eq "-v" or $option eq "--version";

while(<STDIN>) {
	chomp;
	chop($_) if ($_ =~ m/\r$/);
	chomp;
	if( not /^>/ and /\w/ ) { $count{$_}++; }
}

$counter = 0;
foreach my $seq ( sort descSort(keys(%count)) ) {
	$counter++;
	print ">$counter-$count{$seq}\n$seq\n";
}


sub descSort {
	$count{$b} <=> $count{$a};
}


sub version {
print "version 2.0\n"; exit;
}
