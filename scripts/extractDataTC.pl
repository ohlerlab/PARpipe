#!/usr/bin/perl

sub log10 {  my $n = shift;
        return log($n)/log(10);
}

my $fn = $ARGV[0];

print("Summary Stats," . $fn . "\nReads,");
$/ = undef;
open(IN, $fn . ".fastq.processing");
my $t = <IN>;
close(IN);
$/ = "\n";
$t =~ /Total Reads: (\d*)/;
my $reads = $1;
print($reads . "(100%)\nFiltered Total,");
$t =~ /\nReads filtered containing N: (\d*)/;
my $n = $1 + 0;
$t =~ /Too short reads:\s*(\d*)/;
my $short = $1;
my $filtered = $n + $short;
my $pfiltered = ($filtered / $reads) * 100;
printf($filtered . "(%.2f%)\nQuality,", $pfiltered);
my $pn = ($n / $reads) * 100;
printf($n . "(%.2f%)\nLength <20 nts,", $pn);
my $pshort = ($short / $reads) * 100;
printf($short . "(%.2f%)\n3' Adapter,", $pshort);
$t =~ /-a (\w*) /;
my $tadapter = $1;
$t =~ /$tadapter', length \d*, was trimmed (\d*) times.\n/;
$tadapter = $1;
my $ptadapter = ($tadapter / $reads) * 100;
printf($tadapter . "(%.2f%)\n5' Adapter,", $ptadapter);
$t =~ /-g (\w*) /;
my $fadapter = $1;
$t =~ /$fadapter', length \d*, was trimmed (\d*) times.\n/;
$fadapter = $1;
my $pfadapter = ($fadapter / $reads) * 100;
printf($fadapter . "(%.2f%)\n19 nt Marker,", $pfadapter);
$t =~ /, length 19, was trimmed (\d*) times.\n/;
my $nadapter = $1;
my $pnadapter = ($nadapter / $reads) * 100;
printf($nadapter . "(%.2f%)\n24 nt Marker,", $pnadapter);
$t =~ /, length 24, was trimmed (\d*) times.\n/;
my $twadapter = $1;
my $ptwadapter = ($twadapter / $reads) * 100;
printf($twadapter . "(%.2f%)\nPost-Processing,", $ptwadapter);
my $processed = $reads - $n - $short;
my $pprocessed = ($processed / $reads) * 100;
printf($processed . "(%.2f%)\nUniquely Aligned Reads,", $pprocessed);

my %readSummary;

open(IN, $fn . ".readcsv");
$t = <IN>;
@head = split(',', $t);
my @order = (0, 0, 0, 0);
for (my $i = 0; $i <= $#head; ++$i){
 if ($head[$i] =~ /(Aligned to)|(TranscriptLocation)/){
  $order[0] = $i;
 }
 elsif ($head[$i] =~ /mismatch/){
  $order[1] = $i;
 }
 elsif ($head[$i] =~ /end_base/){
  $order[2] = $i;
 }
 elsif ($head[$i] =~ /^count$/){
  $order[3] = $i;
 }
}
while (<IN>){
 chomp $_;
 @read = split(',', $_);
 $readSummary{$read[$order[0]]}{'reads'} += $read[$order[3]];
 $readSummary{$read[$order[0]]}{$read[$order[1]]} += $read[$order[3]];
 $readSummary{$read[$order[0]]}{$read[$order[2]]} += $read[$order[3]];
 if ($read[$order[3]] > 1){
  $readSummary{$read[$order[0]]}{'redundant'} += $read[$order[3]];
  $readSummary{$read[$order[0]]}{'sredundant'}++;
 }
 else{
  $readSummary{$read[$order[0]]}{'unique'}++;
 }
}
close(IN);
my @locations = sort(keys(%readSummary));
my $total = 0;
foreach my $location(@locations){
 if (exists($readSummary{$location}{'reads'})){
  $total += $readSummary{$location}{'reads'};
 }
}
$ptotal = ($total / $processed) * 100;
printf($total . "(%.2f%)\nClusters,", $ptotal);

my %clusterSummary;

open(IN, $fn . ".clusters.csv");
$t = <IN>;
@head = split(',', $t);
@order = (0, 0, 0, 0, 0, 0);
for (my $i = 0; $i <= $#head; ++$i){
 if ($head[$i] =~ /(Aligned to)|(TranscriptLocation)/){
  $order[0] = $i;
 }
 elsif ($head[$i] =~ /T2Cfraction/){
  $order[1] = $i;
 }
 elsif ($head[$i] =~ /ConversionSpecificity/){
  $order[2] = $i;
 }
 elsif ($head[$i] =~ /endG_fraction/){
  $order[3] = $i;
 }
 elsif ($head[$i] =~ /RedundantSeqFraction/){
  $order[4] = $i;
 }
 elsif ($head[$i] =~ /RedundantCopyFraction/){
  $order[5] = $i;
 }
}
while (<IN>){
 chomp $_;
 @cluster = split(',', $_);
 $clusterSummary{$cluster[$order[0]]}{'clusters'}++;
 $clusterSummary{$cluster[$order[0]]}{'T2CF'} += $cluster[$order[1]];
 $clusterSummary{$cluster[$order[0]]}{'ConSpec'} += $cluster[$order[2]];
 $clusterSummary{$cluster[$order[0]]}{'GFrac'} += $cluster[$order[3]];
 $clusterSummary{$cluster[$order[0]]}{'RSeqF'} += $cluster[$order[4]];
 $clusterSummary{$cluster[$order[0]]}{'RCopyF'} += $cluster[$order[5]];
}
close(IN);
$total = 0;
foreach my $location(@locations){
 if (exists($clusterSummary{$location}{'clusters'})){
  $total += $clusterSummary{$location}{'clusters'};
 }
}
print($total . "\nGroups,");

my %groupSummary;

open(IN, $fn . ".groups.csv");
$t = <IN>;
@head = split(',', $t);
@order = (0, 0, 0, 0, 0, 0);
for (my $i = 0; $i <= $#head; ++$i){
 if ($head[$i] =~ /(Aligned to)|(TranscriptLocation)/){
  $order[0] = $i;
 }
 elsif ($head[$i] =~ /T2Cfraction/){
  $order[1] = $i;
 }
 elsif ($head[$i] =~ /ConversionSpecificity/){
  $order[2] = $i;
 }
 elsif ($head[$i] =~ /endG_fraction/){
  $order[3] = $i;
 }
 elsif ($head[$i] =~ /RedundantSeqFraction/){
  $order[4] = $i;
 }
 elsif ($head[$i] =~ /RedundantCopyFraction/){
  $order[5] = $i;
 }
}
while (<IN>){
 chomp $_;
 @group = split(',', $_);
 $groupSummary{$group[$order[0]]}{'groups'}++;
 $groupSummary{$group[$order[0]]}{'T2CF'} += $group[$order[1]];
 $groupSummary{$group[$order[0]]}{'ConSpec'} += $group[$order[2]];
 $groupSummary{$group[$order[0]]}{'GFrac'} += $group[$order[3]];
 $groupSummary{$group[$order[0]]}{'RSeqF'} += $group[$order[4]];
 $groupSummary{$group[$order[0]]}{'RCopyF'} += $group[$order[5]];
}
close(IN);
$total = 0;
foreach my $location(@locations){
 if (exists($groupSummary{$location}{'groups'})){
  $total += $groupSummary{$location}{'groups'};
 }
}
print($total . "\n\n");

print("Reads Level,Total,T2C Fraction,Conv Spec,G Fraction,Red Seq Fraction,Red Copy Fraction\n");
my @mismatches = ("Other_1", "Other_2", "T2C_1", "T2C_2", "T2CandOther_2");
foreach my $location(@locations){
 my @conversions;
 print($location);
 if (!exists($readSummary{$location}{'reads'})){
  print("," . 0 . ",-,-,-,-,-\n");
  next;
 }
 print(",$readSummary{$location}{'reads'}");
 foreach my $mismatch(@mismatches){
  if (exists($readSummary{$location}{$mismatch})){
   push(@conversions, $readSummary{$location}{$mismatch});
  }
  else{
   push(@conversions, 0);
  }
 }
 my $T2C_fraction = (($conversions[2] + $conversions[3])/$readSummary{$location}{'reads'});
 printf(",%.2f", $T2C_fraction);
 if ($conversions[2] + $conversions[3] != 0){
  my $convSpec = log10(($conversions[2] + $conversions[3])/($conversions[0] + $conversions[1] + $conversions[4] + 1));
  printf(",%.2f", $convSpec);
 }
 else{
  print(",0.00");
 }
 my $frac;
 if (exists($readSummary{$location}{'G_End'}) and exists($readSummary{$location}{'Non_G_End'})){
  $frac = $readSummary{$location}{'G_End'}/$readSummary{$location}{'reads'};
 }
 elsif (exists($readSummary{$location}{'G_End'})){
  $frac = 1;
 }
 else{
  $frac = 0;
 }
 printf(",%.2f", $frac);
 if (exists($readSummary{$location}{'sredundant'}) and exists($readSummary{$location}{'unique'})){
  $frac = $readSummary{$location}{'sredundant'}/($readSummary{$location}{'sredundant'} + $readSummary{$location}{'unique'});
 }
 elsif (exists($readSummary{$location}{'sredundant'})){
  $frac = 1;
 }
 else{
  $frac = 0;
 }
 printf(",%.2f", $frac);
 if (exists($readSummary{$location}{'redundant'})){
  $frac = $readSummary{$location}{'redundant'}/$readSummary{$location}{'reads'};
 }
 else{
  $frac = 0;
 }
 printf(",%.2f\n", $frac);
}

print("\n");

print("Cluster Level,Total,T2C Fraction,Conv Spec,G Fraction,Red Seq Fraction,Red Copy Fraction\n");
foreach my $location(@locations){
 print($location);
 if (!exists($clusterSummary{$location}{'clusters'})){
  print("," . 0 . ",-,-,-,-,-\n");
  next;
 }
 print(",$clusterSummary{$location}{'clusters'}");
 my $frac = $clusterSummary{$location}{'T2CF'}/$clusterSummary{$location}{'clusters'};
 printf(",%.2f", $frac);
 my $frac = $clusterSummary{$location}{'ConSpec'}/$clusterSummary{$location}{'clusters'};
 printf(",%.2f", $frac);
 my $frac = $clusterSummary{$location}{'GFrac'}/$clusterSummary{$location}{'clusters'};
 printf(",%.2f", $frac);
 my $frac = $clusterSummary{$location}{'RSeqF'}/$clusterSummary{$location}{'clusters'};
 printf(",%.2f", $frac);
 my $frac = $clusterSummary{$location}{'RCopyF'}/$clusterSummary{$location}{'clusters'};
 printf(",%.2f\n", $frac);
}

print("\n");

print("Group Level,Total,T2C Fraction,Conv Spec,G Fraction,Red Seq Fraction,Red Copy Fraction\n");
foreach my $location(@locations){
 print($location);
 if (!exists($groupSummary{$location}{'groups'})){
  print("," . 0 . ",-,-,-,-,-\n");
  next;
 }
 print(",$groupSummary{$location}{'groups'}");
 my $frac = $groupSummary{$location}{'T2CF'}/$groupSummary{$location}{'groups'};
 printf(",%.2f", $frac);
 my $frac = $groupSummary{$location}{'ConSpec'}/$groupSummary{$location}{'groups'};
 printf(",%.2f", $frac);
 my $frac = $groupSummary{$location}{'GFrac'}/$groupSummary{$location}{'groups'};
 printf(",%.2f", $frac);
 my $frac = $groupSummary{$location}{'RSeqF'}/$groupSummary{$location}{'groups'};
 printf(",%.2f", $frac);
 my $frac = $groupSummary{$location}{'RCopyF'}/$groupSummary{$location}{'groups'};
 printf(",%.2f\n", $frac);
}
