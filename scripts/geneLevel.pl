#!/usr/bin/env perl
use List::Util qw(sum);
use Statistics::Basic qw(median);

my $file = $ARGV[0];
my %geneCount;
my %geneReads;
my %geneT2C;
my %geneCon;
my %geneUniq;
my %geneAlign;
my %geneType;
my %annot;
$annot{"5'utr"} = 0;
$annot{"intron"} = 1;
$annot{"coding"} = 2;
$annot{"3'utr"} = 3;
$annot{"start_codon"} = 4;
$annot{"stop_codon"} = 5;
my $regex = qr/.*(utr|coding|intron)-.*(utr|coding|intron)/;
$annot{$regex} = 6;
open(GTF, "pgrep '^\\w*\\t\\w*\\tgene\\t' $ARGV[1] |");
while(<GTF>){
 chomp;
 $_ =~ /gene_name "([^"]*)"/;
 my $genna = $1;
 $_ =~ /gene_type "([^"]*)"/;
 $geneType{$genna} = $1;
}
close(GTF);
open(IN, $file);
$firstline = <IN>;
my @order = (0, 0, 0, 0, 0, 0);
my @secs = split(",", $firstline);
for (my $i = 0; $i <= $#secs; ++$i){
 if ($secs[$i] =~ /GeneName/){
  $order[0] = $i;
 }
 elsif ($secs[$i] =~ /ReadCount/){
  $order[1] = $i;
 }
 elsif ($secs[$i] =~ /T2Cfraction/){
  $order[2] = $i;
 }
 elsif ($secs[$i] =~ /ConversionSpecificity/){
  $order[3] = $i;
 }
 elsif ($secs[$i] =~ /UniqueReads/){
  $order[4] = $i;
 }
 elsif ($secs[$i] =~ /Aligned to/){
  $order[5] = $i;
 }
}
while(<IN>){
 chomp;
 my @secs = split ",";
 if ($secs[$order[0]] eq ""){
  next;
 }
 $geneCount{$secs[$order[0]]}++;
 push(@{$geneReads{$secs[$order[0]]}}, $secs[$order[1]]);
 push(@{$geneT2C{$secs[$order[0]]}}, $secs[$order[2]]);
 if ($secs[$order[3]] ne "NA"){
  push(@{$geneCon{$secs[$order[0]]}}, $secs[$order[3]]);
 }
 push(@{$geneUniq{$secs[$order[0]]}}, $secs[$order[4]]);
 if (scalar(@{$geneAlign{$secs[$order[0]]}}) == 0){
  push(@{$geneAlign{$secs[$order[0]]}}, (0, 0, 0, 0, 0, 0, 0));
 }
 if (exists($annot{$secs[$order[5]]})){
  $geneAlign{$secs[$order[0]]}[$annot{$secs[$order[5]]}]++;
 }
 if (!exists($geneType{$secs[$order[0]]})){
  $geneType{$secs[$order[0]]} = $secs[$order[5]];
 }
}
close(IN);
print("GeneName,Sites,ReadCountSum,ReadCountMed,T2CFractionSum,T2CFractionMed,ConversionSpecificitySum,ConversionSpecificityMed,UniqueReadsSum,UniqueReadsMed,5'utr,Intron,Exon,3'utr,Start_codon,Stop_codon,Junction,GeneType\n");
foreach $key(sort(keys %geneReads)){
 print("$key," . $geneCount{$key} . "," . sum(@{$geneReads{$key}}) . "," . int(median(@{$geneReads{$key}})) . "," . sum(@{$geneT2C{$key}}) . "," . median(@{$geneT2C{$key}}) . "," . sum(@{$geneCon{$key}}) . "," . median(@{$geneCon{$key}}) . "," . sum(@{$geneUniq{$key}}) . "," . median(@{$geneUniq{$key}}) . "," . $geneAlign{$key}[0] . "," . $geneAlign{$key}[1] . "," . $geneAlign{$key}[2] . "," . $geneAlign{$key}[3] . "," . $geneAlign{$key}[4] . "," . $geneAlign{$key}[5] . "," . $geneAlign{$key}[6] . "," . $geneType{$key} . "\n");
}
