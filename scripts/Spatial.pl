#!/usr/bin/perl
use List::Util qw[min max];
sub version;

if ($#ARGV == -1){
 print("\nProgram:\tSpatial.pl (v1.4)\nAuthor:\tNicholas Jacobs (nickcjacobs\@gmail.com\nSummary:\tAnalyzes PAR-CLIP output from PARalyzer and creates spatial distributions of the data.\n\nUsage:\tperl Spatial.pl [OPTIONS]\n\nOptions:\n\t-g\tgtf file to be used as reference.\n\n\t-a\tCluster file to be analyzed.\n\n\t-b\tOptional second cluster file.\n\n\t-strict\tOnly annotates to transcripts in which we have confidence.\n\n\t-t\tFile containing transcript-level FPKM data from Cufflinks. Optional.\n\n\t-res\tResolution of produced pdf files. Default = 1 (highest resolution), higher numbers for lower resolution. 5 as recommended maximum.\n\n\t-path\tPath for R scripts and packages. Default = ./\n");
 die "\n";
}
my $gtf;
my $fn;
my $fn2;
my $strict = 0;
my $rseq = "";
my %trab;
for ($i = 0; $i <= $#ARGV; ++$i){
 if ($ARGV[$i] eq "-g"){
  $gtf = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-a"){
  $fn = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-b"){
  $fn2 = $ARGV[$i + 1];
  $yn = "y";
 }
 elsif ($ARGV[$i] eq "-strict"){
  $strict = 1;
 }
 elsif ($ARGV[$i] eq "-t"){
  $rseq = $ARGV[$i + 1];
  open(RS, $rseq);
  while(<RS>){
   chomp;
   my @secs = split;
   if ($secs[0] ne "tracking_id"){
    $trab{$secs[0]} = $secs[9];
   }
  }
 }
}
#my @rand = ();
my @Aligned = ();
my @chr = ();
my @reg = ();
my @start = ();
my @end = ();
my @strand = ();
my @exnu = ();
my @genna = ();
my @trid = ();
my @accept = ();
my $genname = "";
my @Strand = ();
my @Chr = ();
my @Start = ();
my @End = ();
my $beg = 1;
if ($yn eq "y"){
 %genco31;
 %genco51;
 %gencoc1;
 %gencoii1;
 %gencole1;
 %gencoli1;
 %genco32;
 %genco52;
 %gencoc2;
 %gencoii2;
 %gencole2;
 %gencoli2;
 @temp31 = ();
 @temp51 = ();
 @tempc1 = ();
 @tempi1 = ();
 @temple1 = ();
 @templi1 = ();
 @temp32 = ();
 @temp52 = ();
 @tempc2 = ();
 @tempi2 = ();
 @temple2 = ();
 @templi2 = ();
 @genli32 = ();
 @genli52 = ();
 @genlic2 = ();
 @genliii2 = ();
 @genlile2 = ();
 @genlili2 = ();
 @genli31 = ();
 @genli51 = ();
 @genlic1 = ();
 @genliii1 = ();
 @genlile1 = ();
 @genlili2 = ();
 @genl31 = ();
 @genl32 = ();
 @genl51 = ();
 @genl52 = ();
 @genlc1 = ();
 @genlc2 = ();
 @genlii1 = ();
 @genlii2 = ();
 @genlle1 = ();
 @genlle2 = ();
 @genlli1 = ();
 @genlli2 = ();
 %genp31;
 %genp32;
 %genp51;
 %genp52;
 %genpc1;
 %genpc2;
 %genpii1;
 %genpii2;
 %genple1;
 %genple2;
 %genpli1;
 %genpli2;
}
else{
 %genco3;
 %genco5;
 %gencoc;
 %gencoii;
 %gencole;
 %gencoli;
 @temp3 = ();
 @temp5 = ();
 @tempc = ();
 @tempi = ();
 @temple = ();
 @templi = ();
 @genli3 = ();
 @genli5 = ();
 @genlic = ();
 @genliii = ();
 @genlile = ();
 @genlili = ();
 @temp3e = ();
 @temp5e = ();
 @tempce = ();
 @tempid = ();
 @templid = ();
 @templed = ();
 @genp3 = ();
 @genp5 = ();
 @genpc = ();
 @genpii = ();
 @genple = ();
 @genpli = ();
}
my @genma3 = ();
my @genma5 = ();
my @genmac = ();
my @genmaii = ();
my @genmale = ();
my @genmali = ();
my @genes3 = ();
my @genes5 = ();
my @genesc = ();
my @genesii = ();
my @genesle = ();
my @genesli = ();
my %stranded;

sub version{
 print "version 1.4\n"; exit;
}

sub addOn{
 if ($_[3]){
  return $_[2] - $_[0];
 }
 else{
  return $_[1] - $_[2];
 }
}

sub manage{
 if ($whic == 0){
  foreach $name (@genesii){
   if (scalar(@{$gencoii{$name}}) > 1){
    my @ind = sort { ${$gencoii{$name}}[$a] <=> ${$gencoii{$name}}[$b] } 0..$#{$gencoii{$name}};
    @{$gencoii{$name}} = @{$gencoii{$name}}[@ind];
    push(@genliii, $gencoii{$name}[1] - $gencoii{$name}[0]);
    if (scalar(@{$gencoii{$name}}) > 2){
     for (my $j = 1; $j < scalar(@{$gencoii{$name}})-1; ++$j){
      if ($gencoii{$name}[$j] - $gencoii{$name}[$j - 1] < $gencoii{$name}[$j + 1] - $gencoii{$name}[$j]){
       $w = 1;
      }
      else{
       $w = 2;
      }
      if ($w == 1){
       push(@genliii, $gencoii{$name}[$j] - $gencoii{$name}[$j - 1]);
      }
      else{
       push(@genliii, $gencoii{$name}[$j + 1] - $gencoii{$name}[$j]);
      }
     }
    }
    push(@genliii, ${$gencoii{$name}}[$#{$gencoii{$name}}] - ${$gencoii{$name}}[$#{$gencoii{$name}} - 1]);
   }
  }
  foreach $name (@genes3){
   if (scalar(@{$genco3{$name}}) > 1){
    my @ind = sort { ${$genco3{$name}}[$a] <=> ${$genco3{$name}}[$b] } 0..$#{$genco3{$name}};
    @{$genco3{$name}} = @{$genco3{$name}}[@ind];
    push(@genli3, $genco3{$name}[1] - $genco3{$name}[0]);
    if (scalar(@{$genco3{$name}}) > 2){
     for (my $j = 1; $j < scalar(@{$genco3{$name}})-1; ++$j){
      if ($genco3{$name}[$j] - $genco3{$name}[$j - 1] < $genco3{$name}[$j + 1] - $genco3{$name}[$j]){
       $w = 1;
      }
      else{
       $w = 2;
      }
      if ($w == 1){
       push(@genli3, $genco3{$name}[$j] - $genco3{$name}[$j - 1]);
      }
      else{
       push(@genli3, $genco3{$name}[$j + 1] - $genco3{$name}[$j]);
      }
     }
    }
    push(@genli3, ${$genco3{$name}}[$#{$genco3{$name}}] - ${$genco3{$name}}[$#{$genco3{$name}} - 1]);
   }
  }
  foreach my $name (@genes5){
   if (scalar(@{$genco5{$name}}) > 1){
    my @ind = sort { ${$genco5{$name}}[$a] <=> ${$genco5{$name}}[$b] } 0..$#{$genco5{$name}};
    @{$genco5{$name}} = @{$genco5{$name}}[@ind];
    push(@genli5, $genco5{$name}[1] - $genco5{$name}[0]);
    if (scalar(@{$genco5{$name}}) > 2){
     for (my $j = 1; $j < scalar(@{$genco5{$name}})-1; ++$j){
      if ($genco5{$name}[$j] - $genco5{$name}[$j - 1] < $genco5{$name}[$j + 1] - $genco5{$name}[$j]){
       $w = 1;
      }
      else{
       $w = 2;
      }
      if ($w == 1){
       push(@genli5, $genco5{$name}[$j] - $genco5{$name}[$j - 1]);
      }
      else{
       push(@genli5, $genco5{$name}[$j + 1] - $genco5{$name}[$j]);
      }
     }
    }
    push(@genli5, ${$genco5{$name}}[$#{$genco5{$name}}] - ${$genco5{$name}}[$#{$genco5{$name}} - 1]);
   }
  }
  foreach my $name (@genesc){
   if (scalar(@{$gencoc{$name}}) > 1){
    my @ind = sort { ${$gencoc{$name}}[$a] <=> ${$gencoc{$name}}[$b] } 0..$#{$gencoc{$name}};
    @{$gencoc{$name}} = @{$gencoc{$name}}[@ind];
    push(@genlic, $gencoc{$name}[1] - $gencoc{$name}[0]);
    if (scalar(@{$gencoc{$name}}) > 2){
     for (my $j = 1; $j < scalar(@{$gencoc{$name}})-1; ++$j){
      if ($gencoc{$name}[$j] - $gencoc{$name}[$j - 1] < $gencoc{$name}[$j + 1] - $gencoc{$name}[$j]){
       $w = 1;
      }
      else{
       $w = 2;
      }
      if ($w == 1){
       push(@genlic, $gencoc{$name}[$j] - $gencoc{$name}[$j - 1]);
      }
      else{
       push(@genlic, $gencoc{$name}[$j + 1] - $gencoc{$name}[$j]);
      }
     }
    }
    push(@genlic, ${$gencoc{$name}}[$#{$gencoc{$name}}] - ${$gencoc{$name}}[$#{$gencoc{$name}} - 1]);
   }
  }
  foreach my $name (@genesle){
   if (scalar(@{$gencole{$name}}) > 1){
    my @ind = sort { ${$gencole{$name}}[$a] <=> ${$gencole{$name}}[$b] } 0..$#{$gencole{$name}};
    @{$gencole{$name}} = @{$gencole{$name}}[@ind];
    push(@genlile, $gencole{$name}[1] - $gencole{$name}[0]);
    if (scalar(@{$gencole{$name}}) > 2){
     for (my $j = 1; $j < scalar(@{$gencole{$name}})-1; ++$j){
      if ($gencole{$name}[$j] - $gencole{$name}[$j - 1] < $gencole{$name}[$j + 1] - $gencole{$name}[$j]){
       $w = 1;
      }
      else{
       $w = 2;
      }
      if ($w == 1){
       push(@genlile, $gencole{$name}[$j] - $gencole{$name}[$j - 1]);
      }
      else{
       push(@genlile, $gencole{$name}[$j + 1] - $gencole{$name}[$j]);
      }
     }
    }
    push(@genlile, ${$gencole{$name}}[$#{$gencole{$name}}] - ${$gencole{$name}}[$#{$gencole{$name}} - 1]);
   }
  }
  foreach my $name (@genesli){
   if (scalar(@{$gencoli{$name}}) > 1){
    my @ind = sort { ${$gencoli{$name}}[$a] <=> ${$gencoli{$name}}[$b] } 0..$#{$gencoli{$name}};
    @{$gencoli{$name}} = @{$gencoli{$name}}[@ind];
    push(@genlili, $gencoli{$name}[1] - $gencoli{$name}[0]);
    if (scalar(@{$gencoli{$name}}) > 2){
     for (my $j = 1; $j < scalar(@{$gencoli{$name}})-1; ++$j){
      if ($gencoli{$name}[$j] - $gencoli{$name}[$j - 1] < $gencoli{$name}[$j + 1] - $gencoli{$name}[$j]){
       $w = 1;
      }
      else{
       $w = 2;
      }
      if ($w == 1){
       push(@genlili, $gencoli{$name}[$j] - $gencoli{$name}[$j - 1]);
      }
      else{
       push(@genlili, $gencoli{$name}[$j + 1] - $gencoli{$name}[$j]);
      }
     }
    }
    push(@genlili, ${$gencoli{$name}}[$#{$gencoli{$name}}] - ${$gencoli{$name}}[$#{$gencoli{$name}} - 1]);
   }
  }
 }
 else{
  foreach my $name (@genesii){
   if (scalar(@{$gencoii1{$name}}) > 0 and scalar(@{$gencoii2{$name}}) > 0){
    for (my $i = 0; $i <= $#{$gencoii1{$name}}; ++$i){
     my @suba = map { ${$gencoii2{$name}}[$_] - ${$gencoii1{$name}}[$i] } 0 .. $#{$gencoii2{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genliii1, $gencoii2{$name}[$i] - $gencoii1{$name}[$mini]);
     }
     else{
      push(@genliii1, $gencoii1{$name}[$i] - $gencoii2{$name}[$mini]);
     }
     push(@genlii1, $genpii1{$name}[$i]);
    }
    for (my $i = 0; $i <= $#{$gencoii2{$name}}; ++$i){
     my @suba = map { ${$gencoii1{$name}}[$_] - ${$gencoii2{$name}}[$i] } 0 .. $#{$gencoii1{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genliii2, $gencoii1{$name}[$i] - $gencoii2{$name}[$mini]);
     }
     else{
      push(@genliii2, $gencoii2{$name}[$i] - $gencoii1{$name}[$mini]);
     }
     push(@genlii2, $genpii2{$name}[$i]);
    }
   }
  }
  foreach my $name (@genes3){
   if (scalar(@{$genco31{$name}}) > 0 and scalar(@{$genco32{$name}}) > 0){
    for (my $i = 0; $i <= $#{$genco31{$name}}; ++$i){
     my @suba = map { ${$genco32{$name}}[$_] - ${$genco31{$name}}[$i] } 0 .. $#{$genco32{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genli31, $genco32{$name}[$i] - $genco31{$name}[$mini]);
     }
     else{
      push(@genli31, $genco31{$name}[$i] - $genco32{$name}[$mini]);
     }
     push(@genl31, $genp31{$name}[$i]);
    }
    for (my $i = 0; $i <= $#{$genco32{$name}}; ++$i){
     my @suba = map { ${$genco31{$name}}[$_] - ${$genco32{$name}}[$i] } 0 .. $#{$genco31{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genli32, $genco31{$name}[$i] - $genco32{$name}[$mini]);
     }
     else{
      push(@genli32, $genco32{$name}[$i] - $genco31{$name}[$mini]);
     }
     push(@genl32, $genp32{$name}[$i]);
    }
   }
  }
  foreach my $name (@genes5){
   if (scalar(@{$genco51{$name}}) > 0 and scalar(@{$genco52{$name}}) > 0){
    for (my $i = 0; $i <= $#{$genco51{$name}}; ++$i){
     my @suba = map { ${$genco52{$name}}[$_] - ${$genco51{$name}}[$i] } 0 .. $#{$genco52{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genli51, $genco52{$name}[$i] - $genco51{$name}[$mini]);
     }
     else{
      push(@genli51, $genco51{$name}[$i] - $genco52{$name}[$mini]);
     }
     push(@genl51, $genp51{$name}[$i]);
    }
    for (my $i = 0; $i <= $#{$genco52{$name}}; ++$i){
     my @suba = map { ${$genco51{$name}}[$_] - ${$genco52{$name}}[$i] } 0 .. $#{$genco51{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genli52, $genco51{$name}[$i] - $genco52{$name}[$mini]);
     }
     else{
      push(@genli52, $genco52{$name}[$i] - $genco51{$name}[$mini]);
     }
     push(@genl52, $genp52{$name}[$i]);
    }
   }
  }
  foreach my $name (@genesc){
   if (scalar(@{$gencoc1{$name}}) > 0 and scalar(@{$gencoc2{$name}}) > 0){
    for (my $i = 0; $i <= $#{$gencoc1{$name}}; ++$i){
     my @suba = map { ${$gencoc2{$name}}[$_] - ${$gencoc1{$name}}[$i] } 0 .. $#{$gencoc2{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genlic1, $gencoc2{$name}[$i] - $gencoc1{$name}[$mini]);
     }
     else{
      push(@genlic1, $gencoc1{$name}[$i] - $gencoc2{$name}[$mini]);
     }
     push(@genlc1, $genpc1{$name}[$i]);
    }
    for (my $i = 0; $i <= $#{$gencoc2{$name}}; ++$i){
     my @suba = map { ${$gencoc1{$name}}[$_] - ${$gencoc2{$name}}[$i] } 0 .. $#{$gencoc1{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genlic2, $gencoc1{$name}[$i] - $gencoc2{$name}[$mini]);
     }
     else{
      push(@genlic2, $gencoc2{$name}[$i] - $gencoc1{$name}[$mini]);
     }
     push(@genlc2, $genpc2{$name}[$i]);
    }
   }
  }
  foreach my $name (@genesle){
   if (scalar(@{$gencole1{$name}}) > 0 and scalar(@{$gencole2{$name}}) > 0){
    for (my $i = 0; $i <= $#{$gencole1{$name}}; ++$i){
     my @suba = map { ${$gencole2{$name}}[$_] - ${$gencole1{$name}}[$i] } 0 .. $#{$gencole2{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genlile1, $gencole2{$name}[$i] - $gencole1{$name}[$mini]);
     }
     else{
      push(@genlile1, $gencole1{$name}[$i] - $gencole2{$name}[$mini]);
     }
     push(@genlle1, $genple1{$name}[$i]);
    }
    for (my $i = 0; $i <= $#{$gencole2{$name}}; ++$i){
     my @suba = map { ${$gencole1{$name}}[$_] - ${$gencole2{$name}}[$i] } 0 .. $#{$gencole1{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genlile2, $gencole1{$name}[$i] - $gencole2{$name}[$mini]);
     }
     else{
      push(@genlile2, $gencole2{$name}[$i] - $gencole1{$name}[$mini]);
     }
     push(@genlle2, $genple2{$name}[$i]);
    }
   }
  }
  foreach my $name (@genesli){
   if (scalar(@{$gencoli1{$name}}) > 0 and scalar(@{$gencoli2{$name}}) > 0){
    for (my $i = 0; $i <= $#{$gencoli1{$name}}; ++$i){
     my @suba = map { ${$gencoli2{$name}}[$_] - ${$gencoli1{$name}}[$i] } 0 .. $#{$gencoli2{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genlili1, $gencoli2{$name}[$i] - $gencoli1{$name}[$mini]);
     }
     else{
      push(@genlili1, $gencoli1{$name}[$i] - $gencoli2{$name}[$mini]);
     }
     push(@genlli1, $genpli1{$name}[$i]);
    }
    for (my $i = 0; $i <= $#{$gencoli2{$name}}; ++$i){
     my @suba = map { ${$gencoli1{$name}}[$_] - ${$gencoli2{$name}}[$i] } 0 .. $#{$gencoli1{$name}};
     my $mini = 0;
     foreach my $j (1..$#suba){
      $mini = $j if abs($suba[$j])<abs($suba[$mini]);
     }
     if ($stranded{$name} eq "+"){
      push(@genlili2, $gencoli1{$name}[$i] - $gencoli2{$name}[$mini]);
     }
     else{
      push(@genlili2, $gencoli2{$name}[$i] - $gencoli1{$name}[$mini]);
     }
     push(@genlli2, $genpli2{$name}[$i]);
    }
   }
  }
 }
}

sub advance{
 if ($Aligned[$_[0]] ne "3'utr" and $Aligned[$_[0]] ne "5'utr" and $Aligned[$_[0]] ne "coding" and $Aligned[$_[0]] ne "intron" and $Aligned[$_[0]] ne "lincRNA_exon" and $Aligned[$_[0]] ne "lincRNA_intron"){
  return;
 }
 my $ub = scalar(@end) - 1;
 my $lb = 0;
 while (1 == 1){
  if ($ub <= $lb + 1){
   return;
  }
  $avg = int(($ub + $lb)/2);
  if ($Aligned[$_[0]] ne "intron" and $Aligned[$_[0]] ne "lincRNA_intron"){
   if ($start[$avg] <= $End[$_[0]] and $start[$avg + 1] > $End[$_[0]] and $Chr[$_[0]] eq $chr[$avg]){
    $y = $avg - 100;
    last;
   }
   elsif ($Chr[$_[0]] lt $chr[$avg] or ($start[$avg] > $End[$_[0]] and $Chr[$_[0]] eq $chr[$avg])){
    $ub = $avg;
   }
   elsif ($Chr[$_[0]] gt $chr[$avg] or ($start[$avg] <= $End[$_[0]] and $Chr[$_[0]] eq $chr[$avg])){
    $lb = $avg;
   }
  }
  else{
   if ($start[$avg] > $End[$_[0]] and $start[$avg - 1] <= $End[$_[0]] and $Chr[$_[0]] eq $chr[$avg]){
    $y = $avg - 100;
    last;
   }
   elsif (($start[$avg] > $End[$_[0]] and $Chr[$_[0]] eq $chr[$avg]) or $Chr[$_[0]] lt $chr[$avg]){
    $ub = $avg;
   }
   elsif (($start[$avg] <= $End[$_[0]] and $Chr[$_[0]] eq $chr[$avg]) or $Chr[$_[0]] gt $chr[$avg]){
    $lb = $avg;
   }
  }
 }
 @indices = ();
 @matches = ();
 @matchname = ();
 $t = 0;
 $counter = -1;
 $me = ($Start[$_[0]] + $End[$_[0]])/2;
 while ($y <= scalar(@end)){
  if (($reg[$y] ne "exon" and $reg[$y] ne "start_codon" and $reg[$y] ne "stop_codon") or $chr[$y] lt $Chr[$_[0]] or !$accept[$y]){
   $y++;
   next;
  }
  elsif ((($start[$y] <= $End[$_[0]] and $End[$_[0]] <= $end[$y]) or ($start[$y] <= $Start[$_[0]] and $Start[$_[0]] <= $end[$y])) and $Strand[$_[0]] eq $strand[$y] and $Chr[$_[0]] eq $chr[$y] and $Aligned[$_[0]] ne "intron" and $Aligned[$_[0]] ne "lincRNA_intron"){
   $ty = $y;
   if (scalar(@matches) == 0){
    $genname = $gena[$y];
    $ty -= 100;
   }
   push(@matches, [$start[$y], $end[$y], $y]);
   push(@matchname, $trid[$y] . "_" . $exnu[$y]);
   $y = $ty;
  }
  elsif (($start[$y] > $End[$_[0]] and $chr[$y] eq $Chr[$_[0]] and $strand[$y] eq $Strand[$_[0]]) or $chr[$y] ne $Chr[$_[0]] or $counter > -1){
   if ((($Aligned[$_[0]] eq "3'utr" and $Strand[$_[0]] eq "-") or ($Aligned[$_[0]] eq "5'utr" and $Strand[$_[0]] eq "+")) and $counter < 150){
    if ((($Aligned[$_[0]] eq "3'utr" and $reg[$y] eq "stop_codon") or ($Aligned[$_[0]] eq "5'utr" && $reg[$y] eq "start_codon")) and $gena[$y] eq $genname){
     push(@indices, [$start[$y], $trid[$y]]);
    }
    $counter++;
   }
   elsif (($Aligned[$_[0]] eq "intron" or $Aligned[$_[0]] eq "lincRNA_intron") and $chr[$y] eq $Chr[$_[0]]){
    $genname = $gena[$y];
    my @rest;
    for (my $i = $y - &min($y, 50); $i <= $y + 50; ++$i){
     if ($gena[$i] eq $genname){
      push(@rest, $i);
     }
    }
    my @begin = ();
    my @matches = ();
    foreach $res (@rest){
     if ($reg[$res] ne "exon"){
      next;
     }
     elsif ($end[$res] < $Start[$_[0]]){
      push(@begin, [$end[$res], $trid[$res], $res, int($exnu[$res])]);
     }
     elsif ($start[$res] > $End[$_[0]]){
      push(@matches, [$start[$res], $trid[$res], $res, int($exnu[$res])]);
     }
    }
    if (scalar(@begin) == 0){
     return;
    }
    my $dis = 0;
    my $te = -1;
    my $tee = 0;
    my $c = 0;
    my $sco = -1;
    for ($i = 0; $i <= $#begin; ++$i){
     for ($j = 0; $j <= $#matches; ++$j){
      if ($matches[$j][1] eq $begin[$i][1]){
       if ($rseq ne ""){
	if ($trab{$matches[$j][1]} > $sco){
         $sco = $trab{$matches[$j][1]};
	 $dis = $matches[$j][0] - $begin[$i][0] - 2;
         $te = $j;
         $tee = $i;
         $c = $matches[$te][2];
        }
       }
       elsif ($matches[$j][0] - $begin[$i][0] - 2 > $dis and ($matches[$j][3] - 1 == $begin[$i][3] or $matches[$j][3] + 1 == $begin[$i][3])){
        $dis = $matches[$j][0] - $begin[$i][0] - 2;
        $te = $j;
        $tee = $i;
        $c = $matches[$te][2];
       }
      }
     }
    }
    if ($te == -1){
     return;
    }
    if ($Aligned[$_[0]] eq "lincRNA_intron"){
     if ($Strand[$_[0]] eq "+"){
      $app = &addOn($begin[$tee][0] + 1, $matches[$te][0] - 1, $me, $beg);
     }
     else{
      $app = &addOn($begin[$tee][0] + 1, $matches[$te][0] - 1, $me, !$beg);
     }
     if ($whic == 0){
      push(@templi, $app);
      push(@templid, $dis);
      push(@genpli, $ID[$_[0]]);
     }
     elsif ($whic == 1){
      push(@templi1, $app / ($matches[$te][0] - $begin[$tee][0]));
     }
     else{
      push(@templi2, $app / ($matches[$te][0] - $begin[$tee][0]));
     }
    }
    else{
     if ($Strand[$_[0]] eq "+"){
      $app = &addOn($begin[$tee][0] + 1, $matches[$te][0] - 1, $me, $beg);
     }
     else{
      $app = &addOn($begin[$tee][0] + 1, $matches[$te][0] - 1, $me, !$beg);
     }
     if ($whic == 0){
      push(@tempi, $app);
      push(@tempid, $dis);
      push(@genpi, $ID[$_[0]]);
     }
     elsif ($whic == 1){
      push(@tempi1, $app / ($matches[$te][0] - $begin[$tee][0]));
     }
     else{
      push(@tempi2, $app / ($matches[$te][0] - $begin[$tee][0]));
     }
    }
    $genna = $trid[$c] . "_" . $exnu[$c];
    if ($whic == 0 and $Aligned[$_[0]] eq "lincRNA_intron" and !(grep(/^$genna$/, keys %gencoli))){
     push(@{$gencoli{$genna}}, $me);
     push(@genmali, $dis);
     push(@genesli, $genna);
     $stranded{$genna}	= $Strand[$_[0]];
    }
    elsif ($whic > 0 and $Aligned[$_[0]] eq "lincRNA_intron" and !(grep(/^$genna$/, keys %gencoli1))){
     if ($whic == 1){
      push(@{$gencoli1{$genna}}, $me);
      @{$gencoli2{$genna}} = ();
      push(@{$genpli1{$genna}}, $_[0]);
      @{$genpli2{$genna}} = ();
     }
     else{
      push(@{$gencoli2{$genna}}, $me);
      @{$gencoli1{$genna}} = ();
      push(@{$genpli2{$genna}}, $_[0]);
      @{$genpli1{$genna}} = ();
     }
     push(@genesli, $genna);
     $stranded{$genna}	= $Strand[$_[0]];
     push(@genmali, [$begin[$tee][0], $matches[$te][0]]);
    }
    elsif ($Aligned[$_[0]] eq "lincRNA_intron"){
     if ($whic == 0){
      push(@{$gencoli{$genna}}, $me);
     }
     elsif ($whic == 1){
      push(@{$gencoli1{$genna}}, $me);
      push(@{$genpli1{$genna}}, $_[0]);
     }
     else{
      push(@{$gencoli2{$genna}}, $me);
      push(@{$genpli2{$genna}}, $_[0]);
     }
    }
    elsif ($whic == 0 and !(grep(/^$genna$/, keys %gencoii))){
     push(@{$gencoii{$genna}}, $me);
     push(@genmaii, $dis);
     push(@genesii, $genna);
     $stranded{$genna}	= $Strand[$_[0]];
    }
    elsif ($whic > 0 and !(grep(/^$genna$/, keys %gencoii1))){
     if ($whic == 1){
      push(@{$gencoii1{$genna}}, $me);
      @{$gencoii2{$genna}} = ();
      push(@{$genpii1{$genna}}, $_[0]);
      @{$genpii2{$genna}} = ();
     }
     else{
      push(@{$gencoii2{$genna}}, $me);
      @{$gencoii1{$genna}} = ();
      push(@{$genpii2{$genna}}, $_[0]);
      @{$genpii1{$genna}} = ();
     }
     push(@genesii, $genna);
     $stranded{$genna}	= $Strand[$_[0]];
     push(@genmaii, [$begin[$tee][0], $matches[$te][0]]);
    }
    else{
     if ($whic == 0){
      push(@{$gencoii{$genna}}, $me);
     }
     elsif ($whic == 1){
      push(@{$gencoii1{$genna}}, $me);
      push(@{$genpii1{$genna}}, $_[0]);
     }
     else{
      push(@{$gencoii2{$genna}}, $me);
      push(@{$genpii2{$genna}}, $_[0]);
     }
    }
    return;
   }
   elsif (scalar(@matches) > 0){
    $c = 0;
    $len = 0;
    $te = -1;
    $tei = 0;
    $sco = -1;
    shift(@matches);
    shift(@matchname);
    while ($c < scalar(@matches)){
     if ($matches[$c][0] <= $Start[$_[0]] and $matches[$c][1] >= $End[$_[0]]){
      my @tri = split("_", $matchname[$c]);
      if ($Aligned[$_[0]] eq "5'utr" or $Aligned[$_[0]] eq "3'utr"){
       if (($Aligned[$_[0]] eq "5'utr" and $Strand[$_[0]] eq "+") or ($Aligned[$_[0]] eq "3'utr" and $Strand[$_[0]] eq "-")){
        for ($i = 0; $i <= $#indices; ++$i){
         if ($indices[$i][1] eq $tri[0]){
          if ($matches[$c][1] > $indices[$i][0]){
	   if ($rseq ne ""){
	    if ($trab{$indices[$i][1]} > $sco and $indices[$i][0] - $matches[$c][0] - 1 > 0){
	     $sco = $trab{$indices[$i][1]};
	     $len = $indices[$i][0] - $matches[$c][0] - 1;
 	     $te = $c;
             $tei = $indices[$i][0];
             $genna = $matchname[$te];
	    }
	   }
           elsif ($indices[$i][0] - $matches[$c][0] - 1 > $len and $End[$_[0]] < $indices[$i][0]){
            $len = $indices[$i][0] - $matches[$c][0] - 1;
            $te = $c;
            $tei = $indices[$i][0];
            $genna = $matchname[$te];
           }
          }
          else{
	   if ($rseq ne ""){
	    if ($trab{$indices[$i][1]} > $sco and $matches[$c][1] - $matches[$c][0] > 0){
	     $sco = $trab{$indices[$i][1]};
	     $len = $matches[$c][1] - $matches[$c][0];
             $te = $c;
             $tei = $matches[$c][1] + 1;
             $genna = $matchname[$te];
	    }
	   }
           elsif ($matches[$c][1] - $matches[$c][0] > $len and $End[$_[0]] < $indices[$i][0]){
            $len = $matches[$c][1] - $matches[$c][0];
            $te = $c;
            $tei = $matches[$c][1] + 1;
            $genna = $matchname[$te];
           }
          }
         }
        }
       }
       else{
        for ($i = 0; $i <= $#indices; ++$i){
         if ($indices[$i][1] eq $tri[0]){
	  if ($rseq ne ""){
	   if ($trab{$indices[$i][1]} > $sco and $matches[$c][1] - $indices[$i][0] - 1 > 0){
	    if ($matches[$c][0] < $indices[$i][0]){
	     $tei = $indices[$i][0];
	     $len = $matches[$c][1] - $indices[$i][0] - 1;
	    }
	    else{
	     $tei = $matches[$c][0] + 1;
	     $len = $matches[$c][1] - $matches[$c][0];
	    }
	    $sco = $trab{$indices[$i][1]};
	    $te = $c;
            $genna = $matchname[$te];
	   }
	  }
          elsif ($matches[$c][1] - $indices[$i][0] - 1 > $len and $Start[$_[0]] > $indices[$i][0]){
           if ($matches[$c][0] < $indices[$i][0]){
            $len = $matches[$c][1] - $indices[$i][0] - 1;
            $tei = $indices[$i][0];
           }
           else{
            $len = $matches[$c][1] - $matches[$c][0];
            $tei = $matches[$c][0] + 1;
           }
           $te = $c;
           $genna = $matchname[$te];
          }
         }
        }
       }
      }
      elsif ($rseq ne ""){
       if ($trab{$tri[0]} > $sco and $matches[$c][1] - $matches[$c][0] > 0){
	$sco = $trab{$tri[0]};
	$len = $matches[$c][1] - $matches[$c][0];
	$te = $c;
	$genna = $matchname[$te];
       }
      }
      elsif ($matches[$c][1] - $matches[$c][0] > $len){
       $len = $matches[$c][1] - $matches[$c][0];
       $te = $c;
       $genna = $matchname[$te];
      }
     }
     $c++;
    }
    if ($te == -1){
     return;
    }
    if ($Aligned[$_[0]] eq "3'utr"){
     if ($Strand[$_[0]] eq "-"){
      $app = &addOn($matches[$te][0], $tei - 1, $me, !$beg);
     }
     else{
      $app = &addOn($tei + 1, $matches[$te][1], $me, $beg);
     }
     $dis = $len;
     if ($whic == 0){
      push(@temp3e, $dis);
      push(@temp3, $app);
      push(@genp3, $ID[$_[0]]);
     }
     elsif ($whic == 1){
      push(@temp31, $app / $dis);
     }
     else{
      push(@temp32, $app / $dis);
     }
     if ($whic == 0 and !(grep(/^$genna$/, keys %genco3))){
      push(@{$genco3{$genna}}, $me);
      push(@genma3, $dis);
      push(@genes3, $genna);
      $stranded{$genna} = $Strand[$_[0]];
     }
     elsif ($whic > 0 and !(grep(/^$genna$/, keys %genco31))){
      if ($whic == 1){
       push(@{$genco31{$genna}}, $me);
       @genco32{$genna} = ();
       push(@{$genp31{$genna}}, $_[0]);
       @genp32{$genna} = ();
      }
      else{
       push(@{$genco32{$genna}}, $me);
       @genco31{$genna} = ();
       push(@{$genp32{$genna}}, $_[0]);
       @genp31{$genna} = ();
      }
      if ($Strand[$_[0]] eq "-"){
       push(@genma3, [$matches[$te][0], $tei - 1]);
      }
      else{
       push(@genma3, [$tei + 1, $matches[$te][1]]);
      }
      push(@genes3, $genna);
      $stranded{$genna} = $Strand[$_[0]];
     }
     else{
      if ($whic == 0){
       push(@{$genco3{$genna}}, $me);
      }
      elsif ($whic == 1){
       push(@{$genco31{$genna}}, $me);
       push(@{$genp31{$genna}}, $_[0]);
      }
      else{
       push(@{$genco32{$genna}}, $me);
       push(@{$genp32{$genna}}, $_[0]);
      }
     }
    }
    elsif ($Aligned[$_[0]] eq "5'utr"){
     if ($Strand[$_[0]] eq "+"){
      $app = &addOn($matches[$te][0], $tei - 1, $me, $beg);
     }
     else{
      $app = &addOn($tei + 1, $matches[$te][1], $me, !$beg);
     }
     $dis = $len;
     if ($whic == 0){
      push(@temp5e, $dis);
      push(@temp5, $app);
      push(@genp5, $ID[$_[0]]);
     }
     elsif ($whic == 1){
      push(@temp51, $app / $dis);
     }
     else{
      push(@temp52, $app / $dis);
     }
     if ($whic == 0 and !(grep(/^$genna$/, keys %genco5))){
      push(@{$genco5{$genna}}, $me);
      push(@genma5, $dis);
      push(@genes5, $genna);
      $stranded{$genna} = $Strand[$_[0]];
     }
     elsif ($whic > 0 and !(grep(/^$genna$/, keys %genco51))){
      if ($whic == 1){
       push(@{$genco51{$genna}}, $me);
       @genco52{$genna} = ();
       push(@{$genp51{$genna}}, $_[0]);
       @genco52{$genna} = ();
      }
      else{
       push(@{$genco52{$genna}}, $me);
       @genco51{$genna} = ();
       push(@{$genp52{$genna}}, $_[0]);
       @genp51{$genna} = ();
      }
      if ($Strand[$_[0]] eq "+"){
       push(@genma5, [$matches[$te][0], $tei-1]);
      }
      else{
       push(@genma5, [$tei+1, $matches[$te][1]]);
      }
      push(@genes5, $genna);
      $stranded{$genna}	= $Strand[$_[0]];
     }
     else{
      if ($whic == 0){
       push(@{$genco5{$genna}}, $me);
      }
      elsif ($whic == 1){
       push(@{$genco51{$genna}}, $me);
       push(@{$genp51{$genna}}, $_[0]);
      }
      else{
       push(@{$genco52{$genna}}, $me);
       push(@{$genp52{$genna}}, $_[0]);
      }
     }
    }
    elsif ($Aligned[$_[0]] eq "coding"){
     if ($Strand[$_[0]] eq "+"){
      $app = &addOn($matches[$te][0], $matches[$te][1], $me, $beg);
     }
     else{
      $app = &addOn($matches[$te][0], $matches[$te][1], $me, !$beg);
     }
     $dis = $matches[$te][1] - $matches[$te][0];
     if ($whic == 0){
      push(@tempce, $dis);
      push(@tempc, $app);
      push(@genpc, $ID[$_[0]]);
     }
     elsif ($whic == 1){
      push(@tempc1, $app / $dis);
     }
     else{
      push(@tempc2, $app / $dis);
     }
     if ($whic == 0 and !(grep(/^$genna$/, keys %gencoc))){
      push(@{$gencoc{$genna}}, $me);
      push(@genmac, $dis);
      push(@genesc, $genna);
      $stranded{$genna}	= $Strand[$_[0]];
     }
     elsif ($whic > 0 and !(grep(/^$genna$/, keys %gencoc1))){
      if ($whic == 1){
       push(@{$gencoc1{$genna}}, $me);
       @gencoc2{$genna} = ();
       push(@{$genpc1{$genna}}, $_[0]);
       @genpc2{$genna} = ();
      }
      else{
       push(@{$gencoc2{$genna}}, $me);
       @gencoc1{$genna} = ();
       push(@{$genpc2{$genna}}, $_[0]);
       @genpc1{$genna} = ();
      }
      push(@genmac, [$matches[$te][0], $matches[$te][1]]);
      push(@genesc, $genna);
      $stranded{$genna}	= $Strand[$_[0]];
     }
     else{
      if ($whic == 0){
       push(@{$gencoc{$genna}}, $me);
      }
      elsif ($whic == 1){
       push(@{$gencoc1{$genna}}, $me);
       push(@{$genpc1{$genna}}, $_[0]);
      }
      else{
       push(@{$gencoc2{$genna}}, $me);
       push(@{$genpc2{$genna}}, $_[0]);
      }
     }
    }
    else{
     if ($Strand[$_[0]] eq "+"){
      $app = &addOn($matches[$te][0], $matches[$te][1], $me, $beg);
     }
     else{
      $app = &addOn($matches[$te][0], $matches[$te][1], $me, !$beg);
     }
     $dis = $matches[$te][1] - $matches[$te][0];
     if ($whic == 0){
      push(@templed, $dis);
      push(@temple, $app);
      push(@genple, $ID[$_[0]]);
     }
     elsif ($whic == 1){
      push(@temple1, $app / $dis);
     }
     else{
      push(@temple2, $app / $dis);
     }
     if ($whic == 0 and !(grep(/^$genna$/, keys %gencole))){
      push(@{$gencole{$genna}}, $me);
      push(@genmale, $dis);
      push(@genesle, $genna);
      $stranded{$genna}	= $Strand[$_[0]];
     }
     elsif ($whic > 0 and !(grep(/^$genna$/, keys %gencole1))){
      if ($whic == 1){
       push(@{$gencole1{$genna}}, $me);
       @gencole2{$genna} = ();
       push(@{$genple1{$genna}}, $_[0]);
       @genple2{$genna} = ();
      }
      else{
       push(@{$gencole2{$genna}}, $me);
       @gencole1{$genna} = ();
       push(@{$genple2{$genna}}, $_[0]);
       @genple1{$genna} = ();
      }
      push(@genmale, [$matches[$te][0], $matches[$te][1]]);
      push(@genesle, $genna);
      $stranded{$genna}	= $Strand[$_[0]];
     }
     else{
      if ($whic == 0){
       push(@{$gencole{$genna}}, $me);
      }
      elsif ($whic == 1){
       push(@{$gencole1{$genna}}, $me);
       push(@{$genple1{$genna}}, $_[0]);
      }
      else{
       push(@{$gencole2{$genna}}, $me);
       push(@{$genple2{$genna}}, $_[0]);
      }
     }
    }
    return;
   }
   else{
    return;
   }
  }
  elsif ((($Aligned[$_[0]] eq "3'utr" and $reg[$y] eq "stop_codon" and $strand[$y] eq "+") or ($Aligned[$_[0]] eq "5'utr" and $reg[$y] eq "start_codon" and $strand[$y] eq "-")) and $gena[$y] eq $genname and $end[$y] < $Start[$_[0]] and $Strand[$_[0]] eq $strand[$y] and $chr[$y] eq $Chr[$_[0]]){
   push(@indices, [$end[$y], $trid[$y]]);
  }
  $y++;
 }
 return;
}

if (substr($gtf, -3) eq ".gz"){
 @gtf = `gunzip -c $gtf`;
}
else{
 open(GTF, $gtf);
 @gtf = <GTF>;
 close(GTF);
}
foreach (@gtf){
 chomp;
 @secs = split("\t");
 if ($secs[0] =~ /\d*|Y|X|MT|chr.*/ and $secs[8] =~ /exon_number "*(\d*)"*/){
  push(@chr, $secs[0]);
  push(@reg, $secs[2]);
  push(@start, $secs[3]);
  push(@end, $secs[4]);
  push(@strand, $secs[6]);
  $secs[8] =~ /exon_number "*(\d*)"*/;
  push(@exnu, $1);
  $secs[8] =~ /gene_name "([^"]*)"/;
  push(@gena, $1);
  $secs[8] =~ /transcript_id "([\w]*[\.]*[\d]*)/;
  push(@trid, $1);
  if ($secs[8] =~ /transcript_status "([^"]*)"/){
   if ($strict and ($1 ne "KNOWN" or $secs[8] =~ "cds_start_NF" or $secs[8] =~ "cds_end_NF" or $secs[8] =~ "mRNA_start_NF" or $secs[8] =~ "mRNA_end_NF")){
    push(@accept, 0);
   }
   else{
    push(@accept, 1);
   }
  }
  else{
   push(@accept, 1);
  }
 }
}
#close(GTF);
s/chr// for @chr;
s/MT/M/ for @chr;
@in = sort { $start[$a] <=> $start[$b] } 0..$#start;
@reg = @reg[@in];
@chr = @chr[@in];
@start = @start[@in];
@strand = @strand[@in];
@end = @end[@in];
@exnu = @exnu[@in];
@gena = @gena[@in];
@trid = @trid[@in];
@accept = @accept[@in];
@in = sort { $chr[$a] cmp $chr[$b] } 0..$#chr;
@reg = @reg[@in];
@chr = @chr[@in];
@start = @start[@in];
@strand = @strand[@in];
@end = @end[@in];
@exnu = @exnu[@in];
@gena = @gena[@in];
@trid = @trid[@in];
@accept = @accept[@in];
@Aligned = ();
@Chr = ();
@Strand = ();
@Start = ();
@End = ();
if ($yn eq "y" or $yn eq "yes"){
 open(F2, $fn2 . ".clusters.csv");
 @f2 = <F2>;
 close(F2);
 @secs = split(",", $f2[0]);
 shift(@f2);
 for ($i = 0; $i <= $#secs; ++$i){
  if ($secs[$i] eq "Chromosome" or $secs[$i] eq "Chr"){
   $chr = $i;
  }
  elsif ($secs[$i] eq "Strand"){
   $str = $i;
  }
  elsif ($secs[$i] eq "Start"){
   $st = $i;
  }
  elsif ($secs[$i] eq "End"){
   $en = $i;
  }
  elsif ($secs[$i] eq "TranscriptLocation" or $secs[$i] eq "Aligned to"){
   $al = $i;
  }
 }
 foreach $i(@f2){
  @secs = split(",", $i);
  push(@Aligned, $secs[$al]);
  push(@Chr, $secs[$chr]);
  push(@Start, $secs[$st]);
  push(@Strand, $secs[$str]);
  push(@End, $secs[$en]);
 }
 @in = sort { $Start[$a] <=> $Start[$b] } 0..$#Start;
 @Aligned = @Aligned[@in];
 @Chr = @Chr[@in];
 @Start = @Start[@in];
 @Strand = @Strand[@in];
 @End = @End[@in];
 s/chr// for @Chr;
 s/MT/M/ for @Chr;
 @in = sort { $Chr[$a] cmp $Chr[$b] } 0..$#Chr;
 @Aligned = @Aligned[@in];
 @Chr = @Chr[@in];
 @Start = @Start[@in];
 @Strand = @Strand[@in];
 @End = @End[@in];
 $whic = 2;
 for ($p = 0; $p <= $#Chr; ++$p){
  &advance($p);
 }
}
@Aligned = ();
@Chr = ();
@Strand = ();
@Start = ();
@End = ();
@ID = ();
open(F1, $fn . ".clusters.csv");
@f1 = <F1>;
close(F1);
@secs = split(",", $f1[0]);
$temp = shift(@f1);
for ($i = 0; $i <= $#secs; ++$i){
 if ($secs[$i] =~ /^Chr/){
  $chr = $i;
 }
 elsif ($secs[$i] =~ /Strand/){
  $str = $i;
 }
 elsif ($secs[$i] =~ /Start/){
  $st = $i;
 }
 elsif ($secs[$i] =~ /End/){
  $en = $i;
 }
 elsif ($secs[$i] =~ /TranscriptLocation/ or $secs[$i] =~ /Aligned to/){
  $al = $i;
 }
 elsif ($secs[$i] =~ /ID/){
  $id = $i;
 }
}
foreach $i(@f1){
 @secs = split(",", $i);
 push(@Aligned, $secs[$al]);
 push(@Chr, $secs[$chr]);
 push(@Start, $secs[$st]);
 push(@Strand, $secs[$str]);
 push(@End, $secs[$en]);
 push(@ID, $secs[$id]);
}
@in = sort { $Start[$a] <=> $Start[$b] } 0..$#Start;
@Aligned = @Aligned[@in];
@Chr = @Chr[@in];
@Start = @Start[@in];
@Strand = @Strand[@in];
@End = @End[@in];
@ID = @ID[@in];
s/chr// for @Chr;
s/M/MT/ for @Chr;
@in = sort { $Chr[$a] cmp $Chr[$b] } 0..$#Chr;
@Aligned = @Aligned[@in];
@Chr = @Chr[@in];
@Start = @Start[@in];
@Strand = @Strand[@in];
@End = @End[@in];
@ID = @ID[@in];
if ($yn eq "y"){
 $whic = 1;
}
else{
 $whic = 0;
}
for ($p = 0; $p <= $#Chr; ++$p){
 &advance($p);
}
&manage();
$" = ",";
if ($yn eq "y"){
 $file = ">tmp_" . $fn . "_" . $fn2 . "_spatial.csv";
}
else{
 $file = ">tmp_" . $fn . "_spatial.csv";
}
open(OUT, $file);
if ($yn eq "y" or $yn eq "yes"){
 my $temp31 = "@temp31";
 my $temp32 = "@temp32";
 my $temp51 = "@temp51";
 my $temp52 = "@temp52";
 my $tempc1 = "@tempc1";
 my $tempc2 = "@tempc2";
 my $tempi1 = "@tempi1";
 my $tempi2 = "@tempi2";
 my $temple1 = "@temple1";
 my $temple2 = "@temple2";
 my $templi1 = "@templi1";
 my $templi2 = "@templi2";
 my $genli31 = "@genli31";
 my $genli32 = "@genli32";
 my $genli51 = "@genli51";
 my $genli52 = "@genli52";
 my $genlic1 = "@genlic1";
 my $genlic2 = "@genlic2";
 my $genliii1 = "@genliii1";
 my $genliii2 = "@genliii2";
 my $genlile1 = "@genlile1";
 my $genlile2 = "@genlile2";
 my $genlili1 = "@genlili1";
 my $genlili2 = "@genlili2";
 my @genma31 = ();
 for (my $i = 0; $i <= $#genma3; $i++){
  push(@genma31, $genma3[$i][0]);
 }
 my @genma32 = ();
 for (my $i = 0; $i <= $#genma3; $i++){
  push(@genma32, $genma3[$i][1]);
 }
 my $genma31 = "@genma31";
 my $genma32 = "@genma32";
 my @genma51 = ();
 for (my $i = 0; $i <= $#genma5; $i++){
  push(@genma51, $genma5[$i][0]);
 }
 my @genma52 = ();
 for (my $i = 0; $i <= $#genma5; $i++){
  push(@genma52, $genma5[$i][1]);
 }
 my $genma51 = "@genma51";
 my $genma52 = "@genma52";
 my @genmac1 = ();
 for (my $i = 0; $i <= $#genmac; $i++){
  push(@genmac1, $genmac[$i][0]);
 }
 my @genmac2 = ();
 for (my $i = 0; $i <= $#genmac; $i++){
  push(@genmac2, $genmac[$i][1]);
 }
 my $genmac1 = "@genmac1";
 my $genmac2 = "@genmac2";
 my @genmaii1 = ();
 for (my $i = 0; $i <= $#genmaii; $i++){
  push(@genmaii1, $genmaii[$i][0]);
 }
 my @genmaii2 = ();
 for (my $i = 0; $i <= $#genmaii; $i++){
  push(@genmaii2, $genmaii[$i][1]);
 }
 my $genmaii1 = "@genmaii1";
 my $genmaii2 = "@genmaii2";
 my @genmale1 = ();
 for (my $i = 0; $i <= $#genmale; $i++){
  push(@genmale1, $genmale[$i][0]);
 }
 my @genmale2 = ();
 for (my $i = 0; $i <= $#genmale; $i++){
  push(@genmale2, $genmale[$i][1]);
 }
 my $genmale1 = "@genmale1";
 my $genmale2 = "@genmale2";
 my @genmali1 = ();
 for (my $i = 0; $i <= $#genmali; $i++){
  push(@genmali1, $genmali[$i][0]);
 }
 my @genmali2 = ();
 for (my $i = 0; $i <= $#genmali; $i++){
  push(@genmali2, $genmali[$i][1]);
 }
 my $genmali1 = "@genmali1";
 my $genmali2 = "@genmali2";
 my @genco31 = ();
 my @genco51 = ();
 my @gencoc1 = ();
 my @gencoii1 = ();
 my @gencole1 = ();
 my @gencoli1 = ();
 my @genco32 = ();
 my @genco52 = ();
 my @gencoc2 = ();
 my @gencoii2 = ();
 my @gencole2 = ();
 my @gencoli2 = ();
 $" = ":";
 foreach $gene(@genes3){
  my $t = "@{$genco31{$gene}}";
  push(@genco31, $gene . ":" . $t);
 }
 foreach $gene(@genes5){
  my $t = "@{$genco51{$gene}}";
  push(@genco51, $gene . ":" . $t);
 }
 foreach $gene(@genesc){
  my $t = "@{$gencoc1{$gene}}";
  push(@gencoc1, $gene . ":" . $t);
 }
 foreach $gene(@genesii){
  my $t = "@{$gencoii1{$gene}}";
  push(@gencoii1, $gene . ":" . $t);
 }
 foreach $gene(@genesle){
  my $t = "@{$gencole1{$gene}}";
  push(@gencole1, $gene . ":" . $t);
 }
 foreach $gene(@genesli){
  my $t = "@{$gencoli1{$gene}}";
  push(@gencoli1, $gene . ":" . $t);
 }
 foreach $gene(@genes3){
  my $t = "@{$genco32{$gene}}";
  push(@genco32, $gene . ":" . $t);
 }
 foreach $gene(@genes5){
  my $t = "@{$genco5{$gene}}";
  push(@genco52, $gene . ":" . $t);
 }
 foreach $gene(@genesc){
  my $t = "@{$gencoc2{$gene}}";
  push(@gencoc2, $gene . ":" . $t);
 }
 foreach $gene(@genesii){
  my $t = "@{$gencoii2{$gene}}";
  push(@gencoii2, $gene . ":" . $t);
 }
 foreach $gene(@genesle){
  my $t = "@{$gencole2{$gene}}";
  push(@gencole2, $gene . ":" . $t);
 }
 foreach $gene(@genesli){
  my $t = "@{$gencoli2{$gene}}";
  push(@gencoli2, $gene . ":" . $t);
 }
 $" = ",";
 my $genco31 = "@genco31";
 my $genco51 = "@genco51";
 my $gencoc1 = "@gencoc1";
 my $gencoii1 = "@gencoii1";
 my $gencole1 = "@gencole1";
 my $gencoli1 = "@gencoli1";
 my $genco32 = "@genco32";
 my $genco52 = "@genco52";
 my $gencoc2 = "@gencoc2";
 my $gencoii2 = "@gencoii2";
 my $gencole2 = "@gencole2";
 my $gencoli2 = "@gencoli2";
 my $ma = &max(scalar(@temp31), scalar(@temp32), scalar(@temp51), scalar(@temp52), scalar(@tempc1), scalar(@tempc2), scalar(@tempi1), scalar(@tempi2), scalar(@temple1), scalar(@temple2), scalar(@templi1), scalar(@templi2), scalar(@genli31), scalar(@genli32), scalar(@genli51), scalar(@genli52), scalar(@genlic1), scalar(@genlic2), scalar(@genlii1), scalar(@genlii2), scalar(@genlile1), scalar(@genlile2), scalar(@genlili1), scalar(@genlili2), scalar(@genma31), scalar(@genma32), scalar(@genma51), scalar(@genma52), scalar(@genmac1), scalar(@genmac2), scalar(@genmaii1), scalar(@genmaii2), scalar(@genmale1), scalar(@genmale2), scalar(@genmali1), scalar(@genmali2), scalar(@genco31), scalar(@genco32), scalar(@genco51), scalar(@genco52), scalar(@gencoc1), scalar(@gencoc2), scalar(@gencoii1), scalar(@gencoii2), scalar(@gencole1), scalar(@gencole2), scalar(@gencoli1), scalar(@gencoli2));
 my $diff = $ma - scalar(@temp31);
 for (my $i = 0; $i < $diff; ++$i){
  $temp31 .= ",";
 }
 my $genl31 = "@genl31";
 my $genl32 = "@genl32";
 my $genl51 = "@genl51";
 my $genl52 = "@genl52";
 my $genlc1 = "@genlc1";
 my $genlc2 = "@genlc2";
 my $genlii1 = "@genlii1";
 my $genlii2 = "@genlii2";
 my $genlle1 = "@genlle1";
 my $genlle2 = "@genlle2";
 my $genlli1 = "@genlli1";
 my $genlli2 = "@genlli2";
 print OUT "$fn\n";
 print OUT "$fn2\n";
 print OUT "$temp31\n";
 print OUT "$temp32\n";
 print OUT "$temp51\n";
 print OUT "$temp52\n";
 print OUT "$tempc1\n";
 print OUT "$tempc2\n";
 print OUT "$tempi1\n";
 print OUT "$tempi2\n";
 print OUT "$temple1\n";
 print OUT "$temple2\n";
 print OUT "$templi1\n";
 print OUT "$templi2\n";
 print OUT "$genli31\n";
 print OUT "$genli32\n";
 print OUT "$genli51\n";
 print OUT "$genli52\n";
 print OUT "$genlic1\n";
 print OUT "$genlic2\n";
 print OUT "$genliii1\n";
 print OUT "$genliii2\n";
 print OUT "$genlile1\n";
 print OUT "$genlile2\n";
 print OUT "$genlili1\n";
 print OUT "$genlili2\n";
 print OUT "$genma31\n";
 print OUT "$genma32\n";
 print OUT "$genma51\n";
 print OUT "$genma52\n";
 print OUT "$genmac1\n";
 print OUT "$genmac2\n";
 print OUT "$genmaii1\n";
 print OUT "$genmaii2\n";
 print OUT "$genmale1\n";
 print OUT "$genmale2\n";
 print OUT "$genmali1\n";
 print OUT "$genmali2\n";
 print OUT "$genco31\n";
 print OUT "$genco32\n";
 print OUT "$genco51\n";
 print OUT "$genco52\n";
 print OUT "$gencoc1\n";
 print OUT "$gencoc2\n";
 print OUT "$gencoii1\n";
 print OUT "$gencoii2\n";
 print OUT "$gencole1\n";
 print OUT "$gencole2\n";
 print OUT "$gencoli1\n";
 print OUT "$gencoli2\n";
 print OUT "$genl31\n";
 print OUT "$genl32\n";
 print OUT "$genl51\n";
 print OUT "$genl52\n";
 print OUT "$genlc1\n";
 print OUT "$genlc2\n";
 print OUT "$genlii1\n";
 print OUT "$genlii2\n";
 print OUT "$genlle1\n";
 print OUT "$genlle2\n";
 print OUT "$genlli1\n";
 print OUT "$genlli2\n";
}
else{
 my $temp3 = "@temp3";
 my $temp5 = "@temp5";
 my $tempc = "@tempc";
 my $tempi = "@tempi";
 my $temple = "@temple";
 my $templi = "@templi";
 my $temp3e = "@temp3e";
 my $temp5e = "@temp5e";
 my $tempce = "@tempce";
 my $tempid = "@tempid";
 my $templed = "@templed";
 my $templid = "@templid";
 my $genli3 = "@genli3";
 my $genli5 = "@genli5";
 my $genlic = "@genlic";
 my $genliii = "@genliii";
 my $genlile = "@genlile";
 my $genlili = "@genlili";
 my $genma3 = "@genma3";
 my $genma5 = "@genma5";
 my $genmac = "@genmac";
 my $genmaii = "@genmaii";
 my $genmale = "@genmale";
 my $genmali = "@genmali";
 my @genco3 = ();
 my @genco5 = ();
 my @gencoc = ();
 my @gencoii = ();
 my @gencole = ();
 my @gencoli = ();
 $" = ":";
 foreach $gene(@genes3){
  my $t = "@{$genco3{$gene}}";
  push(@genco3, $gene . ":" . $t);
 }
 foreach $gene(@genes5){
  my $t = "@{$genco5{$gene}}";
  push(@genco5, $gene . ":" . $t);
 }
 foreach $gene(@genesc){
  my $t = "@{$gencoc{$gene}}";
  push(@gencoc, $gene . ":" . $t);
 }
 foreach $gene(@genesii){
  my $t = "@{$gencoii{$gene}}";
  push(@gencoii, $gene . ":" . $t);
 }
 foreach $gene(@genesle){
  my $t = "@{$gencole{$gene}}";
  push(@gencole, $gene . ":" . $t);
 }
 foreach $gene(@genesli){
  my $t = "@{$gencoli{$gene}}";
  push(@gencoli, $gene . ":" . $t);
 }
 $" = ",";
 my $genco3 = "@genco3";
 my $genco5 = "@genco5";
 my $gencoc = "@gencoc";
 my $gencoii = "@gencoii";
 my $gencole = "@gencole";
 my $gencoli = "@gencoli";
 my $genp3 = "@genp3";
 my $genp5 = "@genp5";
 my $genpc = "@genpc";
 my $genpi = "@genpi";
 my $genple = "@genple";
 my $genpli = "@genpli";
 my $ma = &max(scalar(@temp3), scalar(@temp5), scalar(@tempc), scalar(@tempi), scalar(@temple), scalar(@templi), scalar(@temp3e), scalar(@temp5e), scalar(@tempce), scalar(@tempid), scalar(@templed), scalar(@templid), scalar(@genli3), scalar(@genli5), scalar(@genlic), scalar(@genliii), scalar(@genlile), scalar(@genlili), scalar(@genma3), scalar(@genma5), scalar(@genmac), scalar(@genmaii), scalar(@genmale), scalar(@genmali), scalar(@genco3), scalar(@genco5), scalar(@gencoc), scalar(@gencoii), scalar(@gencole), scalar(@gencoli));
 my $diff = $ma - scalar(@temp3);
 for (my $i = 0; $i < $diff; ++$i){
  $temp3 .= ",";
 }
 print OUT "$fn\n";
 print OUT "$temp3\n";
 print OUT "$temp5\n";
 print OUT "$tempc\n";
 print OUT "$tempi\n";
 print OUT "$temple\n";
 print OUT "$templi\n";
 print OUT "$temp3e\n";
 print OUT "$temp5e\n";
 print OUT "$tempce\n";
 print OUT "$tempid\n";
 print OUT "$templed\n";
 print OUT "$templid\n";
 print OUT "$genli3\n";
 print OUT "$genli5\n";
 print OUT "$genlic\n";
 print OUT "$genliii\n";
 print OUT "$genlile\n";
 print OUT "$genlili\n";
 print OUT "$genma3\n";
 print OUT "$genma5\n";
 print OUT "$genmac\n";
 print OUT "$genmaii\n";
 print OUT "$genmale\n";
 print OUT "$genmali\n";
 print OUT "$genco3\n";
 print OUT "$genco5\n";
 print OUT "$gencoc\n";
 print OUT "$gencoii\n";
 print OUT "$gencole\n";
 print OUT "$gencoli\n";
 print OUT "$genp3\n";
 print OUT "$genp5\n";
 print OUT "$genpc\n";
 print OUT "$genpi\n";
 print OUT "$genple\n";
 print OUT "$genpli\n";
}
