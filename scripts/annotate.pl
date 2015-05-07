#!/usr/bin/perl

use IO::Handle;

if ($#ARGV == -1){
 print("\nProgram:\tannotate.pl (v2.7)\nAuthor:\tNicholas Jacobs (nickcjacobs\@gmail.com)\nSummary:\tAnnotates feature files.\n\nUsage:\tperl annotate.pl [OPTIONS]\n\nOptions:\n\t-p\tFile containing ordered priorities for annotation categories. Optional.\n\n\t-i\tInput file in bed format.\n\n\t-oi\tInput file in groups.csv or clusters.csv format.\n\n\t-bi\tInput file in BAM format.\n\n\t-bed12\tInput file has split reads.\n\n\t-g\tGencode gtf file.\n\n\t-r\tRepeat feature file in bed format. Optional.\n\n\t-s\tFile from UCSC Genome Browser containing repName and repClass data. Required if -r option is used.\n\n\t-all\tWrite out results for all transcripts overlapped by entries in the input file.\n\n\t-strict\tOnly annotates to transcripts in which we have confidence.\n\n\t-a\tFile containing transcript-level FPKM data from Cufflinks. Optional.\n\n\t-c\tAnnotates to the closest protein coding gene within specified distance instead of categories below a threshold (specified by an empty line in the priorities file. Optional.\n");
 die "\n";
}
$" = "";
my $pri;
my $fn1;
my $fn2;
my $fn3;
my $fn4;
my $cluster = 0;
my $bam = 0;
my $bed12 = 0;
my $verbose = 0;
my $strict = 0;
my $repeat = 0;
my %trab;
my $rseq = "";
my $closest = 0;
for ($i = 0; $i <= $#ARGV; ++$i){
 if ($ARGV[$i] eq "-p"){
  $pri = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-i"){
  $fn1 = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-g"){
  $fn2 = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-r"){
  $fn3 = $ARGV[$i + 1];
  $repeat = 1;
 }
 elsif ($ARGV[$i] eq "-s"){
  $fn4 = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-oi"){
  $cluster = 1;
  $fn1 = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-bi"){
  $bam = 1;
  $fn1 = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-other"){
  $cluster = 1;
 }
 elsif ($ARGV[$i] eq "-bam"){
  $bam = 1;
 }
 elsif ($ARGV[$i] eq "-bed12"){
  $bed12 = 1;
 }
 elsif ($ARGV[$i] eq "-all"){
  $verbose = 1;
 }
 elsif ($ARGV[$i] eq "-strict"){
  $strict = 1;
 }
 elsif ($ARGV[$i] eq "-a"){
  $rseq = $ARGV[$i + 1];
  open(RS, $rseq);
  while(<RS>){
   chomp;
   chop($_) if ($_ =~ m/\r$/);
   my @secs = split;
   if ($secs[0] ne "tracking_id"){
    $trab{$secs[0]} = $secs[9];
   }
  }
 }
 elsif ($ARGV[$i] eq "-c"){
  $closest = $ARGV[$i + 1];
 }
}
if (!$cluster){
 if ($bed12){
  print("Chr,Strand,Start,End,Start2,End2,Aligned to,GeneName,AnnotationSource,ID\n");
 }
 else{
  print("Chr,Strand,Start,End,Aligned to,GeneName,AnnotationSource,ID\n");
 }
}
my $ctob;
my %prin;
my @order;
my $clust = "";
my %anns;
my %uend;
my %tend;
my %eend;
my %tstart;
my %estart;
my %tmatch;
my %ematch;
my %umatch;
my $start;
my $end;
my @anli = ("start_codon", "stop_codon", "3'utr-coding", "coding-3'utr", "3'utr-intron", "intron-3'utr", "5'utr-coding", "coding-5'utr", "5'utr-intron", "intron-5'utr", "coding-intron", "intron-coding", "lincRNA_exon-lincRNA_intron", "lincRNA_intron-lincRNA_exon", "coding-coding", "lincRNA_exon-lincRNA_exon");
my $prli;
my $length;
my %clan;
my $break = 0;
my $insert;
my $leave = 0;
open(RF, $fn1);
if ($cluster){
 @order = (0, 0, 0, 0, 0);
 if ($bed12){
  push(@order, (0, 0));
 }
 $secs = readline(RF);
 chomp($secs);
 chop($secs) if ($secs =~ m/\r$/);
 my @secs = split(",", $secs);
 print($secs . ",Aligned to,GeneName,AnnotationSource\n");
 for (my $i = 0; $i <= $#secs; ++$i){
  if ($secs[$i] =~ /Strand/){
   $order[3] = $i;
  }
  elsif ($secs[$i] =~ /Start/){
   $order[0] = $i;
  }
  elsif ($secs[$i] =~ /End/){
   $order[1] = $i;
  }
  elsif ($secs[$i] =~ /ID/){
   $order[2] = $i;
  }
  elsif ($secs[$i] =~ /Start2/){
   $order[4] = $i;
  }
  elsif ($secs[$i] =~ /End2/){
   $order[5] = $i;
  }
 }
}
#my $header;
#if ($bam){
# $header = `samtools view -H -b $fn1`;
#}
while (!eof(RF)){
# if ($bam and $leave){
#  $ctob = $header . "\n";
# }
# else{
  $ctob = "";
# }
 $leave = 1;
 for (my $i = 0; $i < 500000; ++$i){
  if (eof(RF)){
   last;
  }
  $secs = readline(RF);
  chomp($secs);
  chop($secs) if ($secs =~ m/\r$/);
  if ($cluster){
   @secs = split(",", $secs);
   if ($bed12){
    if ($secs[$order[4]] eq ""){
     $ctob .= $secs[0] . "\t" . $secs[$order[0]] . "\t" . $secs[$order[1]] . "\t" . $secs[$order[2]] . "\t1\t" . $secs[$order[3]] . "\t" . $secs[$order[0]] . "\t" . $secs[$order[1]] . "\t255,0,0\t1\t" . ($secs[$order[1]] - $secs[$order[0]]) . "\t0\n";
    }
    else{
     $ctob .= $secs[0] . "\t" . $secs[$order[0]] . "\t" . $secs[$order[5]] . "\t" . $secs[$order[2]] . "\t1\t" . $secs[$order[3]] . "\t" . $secs[$order[0]] . "\t" . $secs[$order[5]] . "\t255,0,0\t2\t" . ($secs[$order[1]] - $secs[$order[0]]) . "," . ($secs[$order[5]] - $secs[$order[4]]) . "\t0," . ($secs[$order[4]] - $secs[$order[0]]) . "\n";
    }
   }
   else{
    $ctob .= $secs[0] . "\t" . $secs[$order[0]] . "\t" . $secs[$order[1]] . "\t" . $secs[$order[2]] . "\t1\t" . $secs[$order[3]] . "\n";
   }
   $prin{$secs[$order[2]]} = $secs;
  }
  else{
   $ctob .= "$secs\n";
  }
 }
 $ctob = substr($ctob, 0, -1);
 if ($pid = open(CHILD, "-|")) {
  my $x = 0;
 }
 else {
  die "cannot fork: $!" unless defined $pid;
  STDOUT->autoflush(1);
  if ($bam){
   if ($bed12){
    open(W, "| intersectBed -wao -s -split -bed -abam stdin -b $fn2");
   }
   else{
    open(W, "| intersectBed -wao -s -bed -abam stdin -b $fn2");
   }
  }
  else{
   if ($bed12){
    open(W, "| intersectBed -wao -s -split -a stdin -b $fn2");
   }
   else{
    open(W, "| intersectBed -wao -s -a stdin -b $fn2");
   }
  }
  print W $ctob;
  close(W);
  exit;
 }
 open(Y, $pri);
 while (<Y>){
  chomp;
  chop($_) if ($_ =~ m/\r$/);
  if ($_ eq "lincRNA"){
   push(@anli, "lincRNA_exon");
   push(@anli, "lincRNA_intron");
  }
  elsif ($_ eq "" and $closest){
   $break = scalar(@anli);
  }
  else{
   push(@anli, $_);
  }
 }
 close(Y);
 if ($closest and !$break){
  $break = scalar(@anli);
 }
 if ($bed12 or $bam){
  $fix = 6;
 }
 if ($bed12 or $bam){
  $insert = "chr1\t1\t10\tGA\t1\t+\t1\t10\t255,0,0\t1\t0\t0\t.\t.\t.\t-1\t-1\t-1\t.\1.\t.\t-1\n";
 }
 else{
  $insert = "chr1\t1\t10\tGA\t1\t+\t.\t.\t.\t-1\t-1\t-1\t.\1.\t.\t-1\n";
 }
 while($_ = (<CHILD> || $leave && $insert)){
  if ($insert eq $_){
   $insert = 0;
  }
  chomp;
  chop($_) if ($_ =~ m/\r$/);
  $_ =~ /gene_type "([^"]*)"/;
  my $genty = $1;
  $_ =~ /transcript_type "([^"]*)"/;
  my $trty = $1;
  $_ =~ /transcript_id "([\w]*\.[\d]*)/;
  my $trid = $1;
  $_ =~ /gene_name "([^"]*)"/;
  my $genna = $1;
  $_ =~ /transcript_status "([^"]*)"/;
  my $trst = $1;
  $_ =~ /exon_number "*(\d*)"*/;
  my $accept = 1;
  if ($strict and ($_ =~ "cds_start_NF" or $_ =~ "cds_end_NF" or $_ =~ "mRNA_start_NF" or $_ =~ "mRNA_end_NF" or $trst ne "KNOWN")){
   $accept = 0;
  }
  @secs = split;
  if ($clust ne $secs[3]){
   my $annu = scalar(@anli) + 1;
   if ($clust ne ""){
    my $ann;
    foreach $key(keys(%anns)){
     if (exists(${$anns{$key}}[0])){
      if ($tend{$key} > 0){
       if ($tmatch{$key} < $length){
        next;
       }
       if ($uend{$key} > 0){
        if ($uend{$key} == $tend{$key}){
         ${$anns{$key}}[0] = "3'utr";
        }
	elsif ($uend{$key} == $eend{$key} and $estart{$key} != $tstart{$key}){
	 $line = `grep $key $fn2 | grep start_codon`;
         @line = split("\t", $line);
	 if ($line[3] > $end){
	  if ($secs[5] eq "+"){
	   ${$anns{$key}}[0] = "5'utr";
	  }
	  else{
	   ${$anns{$key}}[0] = "3'utr";
	  }
	 }
	 else{
          if ($secs[5] eq "+"){
           ${$anns{$key}}[0] = "3'utr";
          }
          else{
           ${$anns{$key}}[0] = "5'utr";
          }
	 }
	}
        else{
         ${$anns{$key}}[0] = "5'utr";
        }
        if ($umatch{$key} < $length){
         if ($ematch{$key} < $length){
          if ($start < $end){
           if ($end > $uend{$key}){
            ${$anns{$key}}[0] = ${$anns{$key}}[0] . "-intron";
           }
           else{
            ${$anns{$key}}[0] = "intron-" .  ${$anns{$key}}[0];
           }
          }
          else{
           if ($end > $uend{$key}){
            ${$anns{$key}}[0] = "intron-" .  ${$anns{$key}}[0];
           }
           else{
            ${$anns{$key}}[0] = ${$anns{$key}}[0] . "-intron";
           }
          }
         }
         else{
          if ($start < $end){
           if ($end > $uend{$key}){
            ${$anns{$key}}[0] = ${$anns{$key}}[0] . "-coding";
           }
           else{
            ${$anns{$key}}[0] = "coding-" .  ${$anns{$key}}[0];
           }
          }
          else{
           if ($end > $uend{$key}){
            ${$anns{$key}}[0] = "coding-" .  ${$anns{$key}}[0];
           }
           else{
            ${$anns{$key}}[0] = ${$anns{$key}}[0] . "-coding";
           }
          }
         }
        }
       }
       elsif ($eend{$key} > 0){
        if ($ematch{$key} < $length){
         if (exists($junction{$key})){
          if (${$anns{$key}}[0] eq "coding"){
           ${$anns{$key}}[0] = "coding-coding";
          }
          else{
           ${$anns{$key}}[0] = "lincRNA_exon-lincRNA_exon";
          }
         }
         elsif ($start < $end){
          if ($end > $eend{$key}){
           if (${$anns{$key}}[0] eq "coding"){
            ${$anns{$key}}[0] = "coding-intron";
           }
           else{
            ${$anns{$key}}[0] = "lincRNA_exon-lincRNA_intron";
           }
          }
          else{
           if (${$anns{$key}}[0] eq "coding"){
            ${$anns{$key}}[0] = "intron-coding";
           }
           else{
            ${$anns{$key}}[0] = "lincRNA_intron-lincRNA_exon"; 
           }
          }
         }
         else{
          if ($end > $eend{$key}){
           if (${$anns{$key}}[0] eq "coding"){
            ${$anns{$key}}[0] = "intron-coding";
           }
           else{
            ${$anns{$key}}[0] = "lincRNA_intron-lincRNA_exon"; 
           }
          }
          else{
           if (${$anns{$key}}[0] eq "coding"){
            ${$anns{$key}}[0] = "coding-intron";
           }
           else{
            ${$anns{$key}}[0] = "lincRNA_exon-lincRNA_intron"; 
           }
          }
         }
        }
       }
      }
      if ($verbose){
       push(@{$clan{$clust}}, ($prli . ${$anns{$key}}[0] . ${$anns{$key}}[1] . $clust . "\n", "," . ${$anns{$key}}[0] . substr(${$anns{$key}}[1], 0, -1) . "\n"));
      }
      else{
       my $index = 0;
       ++$index until $anli[$index] eq ${$anns{$key}}[0] or $index > $#anli;
       if ($index > $#anli){
        push(@anli, ${$anns{$key}}[0]);
       }
       if ($annu > $index){
        $annu = $index;
        $ann = ${$anns{$key}}[0] . ${$anns{$key}}[1];
       }
      }
     }
    }
    if ($annu < scalar(@anli) and !$verbose){
     push(@{$clan{$clust}}, ($annu, $prli . $ann . $clust . "\n", "," . substr($ann, 0, -1) . "\n"));
    }
    if (scalar(@{$clan{$clust}}) == 0){
     if ($verbose){
      push(@{$clan{$clust}}, ($prli . "NoAnnotation,,," . $clust . "\n", ",NoAnnotation,,\n"));
     }
     else{
      push(@{$clan{$clust}}, (1000, $prli . "NoAnnotation,,," . $clust . "\n", ",NoAnnotation,,\n"));
     }
    }
   }
   if ($secs[$#secs] == -1){
    next;
   }
   $clust = $secs[3];
   if ($secs[5] eq "+"){
    $start = $secs[1];
    $end = $secs[2];
    if ($bed12 and $secs[9] == 2){
     @blocks = split(",", $secs[10]);
     $end2 = $end;
     $end = $start + $blocks[0];
     $start2 = $end2 - $blocks[1] + 1;
    }
   }
   else{
    $start = $secs[2];
    $end = $secs[1];
    if ($bed12 and $secs[9] == 2){
     @blocks = split(",", $secs[10]);
     $end2 = $end;
     $end = $start - $blocks[0] + 1;
     $start2 = $end2 + $blocks[1];
    }
   }
   if ($bed12 and $secs[9] == 1){
    $start2 = 0;
    $end2 = 0;
   }
   if ($bed12){
    if ($secs[9] == 2){
     if ($secs[5] eq "+"){
      $prli = $secs[0] . "," . $secs[5] . "," . $start . "," . $end . "," . $start2 . "," . $end2 . ",";
     }
     else{
      $prli = $secs[0] . "," . $secs[5] . "," . $end2 . "," . $start2 . "," . $end . "," . $start . ",";
     }
    }
    else{
     $prli = $secs[0] . "," . $secs[5] . "," . $secs[1] . "," . $secs[2] . ",,,";
    }
   }
   else{
    $prli = $secs[0] . "," . $secs[5] . "," . $secs[1] . "," . $secs[2] . ",";
   }
   %anns = ();
   %uend = ();
   %tend = ();
   %eend = ();
   %tstart = ();
   %estart = ();
   %tmatch = ();
   %ematch = ();
   %umatch = ();
   $length = 0;
   if ($secs[$#secs] == 0){
    next;
   }
  }
  if ($secs[8 + $fix] eq "gene"){
   $length = $secs[$#secs];
   next;
  }
  if (($rseq ne "" and $trab{$trid} <= .001) or !$accept or ($bed12 and ($secs[9] > 2 or ($secs[9] == 2 and ($secs[14] ne "exon" or ($secs[5] eq "+" and $end != $secs[16] and $start2 != $secs[15]) or ($secs[5] eq "-" and $end != $secs[15] and $start2 != $secs[16])))))){
   next;
  }
  if (!exists(${$anns{$trid}}[0])){
   if ($bed12 and $secs[9] == 2){
    push(@{$anns{$trid}}, ("NoAnnotation", ",,,"));
    next;
   }
   elsif ($secs[8 + $fix] eq "transcript" or $secs[8 + $fix] eq "exon" or $secs[8 + $fix] eq "UTR" or $secs[8 + $fix] eq "CDS"){
    if ($genty eq "protein_coding" and $trty eq "protein_coding"){
     push(@{$anns{$trid}}, ("intron", "," . $genna . "," . $secs[7 + $fix] . ","));
    }
    elsif ($trty eq "lincRNA"){
     push(@{$anns{$trid}}, ("lincRNA_intron", "," . $genna . "," . $secs[7 + $fix] . ","));
    }
    else{
     push(@{$anns{$trid}}, ($trty, "," . $genna . "," . $secs[7 + $fix] . ","));
    }
   }
   else{
    push(@{$anns{$trid}}, ($secs[8 + $fix], "," . $genna . "," . $secs[7 + $fix] . ","));
   }
  }
  elsif ($bed12 and $secs[9] == 2){
   if ($trty eq "lincRNA"){
    @{$anns{$trid}} = ("lincRNA_exon-lincRNA_exon", "," . $genna . "," . $secs[13] . ",");
   }
   else{
    @{$anns{$trid}} = ("coding-coding", "," . $genna . "," . $secs[13] . ",");
   }
   next;
  }
  if ($secs[8 + $fix] eq "start_codon"){
   ${$anns{$trid}}[0] = "start_codon";
  }
  elsif ($secs[8 + $fix] eq "stop_codon"){
   ${$anns{$trid}}[0] = "stop_codon";
  }
  elsif ($secs[8 + $fix] eq "transcript" and (${$anns{$trid}}[0] eq "coding" or ${$anns{$trid}}[0] eq "intron" or ${$anns{$trid}}[0] eq "lincRNA_exon" or ${$anns{$trid}}[0] eq "lincRNA_intron")){
   if ($secs[5] eq "+"){
    $tstart{$trid} = $secs[9 + $fix];
    $tend{$trid} = $secs[10 + $fix];
   }
   else{
    $tstart{$trid} = $secs[10 + $fix];
    $tend{$trid} = $secs[9 + $fix];
   }
   $tmatch{$trid} = $secs[$#secs];
  }
  elsif ($secs[8 + $fix] eq "exon" and (${$anns{$trid}}[0] eq "intron" or ${$anns{$trid}}[0] eq "lincRNA_intron")){
   if ($secs[5] eq "+"){
    $estart{$trid} = $secs[9 + $fix];
    $eend{$trid} = $secs[10 + $fix];
   }
   else{
    $estart{$trid} = $secs[10 + $fix];
    $eend{$trid} = $secs[9 + $fix];
   }
   if (${$anns{$trid}}[0] eq "intron"){
    ${$anns{$trid}}[0] = "coding";
   }
   else{
    ${$anns{$trid}}[0] = "lincRNA_exon";
   }
   $ematch{$trid} = $secs[$#secs];
  }
  elsif ($secs[8 + $fix] eq "UTR" and (${$anns{$trid}}[0] eq "coding" or ${$anns{$trid}}[0] eq "intron")){
   if ($secs[5] eq "+"){
    $uend{$trid} = $secs[10 + $fix];
   }
   else{
    $uend{$trid} = $secs[9 + $fix];
   }
   $umatch{$trid} = $secs[$#secs];
  }
 }
 close(CHILD);
 if ($closest){
  if ($bam){
   if ($pid = open(CHILD, "-|")) {
    $bedbam = <CHILD>;
    close(CHILD);
   }
   else{
    die	"cannot fork: $!" unless defined $pid;
    STDOUT->autoflush(1);
    if ($bed12){
     open(W, "| bamToBed -i stdin -bed12");
    }
    else{
     open(W, "| bamToBed -i stdin");
    }
    print W $ctob;
    close(W);
    exit;
   }
  }
  `grep \"transcript_type \\\"protein_coding\\\"\" gencode.v18.annotation.gtf | grep \"gene_type \\\"protein_coding\\\"\" > tmp_protein_coding.gtf`;
  if ($pid = open(CHILD, "-|")) {
   my $x = 0;
  }
  else {
   die "cannot fork: $!" unless defined $pid;
   STDOUT->autoflush(1);
   open(W, "| closestBed -s -D a -io -t first -a stdin -b tmp_protein_coding.gtf");
   if ($bam){
    print W $bedbam;
   }
   else{
    print W $ctob;
   }
   close(W);
   exit;
  }
  @anns = ();
  $prli = "";
  $insert = 1;
  while($_ = (<CHILD> || $leave && $insert)){
   if ($insert eq $_){
    $insert = 0;
   }
   chomp;
   chop($_) if ($_ =~ m/\r$/);
   $_ =~ /gene_name "([^"]*)"/;
   my $genna = $1;
   $_ =~ /transcript_status "([^"]*)"/;
   my $trst = $1;
   @secs = split;
   if ($bed12 and $secs[9] > 1){
    next;
   }
   if ($strict and ($_ =~ "cds_start_NF" or $_ =~ "cds_end_NF" or $_ =~ "mRNA_start_NF" or $_ =~ "mRNA_end_NF" or $trst ne "KNOWN")){
    next;
   }
   if ($break < ${$clan{$secs[3]}}[0] and (($secs[$#secs] < 0 and $secs[$#secs] > -($closest)) or ($secs[$#secs] > 0 and $secs[$#secs] < $closest))){
    if ($bed12){
     if ($secs[9] == 2){
      if ($secs[5] eq "+"){
       $prli = $secs[0] . "," . $secs[5] . "," . $start . "," . $end . "," . $start2 . "," . $end2 . ",";
      }
      else{
       $prli = $secs[0] . "," . $secs[5] . "," . $end2 . "," . $start2 . "," . $end . "," . $start . ",";
      }
     }
     else{
      $prli = $secs[0] . "," . $secs[5] . "," . $secs[1] . "," . $secs[2] . ",,,";
     }
    }
    else{
     $prli = $secs[0] . "," . $secs[5] . "," . $secs[1] . "," . $secs[2] . ",";
    }
    if ($secs[$#secs] < 0){
     @anns = ("Downstream", "," . $genna . "," . $secs[7 + $fix] . ",");
    }
    else{
     @anns = ("Upstream", "," . $genna . "," . $secs[7 + $fix] . ",");
    }
    if ($verbose){
     push(@{$clan{$secs[3]}}, ($prli . $anns[0] . $anns[1] . $secs[3] . "\n", "," . $anns[0] . substr($anns[1], 0, -1) . "\n"));
    }
    else{
     @{$clan{$secs[3]}} = ($break, $prli . $anns[0] . $anns[1] . $secs[3] . "\n", "," . $anns[0] . substr($anns[1], 0, -1) . "\n");
    }
   }
  }
  close(CHILD);
 }
 if ($repeat){
  my %rtype;
  open(INFO, $fn4);
  while(<INFO>){
   chomp($_);
   chop($_) if ($_ =~ m/\r$/);
   @secs = split;
   if ($secs[0] ne "#repname"){
    $rtype{$secs[0]} = $secs[1];
   }
  }
  close(INFO);
  if ($pid = open(CHILD, "-|")) {
   my $x = 0;
  }
  else {
   die "cannot fork: $!" unless defined $pid;
   STDOUT->autoflush(1);
   if ($bam){
    if ($bed12){
     open(W, "| intersectBed -wo -s -split -bed -abam stdin -b $fn3");
    }
    else{
     open(W, "| intersectBed -wo -s -bed -abam stdin -b $fn3");
    }
   }
   else{
    if ($bed12){
     open(W, "| intersectBed -wo -s -split -a stdin -b $fn3");
    }
    else{
     open(W, "| intersectBed -wo -s -a stdin -b $fn3");
    }
   }
   print W $ctob;
   close(W);
   exit;
  }
  $clust = "";
  %anns = ();
  $prli = "";
  $index = 0;
  $insert = 1;
  while($_ = (<CHILD> || $leave && $insert)){
   if ($insert eq $_){
    $insert = 0;
   }
   chomp;
   chop($_) if ($_ =~ m/\r$/);
   @secs = split;
   if ($clust ne $secs[3]){
    my $annu = 1001;
    if ($clust ne ""){
     my $ann;
     $index = 0;
     foreach $key(keys(%anns)){
      if ($verbose){
       push(@{$clan{$clust}}, ($prli . ${$anns{$key}}[0] . ${$anns{$key}}[1] . $clust . "\n", "," . ${$anns{$key}}[0] . substr(${$anns{$key}}[1], 0, -1) . "\n"));
      }
      else{
       ++$index until $anli[$index] eq ${$anns{$key}}[0] or $index > $#anli;
       if ($annu > $index){
        $annu = $index;
        $ann = ${$anns{$key}}[0] . ${$anns{$key}}[1];
       }
      }
     }
     if ($annu < ${$clan{$clust}}[0] and !$verbose){
      @{$clan{$clust}} = ($annu, $prli . $ann . $clust . "\n", "," . substr($ann, 0, -1) . "\n");
     }
    }
    %anns = ();
    $clust = $secs[3];
    if ($bed12){
     if ($secs[9] == 2){
      if ($secs[5] eq "+"){
       $prli = $secs[0] . "," . $secs[5] . "," . $start . "," . $end . "," . $start2 . "," . $end2 . ",";
      }
      else{
       $prli = $secs[0] . "," . $secs[5] . "," . $end2 . "," . $start2 . "," . $end . "," . $start . ",";
      }
     }
     else{
      $prli = $secs[0] . "," . $secs[5] . "," . $secs[1] . "," . $secs[2] . ",,,";
     }
    }
    else{
     $prli = $secs[0] . "," . $secs[5] . "," . $secs[1] . "," . $secs[2] . ",";
    }
   }
   if ($bed12 and $secs[9] > 1){
    next;
   }
   $index = 0;
   ++$index until $anli[$index] eq $rtype{$secs[9 + $fix]} or $index > $#anli;
   if ($index > $#anli){
    @{$anns{$secs[9 + $fix]}} = ("repeat", "," . $secs[9 + $fix] . ",RepeatMasker,");
   }
   else{
    @{$anns{$secs[9 + $fix]}} = ($rtype{$secs[9 + $fix]}, "," . $secs[9 + $fix] . ",RepeatMasker,");
   }
  }
  close(CHILD);
 }
}
close(RF);
if (!$cluster){
 if ($verbose){
  foreach $key(sort(keys(%clan))){
   for (my $i = 0; $i <= scalar(@{$clan{$key}}) - 2; $i += 2){
    print(${$clan{$key}}[$i]);
   }
  }
 }
 else{
  foreach $key(sort(keys(%clan))){
   print(${$clan{$key}}[1]);
  }
 }
}
else{
 if ($verbose){
  foreach $key(sort(keys(%clan))){
   for (my $i = 1; $i < scalar(@{$clan{$key}}); $i += 2){
    print($prin{$key} . ${$clan{$key}}[$i]);
   }
  }
 }
 else{
  foreach $key(sort(keys(%clan))){
   print($prin{$key} . ${$clan{$key}}[2]);
  }
 }
}
