#!/usr/bin/perl

if ($#ARGV == -1){
 print("\n\t-p\tPipeline file.\n\n\t-i\tFASTQ file.\n\n\t-path\tPath for accessory files. Defaults to ./accessory/.\n\n\t-index\tBowtie index. Defaults to hg19.\n\n\t-memory\tMemory used by PARalyzer. Defaults to 8G.\n\n\t-threads\tNumber of threads used by Bowtie. Defaults to 4.\n\n\t-bow\tBowtie parameters. Defaults to -v 1 -m 10 --all --best --strata.\n\n\t-3adapt\t3' adapter sequence. Defaults to TCGTATGCCGTCTTCTGCTTG.\n\n\t-5adapt\t5' adapter sequence. Defaults to AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACGATC.\n\n\t-19mark\t19 nucleotide marker sequence. Defaults to CGTACGCGGGTTTAAACGA.\n\n\t-24mark\t24 nucleotide marker sequence. Defaults to CGTACGCGGAATAGTTTAAACTGT.\n\n\t-g\tGTF file. Defaults to gencode.v19.chr_patch_hapl_scaff.annotation.gtf.\n\n\t-bit\t2bit genome. Defaults to GRCh37.p12.genome.2bit.\n");
 die "\n";
}
my $pipe;
my $fastq;
my $path = "./accessory/";
my $index = "hg19";
my $memory = "8G";
my $threads = "4";
my $bowtie = "-v 1 -m 10 --all --best --strata";
my $tadapt = "TCGTATGCCGTCTTCTGCTTG";
my $fadapt = "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACGATC";
my $nmark = "CGTACGCGGGTTTAAACGA";
my $tmark = "CGTACGCGGAATAGTTTAAACTGT";
my $gtf = "gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz";
my $bit = "GRCh37.p12.genome.2bit";
for (my $i = 0; $i <= $#ARGV; ++$i){
 if ($ARGV[$i] eq "-p"){
  $pipe = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-i"){
  $fastq = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-path"){
  $path = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-index"){
  $index = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-memory"){
  $memory = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-threads"){
  $threads = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-bow"){
  $bowtie = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-3adapt"){
  $tadapt = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-5adapt"){
  $fadapt = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-19mark"){
  $nmark = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-24mark"){
  $tmark = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-g"){
  $gtf = $ARGV[$i + 1];
 }
 elsif ($ARGV[$i] eq "-bit"){
  $bit = $ARGV[$i + 1];
 }
}
open(P, $pipe);
@com = <P>;
close(P);
for ($i = 0; $i <= $#com; ++$i){
 if ($com[$i] =~ /^CUSTOMSCRIPTS=/){
  $com[$i] = "CUSTOMSCRIPTS=\"" . $path . "\"\n";
 }
 elsif ($com[$i] =~ /^BOWTIE_INDEX=/){
  $com[$i] = "BOWTIE_INDEX=\"\${CUSTOMSCRIPTS}" . $index . "\"\n";
 }
 elsif ($com[$i] =~ /^MEMORY_LIMIT=/){
  $com[$i] = "MEMORY_LIMIT=\"" . $memory . "\"\n";
 }
 elsif ($com[$i] =~ /^THREADS=/){
  $com[$i] = "THREADS=\"" . $threads . "\"\n";
 }
 elsif ($com[$i] =~ /^BOWTIE_PARAMS=/){
  $com[$i] = "BOWTIE_PARAMS=\"" . $bowtie . "\"\n";
 }
 elsif ($com[$i] =~ /^THREE_PRIME_ADAPTER_SEQUENCE=/){
  $com[$i] = "THREE_PRIME_ADAPTER_SEQUENCE=\"" . $tadapt . "\"\n";
 }
 elsif ($com[$i] =~ /^FIVE_PRIME_ADAPTER_SEQUENCE=/){
  $com[$i] = "FIVE_PRIME_ADAPTER_SEQUENCE=\"" . $fadapt . "\"\n";
 }
 elsif ($com[$i] =~ /^NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE=/){
  $com[$i] = "NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE=\"" . $nmark . "\"\n";
 }
 elsif ($com[$i] =~ /^TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE=/){
  $com[$i] = "TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE=\"" . $tmark . "\"\n";
 }
 elsif ($com[$i] =~ /^GTF=/){
  $com[$i] = "GTF=\"\${CUSTOMSCRIPTS}" . $gtf . "\"\n";
 }
 elsif ($com[$i] =~ /^BITFILE=/){
  $com[$i] = "BITFILE=\"\${CUSTOMSCRIPTS}" . $bit . "\"\n";
 }
 elsif ($com[$i] =~ /^@/){
  last;
 }
}
$" = "";
$com = "@com";
open(P, ">" . $pipe);
print P $com;
close(P);
`bpipe run -r $pipe $fastq`;
