#!/usr/bin/perl

my $infile = $ARGV[0];
my $seq = '';
my $seq100 = '';
my $name = '';
my $line = '';
open(IN,"$infile")||die"Cannot open $infile!";
open(N100,">N100-M.fasta")||die"Cannot open 'N100-M.fasta'!";
while(<IN>){
 if(/\>(.*)/){
  $name = $1;
  if($seq){
   $seq100 = substr($seq,1,100);
   print N100 $seq100,"\n";
  }
  print N100 ">",$name,"\n";
  $seq = '';
 }
 elsif(/^(\w+)/){
  $line = $1;
  $seq = $seq.$line;
 }
}
if($seq){
  $seq100 = substr($seq,1,100);
  print N100 $seq100,"\n";
 }

close(IN);
close(N100);
################################

system "perl Feature100.pl Pos100AacFrequency N100-M.fasta >sample-Pos.out";
system "perl Feature100.pl Neg100AacFrequency N100-M.fasta >sample-Neg.out";
system "perl Integ.pl sample";

my $data = '';
open(IN1,"sample0.data")||die"Cannot open 'sample0.data'!";
open(OUT1,">sample.data")||die"Cannot open 'sample.data'!";
while(<IN1>){
 if(/(.*)\,Class/){
  $data = $1;
  print OUT1 $data,"\n";
  }
 elsif(/(.*)\,unknown/){
  $data = $1;
  print OUT1 $data,"\n";
  }
}
close(IN1);
close(OUT1);
