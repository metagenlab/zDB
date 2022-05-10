#!/usr/bin/perl
#perl Feature100.pl FREQFILENAME SEQFILENAME >OUTFILE

my $freqfile = $ARGV[0];
my $seqfile = $ARGV[1];
my $pi = '';
my %fre = ();
my $name = '';
my $additive = 0;
my @seq = [];

open(INDEX,"$freqfile")||die"Cannot open $freqfile!";
while(<INDEX>){
 if(/(\S+)\t(\S+)/){
  $pi = $1;
  $fre{$pi} = $2;
 }
}
close(INDEX);

open(SEQ,"$seqfile")||die"Cannot open $seqfile!";
while(<SEQ>){
 if(/\>(\w+)/){
  $name = $1;
  print "$name,";
 }
 else{
  $line = $_;
  chomp($line);
   
  @seq = split(//,$line);
  	for(my $i = 0;$i < 100;$i++){
		 $aa = $seq[$i];
		 $additive = ($i + 1);
		 $pi = "$aa"."$additive";
		 if(!(defined $fre{$pi})){
		 	$fre{$pi} = 0;
		 	}
		 print $fre{$pi},",";
		}
		print "\n";
 }
}
close(SEQ);
