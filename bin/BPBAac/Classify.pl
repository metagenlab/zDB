#!usr/bin/perl

my $cutoff = ($ARGV[0] + 1 - 1);
my $name = '';
my $svm = 0;

open(IN,"values.csv")||die"Cannot open values.csv!";
open(OUT,">FinalResult.csv")||die"Cannot open FinalResult.csv!";
print OUT "Protein", ",","SVM-Value",",","T3S protein or not","\n";
while(<IN>){
 if(/\"(\d+)\"\s\"(.*)\"\s(\S+)/){
  $name = $2;
  $svm = $3 +1 -1;
   if($svm>=$cutoff){
      print OUT $name,",", $svm, ",","YES","\n";
     }
   else{
      print OUT $name,",", $svm, ",","NO","\n";        
     }
 }
}
close(IN);
close(OUT);

system "rm *.data";
system "rm N100-M.fasta";
system "rm *.out";
system "rm values.csv";

