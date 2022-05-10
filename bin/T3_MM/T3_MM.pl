#!/usr/bin/perl
my $infile1 = "log.freq.ratio.txt";
my $infile2 = $ARGV[0];
my $seq = '';
my $seqLength = 0;
my $bi_aa = '';
my $count = 0;
my $sum_log_ratio = 0;
my %log_ratio = ();
my $init = '';
my $flag1 = 0;
my $flag2 = 0;
my $seqName2 = '';

open(IN1,"$infile1")||die"Cannot open $infile1!";
while(<IN1>){
 if(/(\S+)\t(\S+)/){
   $bi_aa = $1;
   $log_ratio{$bi_aa}=$2+1-1;
 }
}
close(IN1);

open(IN_2,"$infile2")||die"Cannot open $infile2!";
open(OUT_1,">tmp.fasta")||die"Cannot open the file!";
open(OUT1,">sum_log_ratio1.txt")||die"Cannot open sum_log_ratio1.txt!";
open(OUT2,">sum_log_ratio2.txt")||die"Cannot open sum_log_ratio2.txt!";
open(OUT3,">1.R")||die"Cannot open 1.R";
open(OUT4,">2.R")||die"Cannot open 2.R";
print OUT1 seqName,",R-value,T3SE\n";
print OUT2 seqName,",R-value,T3SE\n";
print OUT3 "T3SE<-c(";
print OUT4 "nonT3SE<-c(";
while(<IN_2>){
 if(/\>(\S+)/){
  $seqName = $1;
  if($seqName2){
   print OUT_1 ">",$seqName2,"\n";
   print OUT_1 $seq,"\n";
  }
  $seq = '';
  $seqName2 = $seqName;
 }
 elsif(/(\S+)/){
   $seq = $seq.$1;
 }
}
  if($seqName2){
   print OUT_1 ">",$seqName2,"\n";
   print OUT_1 $seq,"\n";
  }
close(IN_2);
close(OUT_1);

open(IN2,"tmp.fasta")||die"Cannot open tmp.fasta!";
while(<IN2>){
 if(/\>(\S+)/){
  $seqName = $1;
  $count++;
 }
 elsif(/(\S+)/){
   $seq = $1;
   if(length($seq)<=101){
     $seqLength = length($seq)-1;
   }
   else{
     $seqLength = 101;
   }
   $init = substr($seq,1,1);
   $sum_log_ratio = $log_ratio{$init};
   for(my $i=1;$i<($seqLength - 1);$i++){
      $bi_aa = substr($seq,$i,2);
      $sum_log_ratio = $sum_log_ratio + $log_ratio{$bi_aa};
   }
   $sum_log_ratio = ($sum_log_ratio * 100)/$seqLength;
   $sum_log_ratio = $sum_log_ratio/30;
   if($sum_log_ratio>=0){print OUT1 $seqName,",",$sum_log_ratio,",","YES","\n";}
   else{print OUT2 $seqName,",",$sum_log_ratio,",","NO","\n";}
   if($sum_log_ratio>=0){
    if($flag1 == 0){
     print OUT3 $sum_log_ratio;
     $flag1 = 1;	
    }
    else{
     print OUT3 ",",$sum_log_ratio;
    }	
   }
   if($sum_log_ratio<0){
    if($flag2 == 0){
     print OUT4 $sum_log_ratio;
     $flag2 = 1;	
    }
    else{
     print OUT4 ",",$sum_log_ratio;
    }	
   }
 }
}
print OUT3 ");\nx<-(pnorm((T3SE-0.28)/0.26)+pnorm((T3SE+0.28)/0.22))/2;\nwrite.table(x,file=\"t3se.txt\");\n";
print OUT4 ");\ny<-1-(pnorm((nonT3SE-0.28)/0.26)+pnorm((nonT3SE+0.28)/0.22))/2;\nwrite.table(y,file=\"nont3se.txt\");\n";
close(IN2);
close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);

system "Rscript 1.R";
system "Rscript 2.R";

my $count = 0;
my $pvalue = 0;
my $seq = '';
my %pvalue = ();
open(IN3,"t3se.txt")||die"Cannot open t3se.txt";
while(<IN3>){
 if(/^(\S+)\s\s*(\S+)/){
 	  $pvalue = $2+1-1;
 	  $count++;
 	  $seq = 'seq'.$count;
 	  $pvalue{$seq} = $pvalue;
 	}	
}
close(IN3);

 $count = -1;
open(IN4,"sum_log_ratio1.txt")||die"Cannot open sum_log_ratio1.txt";
open(OUT5,">final.result.csv")||die"Cannot open final.result.csv";
print OUT5 "seqName,value,T3SE,probability\n";
while(<IN4>){
 if(/(\S+)/){
 	  $ann = $1;
 	  $count++;
 	  $seq = 'seq'.$count;
 	  if(defined $pvalue{$seq}){
 	   print OUT5 $ann,",",$pvalue{$seq},"\n";
 	  }
 	}	
}
close(IN4);

my $count = 0;
open(IN5,"nont3se.txt")||die"Cannot open nont3se.txt";
while(<IN5>){
 if(/^(\S+)\s\s*(\S+)/){
 	  $pvalue = $2+1-1;
 	  $count++;
 	  $seq = 'seq-'.$count;
 	  $pvalue{$seq} = $pvalue;
 	}	
}
close(IN5);

 $count = -1;
open(IN6,"sum_log_ratio2.txt")||die"Cannot open sum_log_ratio2.txt";
while(<IN6>){
 if(/(\S+)/){
 	  $ann = $1;
 	  $count++;
 	  $seq = 'seq-'.$count;
 	  if(defined $pvalue{$seq}){
 	   print OUT5 $ann,",",$pvalue{$seq},"\n";
 	  }
 	}	
}
close(IN6);
close(OUT5);
system "rm 1.R";
system "rm 2.R";
system "rm t3se.txt";
system "rm nont3se.txt";
system "rm sum_log_ratio1.txt";
system "rm sum_log_ratio2.txt";
system "rm tmp.fasta";
