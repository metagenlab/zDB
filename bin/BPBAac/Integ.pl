
my $file = $ARGV[0];
my $infile1 = $file.'-Pos.out';
my $infile2 = $file.'-Neg.out';
my $outfile = "sample0.data";

my $newseq = '';
my $name = '';
my $seq = '';
my %seq = ();
my $count = 0;


open(IN1,"$infile1")||die"Cannot open $infile1!";
while(<IN1>){
	$seq = $_;
	chomp $seq;
	$count++;
	$name = 'name'.$count;
  $seq{$name} = $seq;
}
close(IN1);

open(IN2,"$infile2")||die"Cannot open $infile2!";
open(OUT,">$outfile")||die"Cannot open $outfile!";
$count = 0;
print OUT "Name,paa1,paa2,paa3,paa4,paa5,paa6,paa7,paa8,paa9,paa10,paa11,paa12,paa13,paa14,paa15,paa16,paa17,paa18,paa19,paa20,paa21,paa22,paa23,paa24,paa25,paa26,paa27,paa28,paa29,paa30,paa31,paa32,paa33,paa34,paa35,paa36,paa37,paa38,paa39,paa40,paa41,paa42,paa43,paa44,paa45,paa46,paa47,paa48,paa49,paa50,paa51,paa52,paa53,paa54,paa55,paa56,paa57,paa58,paa59,paa60,paa61,paa62,paa63,paa64,paa65,paa66,paa67,paa68,paa69,paa70,paa71,paa72,paa73,paa74,paa75,paa76,paa77,paa78,paa79,paa80,paa81,paa82,paa83,paa84,paa85,paa86,paa87,paa88,paa89,paa90,paa91,paa92,paa93,paa94,paa95,paa96,paa97,paa98,paa99,paa100,Name0,naa1,naa2,naa3,naa4,naa5,naa6,naa7,naa8,naa9,naa10,naa11,naa12,naa13,naa14,naa15,naa16,naa17,naa18,naa19,naa20,naa21,naa22,naa23,naa24,naa25,naa26,naa27,naa28,naa29,naa30,naa31,naa32,naa33,naa34,naa35,naa36,naa37,naa38,naa39,naa40,naa41,naa42,naa43,naa44,naa45,naa46,naa47,naa48,naa49,naa50,naa51,naa52,naa53,naa54,naa55,naa56,naa57,naa58,naa59,naa60,naa61,naa62,naa63,naa64,naa65,naa66,naa67,naa68,naa69,naa70,naa71,naa72,naa73,naa74,naa75,naa76,naa77,naa78,naa79,naa80,naa81,naa82,naa83,naa84,naa85,naa86,naa87,naa88,naa89,naa90,naa91,naa92,naa93,naa94,naa95,naa96,naa97,naa98,naa99,naa100,Class\n";
while(<IN2>){
	$seq = $_;
	chomp $seq;
	$count++;
	$name = 'name'.$count;
  if(defined $seq{$name}){
   $newseq = $seq{$name}.$seq.'unknown';
   print OUT $newseq,"\n";
  }
}
close(IN2);
close(OUT);
