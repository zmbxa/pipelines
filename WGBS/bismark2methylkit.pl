# BismarkCX2methylkit.pl
#NC_000932.1    3   -   0   0   CHH CAT
#NC_000932.1    4   -   0   0   CHH CCA
#NC_000932.1    5   -   0   0   CHH CCC
#NC_000932.1    6   +   0   0   CG  CGA

# chrbase chr base strnd coverage freqC freqT
open (IN,"$ARGV[0]");
open (OUT1,">${ARGV[0]}_CG_methylkit.txt");
print OUT1 "chrbase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";

while (<IN>) {
   chomp;
@a=split;
$chrbase=$a[0].".".$a[1];
$chr=$a[0];
$base=$a[1];
if ($a[2]=~/-/) { $strand= "R"; 
}
else { $strand ="F";}
$coverage=$a[3]+$a[4];
if ($coverage < 5) {next;}
$freqC= $a[3]/($a[3]+$a[4]) * 100;
$freqT= $a[4]/($a[3]+$a[4]) * 100;
$str=$chrbase."\t".$chr."\t".$base."\t".$strand."\t".$coverage."\t".$freqC."\t".$freqT."\n";
if ($a[5]=~ /CG/) { print OUT1 $str;
} 
else{print $_."\n";}
}