#!/usr/bin/perl -w

# aim :  counting methylated/unmethylated cytosines within genes
# Usage : perl counting_mC_within_genes.pl INPUT1 INPUT2

# INPUT1 : methylation data
# INPUT1 : Chr position ori methylC totalC context
# INPUT1 : context should be X (CG), Y(CHG), Z(CHH)
# INPUT1 : the file should be sorted by Chr and position
# example of INPUT1 : Chr1 109 + 17 19 X
# INPUT2 : tair10_gene_model_startsort.txt
# INPUT2 : start end ID chr ori genetype
# INPUT2 : gene type can be 0 (gene), 1 (te gene), 2 (other)
# INPUT2 : the file should be sorted by Chr and position
# example of INPUT2 : 3631 5899 AT1G01010 1 + 0


my $outfile;
if ($ARGV[0] =~ /(.+)(\.txt)$/){$outfile = $1;} else {die "txt file required: $!"};
$outfile = $outfile.".gene10.txt";
print STDERR "INPUT : @ARGV \n";

my (@cn, @c_p, @c_m, @c_t, @context); 
@cn=(undef,0,0,0,0,0);  #cytosine number initial value
my ($p, $ch, $chro)= (0, 0, 0); 


#reading methylation data (INPUT1)
open (INPUT1,$ARGV[0]) or die "can not :$!";
while (<INPUT1>) {    												   
  if(/^Chr(\d)\s(\d+)\s(.)\s(\d+)\s(\d+)\s(.)/) {
      $chro = $1;
	  $p = ++ $cn[$1];
	  $c_p[$chro][$p] = $2; #cytosine position
	  $c_m[$chro][$p] = $4;  #cytosine methylated
	  $c_t[$chro][$p] = $5; #cytosine total
      $context[$chro][$p] = $6; #context
	  }}
print STDERR "cytosine read done :cytosine num @cn[1,2,3,4,5] \n";
close (INPUT1);

my (@gn, @gs, @ge, @gene, @ori, @te, @lg, $chr, $g);
@gn =(undef,0,0,0,0,0);

#reading gene data (INPUT2)
open (INPUT2,$ARGV[1]) or die "can not :$!";
while (<INPUT2>) {   
  if(/^(\d+)\s(\d+)\s(AT.G.....)\s(\d)\s(.)\s(\d)/){     
    $chr = $4;
	$g = ++$gn[$chr];
    $gs[$chr][$g] = $1;  #start of gene
	$ge[$chr][$g] = $2;  #end of gene
    $gene[$chr][$g] = $3; # gene ID
	$ori[$chr][$g] = $5; #orienation
	$te[$chr][$g] =$6;   #TE or not
	$lg[$chr][$g] = $ge[$chr][$g] - $gs[$chr][$g]; # length of gene;
	}}
print STDERR print "gene read done: @gn[1,2,3,4,5] \n";
close (INPUT2);

#calculate
my ($c, $i, $j) = (1,1,1); 
my ($gtX, $gcX, $gtY, $gcY, $gtZ, $gcZ);

open (OUT, ">", $outfile) or die "can not :$!";
print OUT "gene mX tX mY tY mZ tZ ori te length\n";#hedder

for ($c = 1; $c <= 5; $c++) { 
print STDERR "now starting chromosome $c\n"; 
($i, $j) = (1,1);

while ($c_p[$c][$j] < $gs[$c][$i]) {$j++;} #go forward if the current C position is smaller than the start of the first gene  
 
while ($i <= $gn[$c]) {
   ($gtX, $gcX, $gtY, $gcY, $gtZ, $gcZ) = (0,0,0,0,0,0); #initialize the variables before counting C within the gene

   while ($c_p[$c][$j] <= $ge[$c][$i]) {      # counting mC and tC while C position is within gene
	 if ($context[$c][$j] eq "X")    {$gtX += $c_t[$c][$j]; $gcX += $c_m[$c][$j];}
     elsif ($context[$c][$j] eq "Y") {$gtY += $c_t[$c][$j]; $gcY += $c_m[$c][$j];} 
	 elsif ($context[$c][$j] eq "Z") {$gtZ += $c_t[$c][$j]; $gcZ += $c_m[$c][$j];}
	 $j++; if ($cn[$c] < $j) {last;}
        						   }  #j loop
 print OUT "$gene[$c][$i] $gcX $gtX $gcY $gtY $gcZ $gtZ $ori[$c][$i] $te[$c][$i] $lg[$c][$i]\n"; # output
   $i++; #next gene
   if ($gn[$c] < $i) {last;}
   while ($gs[$c][$i] < $c_p[$c][$j]) {$j--; if ($j == 1) {last;}}  #go back if the current C position is bigger than the start of next gene
   while ($c_p[$c][$j] < $gs[$c][$i]) {$j++; if ($cn[$c] < $j) {last;}} # go forward if the current C position is smaller than the start of next gene  
                  } #i loop
} #c loop

close (OUT);
__END__