#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long < 3 or $long > 6){
    print "Usage: perl $0 <alphaSerFile> <gammaSerFile> <outFile> [<title>] [<xaxis>] [<yaxis>]\n";
    print "Example: perl $0 NAFlex_DDDbsc1_canalOut_alphaW.ser NAFlex_DDDbsc1_canalOut_gammaW.ser patata.AlphaGamma 'Alpha-Gamma Plot' 'Alpha' 'Gamma'\n";
    exit(0);
}

my ($alpha,$gamma,$outFile,$mainTitle,$xlab,$ylab)=@ARGV;

# WARNING: module load R needed to run this script.
my $rexec = "R";

$xlab = "Property X" if (!$xlab);
$ylab = "Property Y" if(!$ylab);
$mainTitle = $outFile if(!$mainTitle);

my $lims;
$lims = "c(0,360)" if ($xlab =~ /amma/ or $ylab =~ /amma/);
$lims = "c(100,350)" if ($xlab !~ /amma/ and $ylab !~ /amma/);

# Pre-processing files
my %hash;
my ($snp,$zero,$len,@arr);
open ALPHA,"$alpha";
while(<ALPHA>){
	chomp;
        # 1    0.00  -85.47  -85.67  -74.24  -59.28  -68.70  -72.18  -66.51  -62.87  -63.10  -72.56  -69.55
	($snp,@arr) = split ' ' if($alpha !~ /alpha/);
	($snp,$zero,@arr) = split ' ' if($alpha =~ /alpha/);
	$len = $#arr;
	$len-- if ($gamma =~ /zeta/);
	for (my $i=0;$i<=$len;$i++) {
		my $col = $i+1;
		push @{$hash{$col}},$arr[$i];
	}
}
close ALPHA;

foreach my $f (sort keys %hash){
	open OUT,">plotAlpha$f.tmp.dat";
	foreach my $v (@{$hash{$f}}){
		$v = ($v<=0) ? $v+360 : $v;
		print OUT "$v\n";
	}
	close OUT;
}

undef(%hash);
open GAMMA,"$gamma";
while(<GAMMA>){
	chomp;
        # 1    0.00  -85.47  -85.67  -74.24  -59.28  -68.70  -72.18  -66.51  -62.87  -63.10  -72.56  -69.55
	($snp,@arr) = split ' ' if($alpha !~ /alpha/);
	($snp,$zero,@arr) = split ' ' if($alpha =~ /alpha/);
	$len = $#arr;
	$len-- if ($gamma =~ /zeta/);
	for (my $i=0;$i<=$len;$i++) {
		my $col = $i+1;
		push @{$hash{$col}},$arr[$i];
	}
}
close GAMMA;

foreach my $f (sort keys %hash){
	open OUT,">plotGamma$f.tmp.dat";
	foreach my $v (@{$hash{$f}}){
		$v = ($v<=0) ? $v+360 : $v;
		print OUT "$v\n";
	}
	close OUT;
}

`cat plotAlpha*.tmp.dat > plotAlphaALL.tmp`;
`cat plotGamma*.tmp.dat > plotGammaALL.tmp`;
`paste plotAlphaALL.tmp plotGammaALL.tmp > AlphaGamma.dat` if($alpha =~ /alpha/);
`paste plotAlphaALL.tmp plotGammaALL.tmp > EpsilonZeta.dat` if($alpha !~ /alpha/);

open RIN,">tempR.in";

print RIN "library(Cairo)\n";
print RIN "library(ggplot2)\n";
print RIN "library(grid)\n";
print RIN "CairoPNG(file=\"$outFile.png\",width=650,height=650)\n";

#print RIN "data <- read.table(\"$alpha\")\n";
#print RIN "data2 <- read.table(\"$gamma\")\n\n";

my @list;
foreach my $f (sort keys %hash){
	print RIN "a$f <- read.table(\"plotAlpha$f.tmp.dat\")\n"; 
	print RIN "g$f <- read.table(\"plotGamma$f.tmp.dat\")\n"; 
	print RIN "d$f <- data.frame(a$f,g$f)\n";
	print RIN "names(d$f)<-c(\"alpha\",\"gamma\")\n";
	push @list,"d$f";
}

my $list = join ',',@list;
print RIN "all <- rbind($list)\n";

print RIN "p <- ggplot()+\n";
print RIN "  geom_point(size=1.8,alpha=0.35,data=all,aes(x=alpha,y=gamma,col=\"Zeta\"),colour=\"red\")+\n";
print RIN "  xlab(\"\\n$xlab (degrees)\")+ylab(\"$ylab (degrees)\\n\")+\n";
print RIN "  scale_x_continuous(expand=c(0,0),limits=$lims)+\n";
print RIN "  scale_y_continuous(expand=c(0,0),limits=$lims)+\n";
print RIN "  coord_cartesian(ylim=$lims)+\n";
print RIN "  theme_grey()+\n";
#print RIN "  theme(plot.margin=unit(c(-1,-1,-1,-1),\"mm\"))+ \n";
print RIN "  theme(axis.text=element_text(size=24),legend.position=\"none\",axis.title=element_text(size=20)) \n";
#print RIN "  labs(x=NULL, y=NULL)\n";

print RIN "p\n";

close RIN;

system ("$rexec < tempR.in --no-save");

`rm *.tmp.dat`;
`rm *.tmp`;
