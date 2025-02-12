#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long < 2 or $long > 5){
    print "Usage: perl $0 <prefixInputFiles> <outFile> [<title>] [<yaxis>] [<xaxis>]\n";
    print "Example: perl $0 dat2plot_1607174998.dat patata.HB 'Hydrogen Bonds Plot' 'Distances (Angstrom)' 'Atom Pairs'\n";
    exit(0);
}

my ($inFile,$outFile,$mainTitle,$ylab,$xlab)=@ARGV;

my $rexec = "/opt/R-2.15.0/bin/R";

$xlab = "" if (!$xlab);
$ylab = "(%)" if(!$ylab);
$mainTitle = $outFile if(!$mainTitle);

open RIN,">tempR.in";

print RIN "library(Cairo)\n";
print RIN "CairoPNG(file=\"$outFile.png\",width=1000,height=600)\n";
#print RIN "par(cex=1,xpd=NA,mgp = c(4,1,0),mar = c(7.0,5.1,4.1,2.1))\n";
print RIN "n <- read.table(\"$inFile\")\n";
print RIN "names(n) <- c('North','East','South','West')\n";
print RIN "n1 <- as.matrix(n)\n";
print RIN "barplot(n1,ylab=\"$ylab\",xlab=\"$xlab\",main=\"$mainTitle\")\n";
print RIN "dev.off()\n";
close RIN;

system ("$rexec < tempR.in --no-save");

