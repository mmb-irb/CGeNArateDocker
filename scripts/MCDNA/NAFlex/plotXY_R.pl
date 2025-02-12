#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 4){
    print "Usage: perl $0 <file.dat> <title> <X_title> <Y_title>\n";
    exit(0);
}

my ($dat,$title,$xtitle,$ytitle)=@ARGV;

#my $rexec = "/opt/R-2.15.0/bin/R";
my $rexec = "R";

print "Title: $title\n";

my $maxY = -999;
my $minY = 999;
my $maxX = -999;
my $minX = 999;
open DAT,"$dat";
while(<DAT>){
	chomp;
	# 2	0.01 
	# 3	0.06  

	my ($n,$value) = split ' ';

	$maxY = $value if ($value > $maxY);
	$minY = $value if ($value < $minY);
	$maxX = $n if ($n > $maxX);
	$minX = $n if ($n < $minX);
	
}
close DAT;

my $scaleY = ($maxY - $minY) / 10;
my $scaleX = int(($maxX - $minX) / 10);
print "MaxY: $maxY, Min: $minY, Scale: $scaleY\n";
print "MaxX: $maxX, Min: $minX, Scale: $scaleX\n";
$maxX += $scaleX;
$minX -= $scaleX;
$maxY += $scaleY;
$minY -= $scaleY;
print "NewMaxY: $maxY, NewMin: $minY\n";
print "NewMaxX: $maxX, NewMin: $minX\n";

my $plot = "gnuplot";
open OUT, ">$plot.tmp";
print OUT "library(Cairo)\n";
print OUT "CairoPNG(file=\"$dat.png\",width=1000,height=400)\n";
print OUT "layout(matrix(c(1,2),ncol=2), widths=c(6,2), heights=c(1,1))\n";
print OUT "par(cex=1,xpd=NA)\n";
print OUT "param <- read.table(\"$dat\",header=F)\n";
print OUT "prm <- param[[1]]\n";
print OUT "prm2 <- param[[2]]\n";
print OUT "plot(prm,prm2,ylab=\"$ytitle\",xlab=\"$xtitle\",col=\"red\",type=\"l\",lty=1,axes=F)\n";
print OUT "axis(2)\n";
print OUT "axis(1)\n";
print OUT "title(main=\"$title\")\n";
#print OUT "legend(10,\"$ytitle\",x=0,xjust=0,yjust=0,col=\"red\",lty=1,lwd=2)\n";
print OUT "dens <- density(prm2)\n";
print OUT "population <- dens\$y * diff(dens\$x)\n"; 		# Converting density to %1
#print OUT "sum(population)\n";					# Sum should be 1.
#print OUT "population <- (dens\$y * diff(dens\$x))*100\n"; 	# If %100 wanted
#print OUT "write.csv(dens\$x, file=\"$dat.density.csv\")\n";
#print OUT "write.csv(population, file=\"$dat.population.csv\")\n";
print OUT "write.table( cbind(dens\$x,population),\"$dat.population.csv\", sep=\",\", col.names=FALSE)\n";
print OUT "plot(population,dens\$x,axes=F,col=\"red\",ylab=\"\",xlab=\"density\")\n";
print OUT "axis(1)\n";
print OUT "axis(4)\n";
print OUT "dev.off()\n";
close OUT;

system ("$rexec < $plot.tmp --no-save");
#system ("rm $plot.tmp");

