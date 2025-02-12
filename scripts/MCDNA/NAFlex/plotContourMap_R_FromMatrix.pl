#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 4){
    print "Usage: perl $0 <inputMatrix> <title> <seqLength> <outFile>\n";
    exit(0);
}

my ($mat,$title,$end,$out)=@ARGV;

my $rexec = "/opt/R-2.15.0/bin/R";

print "Title: $title\n";


my $plot = "gnuplot";
open OUT, ">$plot.tmp";
print OUT "library(Cairo)\n";
print OUT "m <- as.matrix(read.table(\"$mat\"))\n";

print OUT "xdim <- dim(m)\n";
print OUT "xat <- pretty(0:xdim[2])+1\n";
print OUT "yat <- pretty(0:xdim[1])+1\n";
print OUT "pxat <- pretty(0:xdim[2])+0.5\n";
print OUT "pyat <- pretty(0:xdim[1])+0.5\n";

print OUT "CairoPNG(file=\"${out}.png\",width=700,height=700,res=72)\n";

#print OUT "color2D.matplot(MEAN,c(1,0),c(0,0),c(0,1),show.legend=TRUE,\n";
#print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title Mean\",yrev=F,axes=F)\n";
#print OUT "filled.contour(m,xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title\",color.palette = colorRampPalette(c('blue','green', 'yellow','red'), space='rgb'))\n";

print OUT "x <- seq(1,$end,length.out=ncol(m))\n";
print OUT "y <- seq(1,$end,length.out=ncol(m))\n";

print OUT "filled.contour(x,y,m,xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title\",color.palette = colorRampPalette(c('gray','yellow','red'), space='rgb'))\n";

#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "dev.off()\n";
close OUT;

system ("$rexec < $plot.tmp --no-save");
#system ("rm $plot.tmp");
system ("cp $plot.tmp patata.tmp");
