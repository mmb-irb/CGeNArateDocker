#!/usr/bin/perl -w

use strict;

my $long=@ARGV;
if ($long != 5){
    print "Usage: perl $0 <title> <ini> <end> <offset> <outFile>\n";
    exit(0);
}

my ($title,$ini,$end,$offset,$out)=@ARGV;

#my $rexec = "/opt/R-2.15.0/bin/R";
my $rexec = "R";

print "Title: $title\n";

# Matrix number of cols/rows
my $pos = int(($end - $ini) / $offset) + 1;
my $endR = $end - 1;
my $strand = $end/2;

my $plot = "gnuplot";
open OUT, ">$plot.tmp";
print OUT "library(Cairo)\n";
print OUT "par(mar=c(5,5,2,2)+0.1)\n";
print OUT "library(plotrix)\n";
print OUT "library(MASS)\n";
#print OUT "xlabel <- c(1:$pos)\n";
print OUT "MEAN <- matrix(nrow = $pos, ncol = $pos)\n";
print OUT "STDEV <- matrix(nrow = $pos, ncol = $pos)\n";
print OUT "MIN <- matrix(nrow = $pos, ncol = $pos)\n";
print OUT "MAX <- matrix(nrow = $pos, ncol = $pos)\n";
print OUT "p <- $pos-1\n";
print OUT "pp <- $pos\n";
print OUT "conti <- $ini\n";
print OUT "for (i in 1:p){\n";
print OUT "\tMEAN[i,i]= 0.0\n";
print OUT "\tSTDEV[i,i]= 0.0\n";
print OUT "\tMIN[i,i]= 0.0\n";
print OUT "\tMAX[i,i]= 0.0\n";
print OUT "\tpatata <- i+1\n";
print OUT "\tcontj <- $ini + (i * $offset)\n";
print OUT "\tfor (j in patata:pp){\n";
#print OUT "\t\tindex1 <- i + $conti\n";
#print OUT "\t\tindex2 <- j + $contj\n";
print OUT "\t\tmyfile <- paste(conti,\"-\",contj,\".dat\",sep=\"\")\n";
print OUT "\t\tcontj <- contj + $offset\n";
print OUT "\tif (file.exists(myfile)){\n";
print OUT "\t\tinteract <- read.table(myfile)\n";
print OUT "\t\tMEAN[i,j]= sapply (interact[2],mean)\n";
print OUT "\t\tMEAN[j,i]= sapply (interact[2],mean)\n";
print OUT "\t\tSTDEV[i,j]= sapply (interact[2],sd)\n";
print OUT "\t\tSTDEV[j,i]= sapply (interact[2],sd)\n";
print OUT "\t\tMIN[i,j]= sapply (interact[2],min)\n";
print OUT "\t\tMIN[j,i]= sapply (interact[2],min)\n";
print OUT "\t\tMAX[i,j]= sapply (interact[2],max)\n";
print OUT "\t\tMAX[j,i]= sapply (interact[2],max)\n";
print OUT "\t}\n";
#print OUT "\t\tMEAN[i,j]=mean(interact[2])\n";
#print OUT "\t\tMEAN[j,i]=mean(interact[2])\n";
#print OUT "\t\tSTDEV[i,j]=sd(interact[2])\n";
#print OUT "\t\tSTDEV[j,i]=sd(interact[2])\n";
#print OUT "\t\tMIN[i,j]=min(interact[2])\n";
#print OUT "\t\tMIN[j,i]=min(interact[2])\n";
#print OUT "\t\tMAX[i,j]=max(interact[2])\n";
#print OUT "\t\tMAX[j,i]=max(interact[2])\n";
print OUT "\t}\n";
print OUT "\tconti <- conti + $offset\n";
print OUT "}\n";
print OUT "MEAN[j,j]= 0.0\n";
print OUT "STDEV[j,j]= 0.0\n";
print OUT "MIN[j,j]= 0.0\n";
print OUT "MAX[j,j]= 0.0\n";

print OUT "xdim <- dim(MEAN)\n";
print OUT "xat <- pretty(0:xdim[2]) * $offset + $ini\n";
print OUT "yat <- pretty(0:xdim[1]) * $offset + $ini\n";
print OUT "pxat <- pretty(0:xdim[2])+0.5\n";
print OUT "pyat <- pretty(0:xdim[1])+0.5\n";

print OUT "legendExt <- read.table('seqR.info',colClasses='character')\n";

print OUT "CairoPNG(file=\"${out}MEAN.dat.png\",width=700,height=700,res=72)\n";
print OUT "color2D.matplot(MEAN,c(1,1),c(0,1),c(0,1),show.legend=TRUE,\n";
print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title MEAN\",yrev=F,axes=F)\n";
print OUT "write.matrix(MEAN,\"${out}MEAN.dat\")\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "abline(v=$strand,lwd=2,col='blue')\n";
print OUT "abline(h=$strand,lwd=2,col='blue')\n";
#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "CairoPNG(file=\"${out}MIN.dat.png\",width=700,height=700,res=72)\n";
print OUT "color2D.matplot(MIN,c(1,1),c(0,1),c(0,1),show.legend=TRUE,\n";
print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title MIN\",yrev=F,axes=F)\n";
print OUT "write.matrix(MIN,\"${out}MIN.dat\")\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "abline(v=$strand,lwd=2,col='blue')\n";
print OUT "abline(h=$strand,lwd=2,col='blue')\n";
#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "CairoPNG(file=\"${out}MAX.dat.png\",width=700,height=700,res=72)\n";
print OUT "color2D.matplot(MAX,c(1,1),c(0,1),c(0,1),show.legend=TRUE,\n";
print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title MAX\",yrev=F,axes=F)\n";
print OUT "write.matrix(MAX,\"${out}MAX.dat\")\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "abline(v=$strand,lwd=2,col='blue')\n";
print OUT "abline(h=$strand,lwd=2,col='blue')\n";
#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "CairoPNG(file=\"${out}STDEV.dat.png\",width=700,height=700,res=72)\n";
print OUT "color2D.matplot(STDEV,c(1,1),c(0,1),c(0,1),show.legend=TRUE,\n";
print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title STDEV\",yrev=F,axes=F)\n";
print OUT "write.matrix(STDEV,\"${out}STDEV.dat\")\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "abline(v=$strand,lwd=2,col='blue')\n";
print OUT "abline(h=$strand,lwd=2,col='blue')\n";
#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "dev.off()\n";
close OUT;

system ("$rexec < $plot.tmp --no-save");
#system ("rm $plot.tmp");

