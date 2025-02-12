#!/usr/bin/perl -w

# Taken from Ignacio Faustino, Nov 2012

use strict;

my $long=@ARGV;
if ($long != 3){
    print "Usage: perl $0 <title> <seqLength> <outFile>\n";
    exit(0);
}

my ($title,$end,$out)=@ARGV;

#my $rexec = "/opt/R-2.15.0/bin/R";
my $rexec = "R";

print "Title: $title\n";

my $strand = ($end/2) + 0.5;

my $plot = "gnuplot";
open OUT, ">$plot.tmp";
print OUT "library(Cairo)\n";
#print OUT "par(mar=c(5,5,2,2)+0.1)\n";
print OUT "par(mar=c(4,4,4,4)+0.1)\n";
print OUT "library(plotrix)\n";
print OUT "library(MASS)\n";
#print OUT "xlabel <- c(1:$pos)\n";
print OUT "MEAN <- data.frame()\n";
print OUT "MAX <- data.frame()\n";
print OUT "MIN <- data.frame()\n";
print OUT "STDEV <- data.frame()\n";
print OUT "ZSCORE <- data.frame()\n";
print OUT "for (j in 0:$end){\n";
print OUT "\tmyfile <- paste(\"interact_\",j,sep=\"\")\n";
print OUT "\tinteract <- read.table(myfile)\n";
print OUT "\th <- j+1\n";
print OUT "\tfor (i in 1:length(interact)){\n";
print OUT "\t\tMEAN[i,h]= sapply (interact[i],mean)\n";
print OUT "\t\tMAX[i,h]= sapply (interact[i],max)\n";
print OUT "\t\tMIN[i,h]= sapply (interact[i],min)\n";
print OUT "\t\tSTDEV[i,h]= sapply (interact[i],sd)\n";
print OUT "\t}\n";
print OUT "}\n";

#print OUT "for (j in 0:$end){\n";
#print OUT "\tmyfile <- paste(\"interact_\",j,sep=\"\")\n";
#print OUT "\tinteract <- read.table(myfile)\n";
#print OUT "\th <- j+1\n";
#print OUT "\tfor (i in 1:length(interact)){\n";
#print OUT "\t\tZSCORE[i,h]= (interact[i] - MEAN[i,h]) / STDEV[i,h]\n";
#print OUT "\t}\n";
#print OUT "}\n";

# Standardizing Standard Deviation (Relative Standard Deviation, Coefficient of variation)
print OUT "for (a in 0:$end+1){ for (b in 0:$end+1){ ZSCORE[a,b] <- abs(STDEV[a,b]/MEAN[a,b]);}}\n";

# Replacing diagonal with zeros.
print OUT "for (a in 0:$end+1){ MEAN[a,a] <- 0.0; MIN[a,a] <- 0.0; MAX[a,a] <- 0.0; STDEV[a,a] <- 0.0; ZSCORE[a,a] <- 0.0;}\n";

print OUT "xdim <- dim(MEAN)\n";
print OUT "xat <- pretty(0:xdim[2])+1\n";
print OUT "yat <- pretty(0:xdim[1])+1\n";
print OUT "pxat <- pretty(0:xdim[2])+0.5\n";
print OUT "pyat <- pretty(0:xdim[1])+0.5\n";

print OUT "legendExt <- read.table('../seqR.info',colClasses='character')\n";

#print OUT "color2D.matplot(MEAN,c(1,1,0),c(0,1,0),c(0,1,1),show.legend=TRUE,\n";
#print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title Mean\",yrev=F,axes=F)\n";

print OUT "CairoPNG(file=\"${out}MEAN.dat.png\",width=700,height=700,res=72)\n";
print OUT "cellcolors <- matrix(NA,xdim[1],xdim[1])\n";
print OUT "cellcolors[MEAN<0] <- color.scale(MEAN[MEAN<0],c(1,1),c(0,1),c(0,1))\n";
print OUT "cellcolors[MEAN>=0] <- color.scale(MEAN[MEAN>=0],c(1,0),c(1,0),c(1,1))\n";
print OUT "color2D.matplot(MEAN,cellcolors=cellcolors,xlab=\"Sequence\",show.legend=F,\n";
print OUT "ylab=\"Sequence\",main=\"HB/Stacking Energy Contact Map MEAN (Kcal / mol)\",yrev=F,axes=F,show.values=F)\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "color.legend(0,-3,5,-2.5,legend=c(round(min(MEAN)),0,round(max(MEAN))),\n";
print OUT "seq(0,50),cex=0.8,rect.col=color.scale(seq(0,50),c(1,1,0),c(0,1,0),c(0,1,1)) )\n";
print OUT "write.matrix(MEAN,\"${out}MEAN.dat\")\n";
print OUT "abline(v=$strand,lwd=2,col='black')\n";
print OUT "abline(h=$strand,lwd=2,col='black')\n";

#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";
# axis(1,xat,at=0.5:1:xdim[1],las=2,labels=colnames(interact))
# axis(2,xat,at=0.5:1:xdim[2],las=2,labels=rev(colnames(interact)))
# axis(1,xat,at=0.5:1:xdim[1],las=2,labels=nops)
# axis(2,xat,at=0.5:1:xdim[1],las=2,labels=nops)

#print OUT "color2D.matplot(MAX,c(1,0),c(0,0),c(0,1),show.legend=TRUE,\n";
#print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title Max\",yrev=F,axes=F)\n";

print OUT "CairoPNG(file=\"${out}MAX.dat.png\",width=700,height=700,res=72)\n";
print OUT "cellcolors[MAX<0] <- color.scale(MAX[MAX<0],c(1,1),c(0,1),c(0,1))\n";
print OUT "cellcolors[MAX>=0] <- color.scale(MAX[MAX>=0],c(1,0),c(1,0),c(1,1))\n";
print OUT "color2D.matplot(MAX,cellcolors=cellcolors,xlab=\"Sequence\",show.legend=F,\n";
print OUT "ylab=\"Sequence\",main=\"HB/Stacking Energy Contact Map MAX (Kcal / mol)\",yrev=F,axes=F,show.values=F)\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "color.legend(0,-3,5,-2.5,legend=c(round(min(MAX)),0,round(max(MAX))),\n";
print OUT "seq(0,50),cex=0.8,rect.col=color.scale(seq(0,50),c(1,1,0),c(0,1,0),c(0,1,1)) )\n";
print OUT "write.matrix(MAX,\"${out}MAX.dat\")\n";
print OUT "abline(v=$strand,lwd=2,col='black')\n";
print OUT "abline(h=$strand,lwd=2,col='black')\n";

#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

#print OUT "color2D.matplot(MIN,c(1,0),c(0,0),c(0,1),show.legend=TRUE,\n";
#print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title Min\",yrev=F,axes=F)\n";

print OUT "CairoPNG(file=\"${out}MIN.dat.png\",width=700,height=700,res=72)\n";
print OUT "cellcolors[MIN<0] <- color.scale(MIN[MIN<0],c(1,1),c(0,1),c(0,1))\n";
print OUT "cellcolors[MIN>=0] <- color.scale(MIN[MIN>=0],c(1,0),c(1,0),c(1,1))\n";
print OUT "color2D.matplot(MIN,cellcolors=cellcolors,xlab=\"Sequence\",show.legend=F,\n";
print OUT "ylab=\"Sequence\",main=\"HB/Stacking Energy Contact Map MIN (Kcal / mol)\",yrev=F,axes=F,show.values=F)\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "color.legend(0,-3,5,-2.5,legend=c(round(min(MIN)),0,round(max(MIN))),\n";
print OUT "seq(0,50),cex=0.8,rect.col=color.scale(seq(0,50),c(1,1,0),c(0,1,0),c(0,1,1)) )\n";
print OUT "write.matrix(MIN,\"${out}MIN.dat\")\n";
print OUT "abline(v=$strand,lwd=2,col='black')\n";
print OUT "abline(h=$strand,lwd=2,col='black')\n";

#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "CairoPNG(file=\"${out}STDEV.dat.png\",width=700,height=700,res=72)\n";
print OUT "cellcolors[STDEV<0] <- color.scale(STDEV[STDEV<0],c(1,1),c(0,1),c(0,1))\n";
print OUT "cellcolors[STDEV>=0] <- color.scale(STDEV[STDEV>=0],c(1,0),c(1,0),c(1,1))\n";
print OUT "color2D.matplot(STDEV,cellcolors=cellcolors,xlab=\"Sequence\",show.legend=F,\n";
print OUT "ylab=\"Sequence\",main=\"HB/Stacking Energy Contact Map STDEV (Kcal / mol)\",yrev=F,axes=F,show.values=F)\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "color.legend(0,-3,5,-2.5,legend=c(round(min(STDEV)),0,round(max(STDEV))),\n";
print OUT "seq(0,50),cex=0.8,rect.col=color.scale(seq(0,50),c(1,1,0),c(0,1,0),c(0,1,1)) )\n";
print OUT "write.matrix(STDEV,\"${out}STDEV.dat\")\n";
print OUT "abline(v=$strand,lwd=2,col='black')\n";
print OUT "abline(h=$strand,lwd=2,col='black')\n";

#print OUT "color2D.matplot(STDEV,c(0,1),c(0,0),c(1,0),show.legend=TRUE,\n";
#print OUT "xlab=\"Sequence\",ylab=\"Sequence\",main=\"$title Stdev\",yrev=F,axes=F)\n";
#print OUT "axis(1, xat, at = pxat)\n";
#print OUT "axis(2, yat, at = pyat)\n";

print OUT "CairoPNG(file=\"${out}RSD.dat.png\",width=700,height=700,res=72)\n";
print OUT "cellcolors[ZSCORE<0] <- color.scale(ZSCORE[ZSCORE<0],c(1,1),c(0,1),c(0,1))\n";
print OUT "cellcolors[ZSCORE>=0] <- color.scale(ZSCORE[ZSCORE>=0],c(1,0),c(1,0),c(1,1))\n";
print OUT "color2D.matplot(ZSCORE,cellcolors=cellcolors,xlab=\"Sequence\",show.legend=F,\n";
print OUT "ylab=\"Sequence\",main=\"HB/Stacking Energy Contact Map RSD (Kcal / mol)\",yrev=F,axes=F,show.values=F)\n";
print OUT "axis(1, xat, at = 0.5:1:xdim[1],las=2,labels=legendExt)\n";
print OUT "axis(2, yat, at = 0.5:1:xdim[2],las=2,labels=legendExt)\n";
print OUT "color.legend(0,-3,5,-2.5,legend=c(round(min(ZSCORE)),0,round(max(ZSCORE))),\n";
print OUT "seq(0,50),cex=0.8,rect.col=color.scale(seq(0,50),c(1,1,0),c(0,1,0),c(0,1,1)) )\n";
print OUT "write.matrix(ZSCORE,\"${out}RSD.dat\")\n";
print OUT "abline(v=$strand,lwd=2,col='black')\n";
print OUT "abline(h=$strand,lwd=2,col='black')\n";

print OUT "dev.off()\n";
close OUT;

system ("$rexec < $plot.tmp --no-save");
#system ("rm $plot.tmp");

