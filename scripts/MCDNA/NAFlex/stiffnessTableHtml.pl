#!/usr/bin/perl -w
use strict;

my $long=@ARGV;
if ($long != 2){
    print "Usage: perl $0 <input_stiffness.dat> <output_stiffness.html>\n";
    exit(0);
}

my ($file,$out) = @ARGV;

my $title = $file;
$title =~s/\.cte//g;
$title = uc($title);

my @hp = ('Shift','Slide','Rise','Tilt','Roll','Twist');

my $code = "<h1>$title</h1>\n";
$code .= "<table cellpadding='15' align='left' border='0'>\n";
$code .="<tr align='center' style='background-color:#dcdcdc'><td style='background-color:#ffffff'></td><td>Shift</td><td>Slide</td><td>Rise</td><td>Tilt</td><td>Roll</td><td>Twist</td></tr>\n";

# Getting Min & Max Values
my $max = -999;
my $min = 999;

open TABLE,"$file";
while (<TABLE>){ 
	my @array  = split ' ';
	foreach my $value (@array){
		if($value){
			$max = $value if($value > $max);
			$min = $value if($value < $min);
                }
        }
}
close TABLE;

# Writting Table
my $i = 0;
open TABLE,"$file";
while (<TABLE>){       
	my @array  = split ' ';
	$code .= "<tr>\n";
	$code .= "<td style='background-color:#dcdcdc'>$hp[$i]</td>\n";
	foreach my $value (@array){
		if($value){
			$value = sprintf("%8.3f",$value);
			my $rgb = tableColor($value,$min);
			$code .= "<td style='background-color: rgb($rgb)'>$value</td>\n";
		}
	}
	$i++;
	$code .= "</tr>\n";
}
$code .= "</table>\n";
close(TABLE);

$code .= "<br style='clear: both;'/>\n<p align='left' style='float: left';><b>Units:<br/> </b><i>Diagonal Shift/Slide/Rise in kcal/(mol*&Aring;&sup2;), Diagonal Tilt/Roll/Twist in kcal/(mol*degree&sup2;)<br/>\n";
$code .= "Out of Diagonal: Shift/Slide/Rise in kcal/(mol*&Aring;), Out of Diagonal Tilt/Roll/Twist in kcal/(mol*degree)</i></p>\n";


open OUT,">$out";
print OUT "$code";
close OUT;

#####################################

sub tableColor {
	my ($number, $min) = @_;

	my $num = ($number - $min) / 10;
	my $percent = 255 * $num;

	my $r = 255;
	my $g = sprintf("%d",240 - $percent);
	my $b = sprintf("%d",240 - $percent);

	return "$r,$g,$b";
}

