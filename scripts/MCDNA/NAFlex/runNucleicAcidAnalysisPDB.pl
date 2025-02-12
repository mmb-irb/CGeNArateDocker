#!/usr/bin/perl 

use lib "/app/Scripts/MCDNA/NAFlex";
use NA_MDWeb_Common;
use Switch;
use strict;

require "NA_MDWeb_curvesPDB.pl";
require "NA_MDWeb_stiffnessPDB.pl";
require "NA_MDWeb_DistContactMapsProtsPDB.pl";

my $long=@ARGV;
if ($long < 4){
    print "Usage: perl $0 <pdbFile> <workdir> <output> <operation>\n";
    print "Operations:\n";
    print "\tALL\n\tCurves\n\tStiffness\n\tDistanceContactMaps\n";
    exit(0);
}
my ($mypdb, $workdir, $prefix, $operation) =@ARGV;

# 3-letter 20-Standard Residue Codes
my %residue_codes_std = (
ALA => 1,ARG => 1,ASN => 1,ASP => 1,
CYS => 1,CYX => 1,GLN => 1,GLU => 1,
GLY => 1,HIS => 1,HIE => 1,HID => 1,
HIP => 1,ILE => 1,LEU => 1,LYS => 1,
MET => 1,PHE => 1,PRO => 1,SER => 1,
THR => 1,TRP => 1,TYR => 1,VAL => 1,
WAT => 1,HOH => 1,'K+'=>1,'Cl-'=>1
);

my $dirScripts = "/app/Scripts/MCDNA/NAFlex";

# AMBER folder and Environmental Vars.
#my $achome = "/home/MCDNA/amber16";
#$ENV{'AMBERHOME'} = "/home/MCDNA/amber16";

# setenv LD_LIBRARY_PATH /opt/libgfortran1/lib/
#$ENV{'LD_LIBRARY_PATH'} = "$achome/lib:/opt/software/lib:/usr/lib/atlas:$ENV{'LD_LIBRARY_PATH'}";
#print "$ENV{'LD_LIBRARY_PATH'}\n";

#$ENV{'PATH'} = "$achome/bin:/opt/software/bin:$ENV{'PATH'}";
#print "$ENV{'PATH'}\n";

print "Workdir: $workdir\n";

$prefix = 'MCDNA_prova' if (!$prefix);

my $pdbFileIn = $mypdb;

mkdir("$workdir/$prefix") if(! -e "$workdir/$prefix");
mkdir("$workdir/$prefix/PDB") if(! -e "$workdir/$prefix/PDB");

`cp $pdbFileIn $workdir/$prefix/PDB`;

$pdbFileIn =~ s{.*/}{}; 

chdir("$workdir/$prefix/PDB");

if (-s "$pdbFileIn"){
	my $catpdb=''; 
	$catpdb=`zcat $pdbFileIn` if ($pdbFileIn =~ /gz$/);
	$catpdb=`cat $pdbFileIn` if ($pdbFileIn !~ /gz$/);
	my @models = &processPDB_MODELS($catpdb);
        for (my $model=1;$model<=$#models;$model++){

                my $newPdb = $models[$model];

		# Filtering and Adding chains:
		open TMP,">pdbTmp.addingChains.pdb";
		print TMP $newPdb;
		close TMP;

		print "perl $dirScripts/addChain.pl pdbTmp.addingChains.pdb pdbTmp.addedChains.pdb\n";
		`perl $dirScripts/addChain.pl pdbTmp.addingChains.pdb pdbTmp.addedChains.pdb`;

		$newPdb =  `cat pdbTmp.addedChains.pdb` if(-s "pdbTmp.addedChains.pdb"); 

                my @copies = &processPDB_CHAINS($newPdb);

                my $count = 1;
                print "Number of NA copies: $#copies\n";
                for(my $i=1;$i<=$#copies;$i+=2){
                        my $j = $i+1;
                        my $copy1 = $copies[$i];
                        my $copy2 = $copies[$j];

                        if(!$copy2){
                                print "FOUND just 1 copy in PDB\n";
                                next;
                        }

                        open NEW,">$pdbFileIn.$model.cp$count.pdb";
                        print NEW $copy1;
                        print NEW $copy2;
                        close NEW;

                        print "\t$pdbFileIn.$model.cp$count.pdb\n";
                        $count ++;
                }
        }
	if(-s "$pdbFileIn.1.cp1.pdb"){
		print "GOOD_ONE: $pdbFileIn.1.cp1.pdb\n";
	}
	else{
		print "Sorry, no proper double-helix DNA/RNA found in PDB $pdbFileIn...";
		exit;
	}
}
else{
	print "Sorry, we had problems getting the pdb file...\n";
	exit;
}
		
my $oldpdb = "$pdbFileIn.1.cp1.pdb";
my $pdb = "$pdbFileIn.renum.pdb";
&renumberPDB($oldpdb,$pdb);

my $pdb = `readlink -f $pdb`;
chomp($pdb);
print "PDB: $pdb\n";
my $top = $pdb;
$top =~s/.pdb/.top/g;
`$dirScripts/prepTOP.pl $pdb > $top`;

my $NAFleWorkDir = "$workdir/$prefix/PDB";

switch ($operation) {

        case 'Curves' {

		# case "Curves"
		print "Executing Curves Analysis (using Curves+ program).\n\n";
		print "run_curvesPDB($NAFleWorkDir,$prefix,$pdb)\n";
		run_curvesPDB($NAFleWorkDir,$prefix,$pdb);
	}

        case 'Stiffness' {

		# case "Stiffness"
		print "Executing Stiffness Analysis (using Curves+ program).\n\n";
		print "run_stiffnessPDB($NAFleWorkDir,$prefix,$pdb)\n";
		run_stiffnessPDB($NAFleWorkDir,$prefix,$pdb);
	}

	case 'DistanceContactMaps' {

		# case "DistanceContactMaps"
		print "Executing Distance Contact-Maps Analysis.\n\n";
		#print "run_distContactMapsPDB($NAFleWorkDir,$prefix,$pdb)\n";
		#run_distContactMapsPDB($NAFleWorkDir,$prefix,$pdb,$top);
		run_distContactMapsPDB($NAFleWorkDir,$prefix,$mypdb,$top);

	}

	case "ALL" {

		print "Executing Curves Analysis (using Curves+ program).\n\n";
		print "run_curves($NAFleWorkDir,$prefix,$pdb,$top)\n";
		if (! -e "CURVES"){
			run_curvesPDB($NAFleWorkDir,$prefix,$pdb,$top);
		}
		else{
			print "Curves already computed, skipping...\n";
		}

                print "Executing Distance Contact-Maps Analysis.\n\n";
                #print "run_distContactMapsPDB($NAFleWorkDir,$prefix,$pdb)\n";
                if (! -e "CONTACTS"){
                        run_distContactMaps($NAFleWorkDir,$prefix,$pdb,$top);
                }
                else{
                        print "Distance Contact Maps already computed, skipping...\n";
                }
	}

	else {

	    print "Usage: perl $0 <pdbFile> <topFile> <crdFile> <output> <operation>\n";
	    print "Operations:\n";
	    print "\tALL\n\tCurves\n\tDistanceContactMaps\n";

	}
}

########################################################################

sub processPDB_MODELS {
        my $pdb = shift;

        my @outpdb;

        my $model = 1;
        foreach (split "\n",$pdb){
                if (/^ENDMDL/){
                        $model++;
                }

                next if(!(/^ATOM/ or /^HETATM/ or /^TER/));

		if (/^TER/){
			$outpdb[$model].="$_\n";
		}

                my $at=substr($_,12,4);
                $at=~s/ //g;
                my $resn=substr($_,22,5);
                #$resn +=0;
                $resn=~s/ //g;
                my $res=substr($_,17,3);
                $res=~s/ //g;
                my $alt=substr($_,16,1);

                if($at!~/^H/ and $at!~/^\dH/ and $alt =~ /[ A]/){
                        #if($nucleic_codes_std{$res}){
                        if(!$residue_codes_std{$res} and $res ne 'HOH'){
                                $outpdb[$model].="$_\n";
                        }
                }
        }
        return @outpdb;
}

sub processPDB_CHAINS {
        my $pdb = shift;

        my @outpdb;

        my $chAnt = '';
        my $copy = 0;
        foreach (split "\n",$pdb){
                next if(!(/^ATOM/ or /^HETATM/ or /^TER/));

                my $at=substr($_,12,4);
                $at=~s/ //g;
                my $resn=substr($_,22,5);
                #$resn +=0;
                $resn=~s/ //g;
                my $res=substr($_,17,3);
                $res=~s/ //g;
                my $alt=substr($_,16,1);
                my $ch=substr($_,21,1);

                if ($ch ne $chAnt){
                        $copy++;
                }

                $outpdb[$copy].="$_\n";
                $chAnt = $ch;
        }
        return @outpdb;
}

# renumberPDB:
sub renumberPDB{
        my ($input,$output)=@_;

        my $cont=0;
        my $codeant='';
        my %output;
        open (OUT, ">$output");
        open FILE,"$input";
        while(<FILE>){
            next if (!(/^ATOM/ || /^TER/ || /^HETATM/));
            next if (/HOH/);

            my $mon=substr($_,17,3);
            $mon=~s/ //g;
            my $ch=substr($_,21,1);
            $ch=~s/ //g;
            my $resn=substr($_,22,5);
            $resn=~s/ //g;
            my $at=substr($_,12,4);
            $at=~s/ //g;

            my $code = "$mon$ch$resn";

                if($codeant ne $code){
                        $cont++;
                }

                my $newResn = sprintf("%4d",$cont);
                my $newLine = substr($_,0,22).$newResn." ".substr($_,27);

		print OUT $newLine;

                $codeant = $code;
        }
        close FILE;
	close OUT;
}

