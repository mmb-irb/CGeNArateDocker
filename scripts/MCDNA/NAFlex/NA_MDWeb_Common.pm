#!/usr/bin/perl -w

package NA_MDWeb_Common;
use base 'Exporter';
use strict;

# PCAsuite Environment
$ENV{'LD_LIBRARY_PATH'} = "/home/NAFlex/lib:$ENV{'LD_LIBRARY_PATH'}";

our @EXPORT    = qw(
	%nucleic_codes %nucleic_codes_charmm %bases
	nucleicType checkPDB exeCanal exeCurves exeCurvesPDB getAvg
	$softdir $scriptsDir 
	$curvesDir $curves $canal
	$HelicalParamsTable $HelicalBaseParamsTable
	%translation
	$stiffnessTable $amberexe $noes $jcoupling
	$pcaunzip $pcazip $pczdump $preptop 
);

# Auxiliary Folders
our $scriptsDir = "/app/Scripts/MCDNA/NAFlex";
our $softdir = "/app/Scripts/MCDNA/NAFlex/soft";
our $amberexe = "/home/MCDNA/amber16/bin";
#our $curvesDir = "$softdir/Curves3.0";
our $curvesDir = "/opt/conda/envs/glimps_env/.curvesplus";
#our $curves = "$curvesDir/Cur+";
our $curves = "Cur+";
#our $canal = "$curvesDir/canal";
our $canal = "Canal";
our $HelicalBaseParamsTable = "$scriptsDir/BasePairParamsAvg.table";
our $HelicalParamsTable = "$scriptsDir/BasePairStepParamsAvg.table";
our $stiffnessTable = "$scriptsDir/BasePairParamsStiffness.table";
our $noes = "$softdir/nmr/NOEsClassification";
our $jcoupling = "$softdir/nmr/Jcouplings";
#our $pcaunzip = "$softdir/pcazip/pcaunzip64"; # 64 bits version with netcdf libraries
our $pcaunzip = "pcaunzip"; 
#our $pcazip = "$softdir/pcazip/pcazip64";  # 64 bits version with netcdf libraries
our $pcazip = "pcazip"; 
#our $pczdump = "$softdir/pcazip/pczdump";
our $pczdump = "pczdump";
our $preptop = "$scriptsDir/prepTOP.pl";

our %translation;
$translation{'Shear'} = 1;
$translation{'Stretch'} = 1;
$translation{'Stagger'} = 1;
$translation{'Shift'} = 1;
$translation{'Slide'} = 1;
$translation{'Rise'} = 1;
$translation{'shear'} = 1;
$translation{'stretch'} = 1;
$translation{'stagger'} = 1;
$translation{'shift'} = 1;
$translation{'slide'} = 1;
$translation{'rise'} = 1;

# 3-letter 20-Standard Residue Codes
our %residue_codes_std = (
        ALA => 1,
        ARG => 1,
        ASN => 1,
        ASP => 1,
        CYS => 1,
        GLN => 1,
        GLU => 1,
        GLY => 1,
        HIS => 1,
        ILE => 1,
        LEU => 1,
        LYS => 1,
        MET => 1,
        PHE => 1,
        PRO => 1,
        SER => 1,
        THR => 1,
        TRP => 1,
        TYR => 1,
        VAL => 1,
	NN => 1
);

our %residue_codes_charmm = (
                %residue_codes_std,
        HSE => 1,       # CHARMM
        HSD => 1,       # CHARMM
        HSP => 1,       # CHARMM
        LSN => 1,       # CHARMM
        ASPP => 1,      # CHARMM
        GLUP => 1       # CHARMM
);

# 3-letter AMBER Residue Codes
our %residue_codes_amber = (
                %residue_codes_std,
        CYX => 1,       # AMBER
        CYM => 1,       # AMBER
        CYN => 1,       # AMBER
        LYN => 1,       # AMBER
        LYP => 1,       # AMBER
        ASH => 1,       # AMBER
        GLH => 1,       # AMBER
        HIE => 1,       # AMBER
        HID => 1,       # AMBER
        HIP => 1        # AMBER
);

# 3/4-letter GMX Residue Codes
our %residue_codes_gmx = (
                %residue_codes_std,
        HISA => 1,      # GMX (Gromos, opls)
        HISB => 1,      # GMX (Gromos, opls)
        HISH => 1,      # GMX (Gromos, opls)
        LYSH => 1,      # GMX (Gromos, opls)
	CYS1 => 1,      # GMX (Gromos, opls)
        CYS2 => 1,      # GMX (Gromos, opls)
        CYSH => 1,      # GMX (Gromos, opls)
        ASPH => 1,      # GMX (Gromos, opls)
        GLUH => 1,      # GMX (Gromos, opls)
        ARGN => 1       # GMX (Gromos, opls)
);

# 3/4-letter Residue Codes (ALL)
our %residue_codes = (%residue_codes_std, %residue_codes_charmm, %residue_codes_amber, %residue_codes_gmx);

# Nucleic Bases
our %nucleic_codes_std = (
        "A" => 1,
        "T" => 1,
        "C" => 1,
        "G" => 1,
        "U" => 1,
        "DA" => 1,
        "DT" => 1,
        "DC" => 1,
        "DG" => 1,
        "DH" => 1,
        "DF" => 1,
        "RU" => 1,
        "RC" => 1,
        "RA" => 1,
        "RG" => 1,
        "A5" => 1,
        "A3" => 1,
        "C5" => 1,
        "C3" => 1,
        "T5" => 1,
        "T3" => 1,
        "G5" => 1,
        "G3" => 1,
        "U5" => 1,
        "U3" => 1,
        "DA5" => 1,
        "DA3" => 1,
        "DC5" => 1,
        "DC3" => 1,
        "DT5" => 1,
        "DT3" => 1,
        "DG5" => 1,
        "DG3" => 1,
        "RU5" => 1,
        "RU3" => 1,
        "RC5" => 1,
        "RC3" => 1,
        "RA5" => 1,
        "RA3" => 1,
        "RG5" => 1,
        "RG3" => 1
);

# Nucleic Bases Charmm
our %nucleic_codes_charmm = (
        "ADE" => 1,
        "THY" => 1,
        "CYT" => 1,
        "GUA" => 1,
        "URA" => 1
);

# Nucleic Bases Codes (ALL)
our %nucleic_codes = (%nucleic_codes_std, %nucleic_codes_charmm);

# Nucleic Bases Atoms
our %bases = (
        "C3*" => 1,
        "O3*" => 1,
        "C4*" => 1,
        "C5*" => 1,
        "05*" => 1,
        "C3'" => 1,
        "O3'" => 1,
        "C4'" => 1,
        "C5'" => 1,
        "05'" => 1,
        "P" => 1
);

# 3-letter to 1-letter Nucleic Codes
my %oneLetter = (
        "ADE" => 'A',
        "THY" => 'T',
        "CYT" => 'C',
        "GUA" => 'G',
        "URA" => 'U'
);

sub nucleicType{

        my ($pdb) = @_;

        my $type = 'DNA';
        open IN,"$pdb";
        while(<IN>){
                chomp;

                if(/^ATOM/){

                        my $mon = substr($_,17,3);
                        $mon =~s/ //g;
                        my $at = substr($_,12,4);
                        $at =~s/ //g;

                        if(!$nucleic_codes{$mon}){
                                #$type = 'Prot';
                        }

                        if ($at eq "O2'" or $at eq "O2*"){
                                $type = 'RNA';
                        }
                }
        }
        close IN;

        return $type;
}

# checkPDB
# Getting info about chains, segments and residues.
sub checkPDB{
        my ($input) = @_;

        my $model = 0;
        my $chAnt = '';
        my $resAnt = '';
        my $resnAnt = '';
        my %residues;
        my @firsts;
        my @lasts;
        my $seq1 = '';
        my $seq2 = '';
        my %done;
        my $first = 1;
        my $strand = 0;

 
	# Adding Chains if needed
	`cp $input ./input.pdb`;
	my $isThereAchain = `grep ATOM input.pdb | cut -c22-22 | sort -u`;
	
	if ($isThereAchain == ''){
		`cp input.pdb input.pdb.tmp`;
		`perl $scriptsDir/addChain2.pl input.pdb.tmp input.pdb`;
	}

        open PDB,"input.pdb";
        while (<PDB>){
            chomp;

            if (/^MODEL/){
                last if($model);
                $model = substr($_,6,10)+0;
            }
            #next if (!(/^ATOM/ || /^TER/ || /^HETATM/));
            next if (!(/^ATOM/ || /^HETATM/));

	    next if (/WAT/ || /HOH/ || /Na+/ || /Cl-/ || /K+/);

            my $mon=substr($_,17,3);
            $mon =~ /([ACTGUHF])/;
            my $base = $1;
            #if(!$base){
            #    $mon =~ /([ACTGU])/;
            #    $base = $1;
            #}

            $mon=~s/ //g;

	    if(!$nucleic_codes_std{$mon}){
		$base = "X";
	    }
  	    if($nucleic_codes_charmm{$mon}){
		$base = $oneLetter{$mon};
	    }
	    if ($residue_codes{$mon}){
		$base = "P";	
	    }

	    # WLC Coarse-Grained Bead Residue Name
	    #if($mod eq 'Mid'){
		#$base = "B";	
	    #}

            my $ch=substr($_,21,1);
            $ch=~s/ //g;
            my $resn=substr($_,22,5);
            $resn=~s/ //g;
            my $at=substr($_,12,4);
            $at=~s/ //g;

            my $code = "$mon$ch$resn";

                # Saving Sequence (only first strand).
                #$seq1.="$base" if(!$done{$code} and !$strand and !/^TER/ and $base ne "P");
                #$seq2.="$base" if(!$done{$code} and $strand and !/^TER/ and $base ne "P");
                $seq1.="$base" if(!$done{$code} and !$strand and $base ne "P");
                $seq2.="$base" if(!$done{$code} and $strand and $base ne "P");
                $done{$code} = 1;

                if(!$ch){
                        $ch = "@";
                }

                # Terminal Residues/Bases
                #if(/^TER/ or ( ($ch ne $chAnt) and ($chAnt ne '') ) ){
                if( ($ch ne $chAnt) and ($chAnt ne '') ){
                        $residues{$chAnt}->{"$resAnt:$resnAnt"} = 2 if ($resAnt ne '');
                        $residues{$ch}->{"$mon:$resn"} = 1 if($mon);
                        if(!$strand){
                                $strand = 1; # if (!$first);
                                $seq2 = $seq2.substr($seq1,-1);
                                $seq1 = substr($seq1,0,-1);
                        }
                }
                if($first or $mon =~ /5$/ or $at eq 'H5T'){
                        $residues{$ch}->{"$mon:$resn"} = 1;
                        $first = 0;
                }
                if($mon =~ /3$/ or $at eq 'H3T'){
                        $residues{$ch}->{"$mon:$resn"} = 2;
                        $strand = 1;
                }

                $resAnt = $mon;
                $resnAnt = $resn;
                $chAnt = $ch;
        }
        close PDB;

        # Last Terminal Residue/Base
        $residues{$chAnt}->{"$resAnt:$resnAnt"} = 2;

        foreach my $ch (sort keys %residues){
                foreach my $code (sort keys %{$residues{$ch}}){
                        my $v = $residues{$ch}->{$code};
			print "First: $ch-$code\n" if($v == 1);
			print "Lasts: $ch-$code\n" if($v == 2);
                        push @firsts,"$ch-$code" if($v == 1);
                        push @lasts,"$ch-$code" if($v == 2);
                }
        }

	@firsts = sort sortCh @firsts;
	@lasts = sort sortCh @lasts;

        $seq2 = reverse $seq2;

        return ($seq1,$seq2,\@firsts,\@lasts);
}

# sortCh
# Sorting chain
sub sortCh {
	# @-A5:19
	my ($ch1,$res1) = split ':',$a;
	my ($ch2,$res2) = split ':',$b;
	return $res1 <=> $res2;
}

# exeCurves
# Executing Curves program from Curves+ package
sub exeCurves {
        my ($inputCrd,$inputTop,$output,$iniSnap,$endSnap,$jumpSnap,$strand1_init,$strand1_end,$strand2_init,$strand2_end) = @_;

        # Building Curves.in configuration file.
        open IN,">curves.in";
        print IN "\&inp\n";
        print IN "  file=$inputCrd,\n";
        print IN "  lis=$output,\n";
        print IN "  lib=$curvesDir/standard,\n";
        print IN "  ftop=$inputTop,\n";
        print IN "  itst=$iniSnap,itnd=$endSnap,itdel=$jumpSnap,\n";
        #print IN "  circ=T,\n";
        print IN "\&end\n";
        print IN "2 1 -1 0 0\n";
        print IN "$strand1_init:$strand1_end\n";
        print IN "$strand2_end:$strand2_init\n";
        close IN;

        # Executing Curves+
        `$curves < curves.in &> curves.log`;
}

# exeCurvesPDB
# Executing Curves program from Curves+ package from just a PDB
sub exeCurvesPDB {
        my ($inputPDB,$output,$strand1_init,$strand1_end,$strand2_init,$strand2_end) = @_;

        # Building Curves.in configuration file.
        open IN,">curves.in";
        print IN "\&inp\n";
        print IN "  file=$inputPDB,\n";
        print IN "  lis=$output,\n";
        print IN "  lib=$curvesDir/standard,\n";
        print IN "\&end\n";
        print IN "2 1 -1 0 0\n";
        print IN "$strand1_init:$strand1_end\n";
        print IN "$strand2_end:$strand2_init\n";
        close IN;

        # Executing Curves+
        `$curves < curves.in &> curves.log`;
}

# exeCanal
# Executing Canal program from Curves+ Package
sub exeCanal {
        my ($outputCanal,$outputCurves,$seq,$histo,$series,$corr) = @_;

        # Building Curves.in configuration file.
        open IN,">canal.in";
        print IN "\&inp\n";
        print IN "  lis=$outputCanal,\n";
        print IN "  lev1=0,lev2=0,\n";
        print IN "  histo=.t.,\n" if ($histo);
        print IN "  series=.t.,\n" if ($series);
        print IN "  corr=.t.,\n" if ($corr);
        print IN "\&end\n";
        print IN "$outputCurves $seq\n";
        close IN;

        # Executing Curves+
        `$canal < canal.in &> canal.log`;
}

# getAvg
# Computing Mean and Stdev from a set of values in a 2-rows file.
sub getAvg{

        my ($file) = @_;

        my $mean = 0;
        my $stdev = 0;
        my %info;

        open FILE,"$file";
        while(<FILE>){
                # 1     0.42
                # 2     1.25
                chomp;
                my ($snapshot,$value,$rest) = split ' ';
                $info{$snapshot} = $value;
                $mean += $value;
        }
        close FILE;

        my $cont = keys %info;
        $mean /= $cont if($cont);

        foreach my $snp (sort keys %info){
                my $v = $info{$snp};
                my $add = $v - $mean;
                my $add2 = $add * $add;
                $stdev += $add2;
        }
        $stdev /= $cont if($cont);
        $stdev = sqrt($stdev);

        return ($mean,$stdev);
}

