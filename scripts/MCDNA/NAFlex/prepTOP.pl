#!/usr/bin/perl
use Data::Dumper;
use strict;
#
my $atomList;
my %RESIDUES;
my $residueList=[];
my %MASS = ('H' => 1.01, 'C' => 12.0, 'O'=>16.0, 'N'=>14.0, 'S' => 32.0, 'P' => 31.0);
while (<>) {
	next if (!/^ATOM/ && !/^HETATM/);
	my ($atom) = readPDBLine($_);
	push @$atomList,$atom;
	if (!$RESIDUES{$atom->{residueId}}) {
		my $newRes = {
			id => $atom->{residueId}, 
			firstAt => $#$atomList+1,
			name => substr($atom->{residueId},0,3),
			chain => substr($atom->{residueId},4,1),
			number => substr($atom->{residueId},5,4),
			icode => substr($atom->{residueId},9,1)			
			};
		$RESIDUES{$atom->{residueId}} = 1;
		push @$residueList, $newRes
	}
}
my ($listTypes, $types) = getTypes($atomList);
my $NATOM = scalar @$atomList;
my $NTYPES = scalar @$listTypes;
my $NRES = scalar @$residueList;
#
print "%VERSION  VERSION_STAMP = V0001.000  DATE = 09/20/06  19:23:43                  
%FLAG TITLE                                                                     
%FORMAT(20a4)
                                                                   
%FLAG POINTERS                                                                  
%FORMAT(10I8)                                                                   
";
my $format1="%8d" x 10 ."\n";
printf $format1, $NATOM, $NTYPES, 0, 0, 0, 0, 0, 0, 0, 0; 
printf $format1, 0, $NRES, 0, 0, 0, 0, 0, 0, 0, 0;
printf $format1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
printf "%8d\n", 0;
print "%FLAG ATOM_NAME                                                                 
%FORMAT(20a4)
";
my $txt='';
foreach my $at (@$atomList) {
	$txt.= substr($at->{name}."    ",0,4);
}                                           
printGroup($txt);                       
print "%FLAG CHARGE                                                                    
%FORMAT(5E16.8)
";
$txt='';
foreach my $at (@$atomList) {
	$txt .= sprintf "%16.8f",$at->{charge}*18.2223;
}
printGroup($txt);
print "%FLAG MASS                                                                      
%FORMAT(5E16.8)         
";                                                        
$txt='';
foreach my $at (@$atomList) {
	$txt .= sprintf "%16.8f",$MASS{substr($at->{name},0,1)};
}
printGroup($txt);
print"%FLAG ATOM_TYPE_INDEX                                                           
%FORMAT(10I8)
";
$txt='';
foreach my $at (@$atomList) {
	$txt .= sprintf "%8d", $types->{$at->{type}}
}                                                                  
printGroup($txt);
print "%FLAG NUMBER_EXCLUDED_ATOMS                                                     
%FORMAT(10I8)
";
printConstArray ("%8d",0, $NATOM);
print "%FLAG NONBONDED_PARM_INDEX                                                      
%FORMAT(10I8)
";
printConstArray ("%8d",0,$NTYPES * $NTYPES);
print "%FLAG RESIDUE_LABEL                                                             
%FORMAT(20a4)
";
$txt = "";
foreach my $res (@$residueList) {
	$txt .= substr ($res->{name}."    ",0,4);
}
printGroup($txt);
print "%FLAG RESIDUE_POINTER                                                           
%FORMAT(10I8)
";
$txt='';
foreach my $res (@$residueList) {
	$txt .= sprintf "%8d", $res->{firstAt};
}
printGroup($txt);
print "%FLAG BOND_FORCE_CONSTANT                                                       
%FORMAT(5E16.8)

%FLAG BOND_EQUIL_VALUE                                                          
%FORMAT(5E16.8)                                                                 

%FLAG ANGLE_FORCE_CONSTANT                                                      
%FORMAT(5E16.8)                                                                 

%FLAG ANGLE_EQUIL_VALUE                                                         
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_FORCE_CONSTANT                                                   
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_PERIODICITY                                                      
%FORMAT(5E16.8)                                                                 

%FLAG DIHEDRAL_PHASE                                                            
%FORMAT(5E16.8)                                                                 

%FLAG SOLTY
%FORMAT(5E16.8)                                                                 

";
printConstArray("%16.8f",0,$NTYPES);
print "%FLAG LENNARD_JONES_ACOEF                                                       
%FORMAT(5E16.8)

";
printConstArray("%16.8f",0,$NTYPES*($NTYPES+1)/2);
print "%FLAG LENNARD_JONES_BCOEF                                                       
%FORMAT(5E16.8)

";
printConstArray("%16.8f",0,$NTYPES*($NTYPES+1)/2);
print "%FLAG BONDS_INC_HYDROGEN                                                        
%FORMAT(10I8)                                                                   

%FLAG BONDS_WITHOUT_HYDROGEN                                                    
%FORMAT(10I8)                                                                   

%FLAG ANGLES_INC_HYDROGEN                                                       
%FORMAT(10I8)                                                                   

%FLAG ANGLES_WITHOUT_HYDROGEN                                                   
%FORMAT(10I8)                                                                   

%FLAG DIHEDRALS_INC_HYDROGEN                                                    
%FORMAT(10I8)                                                                   

%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                
%FORMAT(10I8)                                                                   

%FLAG EXCLUDED_ATOMS_LIST                                                       
%FORMAT(10I8)                                                                   

%FLAG HBOND_ACOEF                                                               
%FORMAT(5E16.8)                                                                 

%FLAG HBOND_BCOEF                                                               
%FORMAT(5E16.8)                                                                 

%FLAG HBCUT                                                                     
%FORMAT(5E16.8)                                                                 

";
print "%FLAG AMBER_ATOM_TYPE 
%FORMAT(20a4)

";
$txt='';
foreach my $at (@$atomList) {
	$txt .= substr($at->{type}."    ",0,4);
}
printGroup($txt);
print "%FLAG TREE_CHAIN_CLASSIFICATION                                                 
%FORMAT(20a4)

";
printConstArray("%s","E   ",$NATOM);
print "%FLAG JOIN_ARRAY                                                           
%FORMAT(10I8)

";
printConstArray("%8d",0,$NATOM);
print "%FLAG IROTAT                                                                    
%FORMAT(10I8)

";
printConstArray("%8d",0,$NATOM);
#print "%FLAG RADII                                                                     
#%FORMAT(5E16.8)                                                                 
#
#%FLAG SCREEN                                                                    
#%FORMAT(5E16.8)                                                                 
#
#";


sub readPDBLine {
	my ($line) = shift;
	my $newAt = {};
	my $newRes = {};
	$newAt->{name}=substr($line, 12,4);
	$newAt->{name}=~ s/ //g;
#Canvi d'ordre per atoms que comencin per numero
	$newAt->{name}=~ s/^([1-9])(.*)/$2$1/;
	$newAt->{altloc}=substr($line, 16,1);
	$newAt->{residueId}=substr($line,17,9);
	$newAt->{x}=substr($line,30,8);
	$newAt->{y}=substr($line,38,8);
	$newAt->{z}=substr($line,46,8);
	$newAt->{occ}=substr($line, 54,6);
	$newAt->{Bfact}=substr($line, 60,6);
	$newAt->{charge}=$newAt->{occ};
	$newAt->{type}=substr($newAt->{name},0,1);
	return $newAt;
}
	
	
sub getTypes {
	my $list = shift;
	my $typ = {};
	foreach my $at (@$list) {
		$typ->{$at->{type}}=1;
	}
	my $ityp = {};
	my $ltyp = [];
	foreach my $t (keys(%$typ)) {
		push @$ltyp, $t;
		$ityp->{$t} = $#$ltyp;
	}
	return ($ltyp,$ityp);
}
	
sub printGroup {
	my $txt = shift;
	$txt =~s/(.{80})/$1\n/g;
	print "$txt\n";                       
}

sub printConstArray {
	my ($format, $val, $num) = @_;
	my $txt='';
	foreach my $i (1 .. $num) {
		$txt .= sprintf $format, $val
	}
	printGroup ($txt);
}
