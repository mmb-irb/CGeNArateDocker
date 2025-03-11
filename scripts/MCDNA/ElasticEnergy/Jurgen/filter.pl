
my $list = `ls --color=never *.helprms`;

my $bad = 0;
foreach my $file (split '\n',$list){
	my $bad = 0;
	my $count = 0;
	open FILE,"$file";
	while(<FILE>){
		chomp;
		my ($shift,$slide,$rise,$tilt,$roll,$twist) = split ' ';
		$bad = 1 if ($twist > 40 or $twist<20);
		$count++;
	}
	close FILE;
	print "$file $count\n" if (!$bad && $count);	
}

