foreach f ( `ls --color=never *.helprms` )
echo $f
perl -ane {'print "$_" if ($F[5] > 40 or $F[5]<20)'} $f
end
