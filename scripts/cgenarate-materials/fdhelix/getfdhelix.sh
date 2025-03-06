#!/bin/bash

#first, compile: "bash compilefdhelix.sh"

#if sequence is short: "./fdhelix.exe abdna AAAAAAAAAATTTTTTTTTTATATATATATAAAAAAAAAATTTTTTTTTT | grep C1 > ATAT.pdb"
#if sequence is in file: "bash getfdhelix.sh seqfile.txt ATAT.pdb"

var=$(cat $1)

#./fdhelix.exe abdna $var | grep C1 > $2
#fdhelix.exe abdna $var | grep C1 > $2
fdhelix.exe abdna $var | egrep 'C1|TER' > $2
#fdhelix.exe abdna AAAA | grep C1 > $2
