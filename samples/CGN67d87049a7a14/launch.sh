#!/bin/csh
cd /mnt/Web/CGN67d87049a7a14

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_Prots_all.pl /mnt/Web/CGN67d87049a7a14/inputSequence.txt 50 3 proteins.mcdna.in 0 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 3 0 1

echo "## Execute end-of-work routines ##"
wget -q -O- "https://mmb.irbbarcelona.org/CGNAW/end/jobs/67d87049a7f7a2.49210743"
