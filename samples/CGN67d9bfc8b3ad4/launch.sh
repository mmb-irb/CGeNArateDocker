#!/bin/csh
cd /mnt/Web/CGN67d9bfc8b3ad4

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_Prots_all.pl /mnt/Web/CGN67d9bfc8b3ad4/inputSequence.txt 50 3 proteins.mcdna.in 1 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 3 1 1

echo "## Execute end-of-work routines ##"
wget -q -O- "http://localhost/public/end/jobs/67d9bfc8b5dd99.99963149"
