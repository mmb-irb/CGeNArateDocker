#!/bin/csh
cd /mnt/Web/CGN67d9b89285e4a

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_all_new.pl /mnt/Web/CGN67d9b89285e4a/inputSequence.txt 100 1 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 1 1 1

echo "## Execute end-of-work routines ##"
wget -q -O- "http://localhost/public/end/jobs/67d9b89286d769.61804502"
