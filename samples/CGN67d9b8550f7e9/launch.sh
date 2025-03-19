#!/bin/csh
cd /mnt/Web/CGN67d9b8550f7e9

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_all_new.pl /mnt/Web/CGN67d9b8550f7e9/inputSequence.txt 100 0 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 1 0 1

echo "## Execute end-of-work routines ##"
wget -q -O- "http://localhost/public/end/jobs/67d9b855102df1.14756966"
