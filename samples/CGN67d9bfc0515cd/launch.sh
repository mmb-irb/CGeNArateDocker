#!/bin/csh
cd /mnt/Web/CGN67d9bfc0515cd

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_Prots_all.pl /mnt/Web/CGN67d9bfc0515cd/inputSequence.txt 50 3 proteins.mcdna.in 0 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 3 0 1

echo "## Execute end-of-work routines ##"
wget -q -O- "http://localhost/public/end/jobs/67d9bfc0534745.72291737"
