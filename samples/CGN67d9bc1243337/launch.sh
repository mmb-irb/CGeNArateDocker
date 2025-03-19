#!/bin/csh
cd /mnt/Web/CGN67d9bc1243337

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_circular_all.pl /mnt/Web/CGN67d9bc1243337/inputSequence.txt 10 5 -1 25000000 0 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 2 0 1

echo "## Execute end-of-work routines ##"
wget -q -O- "http://localhost/public/end/jobs/67d9bc12443dc0.80753534"
