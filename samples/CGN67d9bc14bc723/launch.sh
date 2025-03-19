#!/bin/csh
cd /mnt/Web/CGN67d9bc14bc723

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_circular_all.pl /mnt/Web/CGN67d9bc14bc723/inputSequence.txt 10 5 -1 25000000 1 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 2 1 1

echo "## Execute end-of-work routines ##"
wget -q -O- "http://localhost/public/end/jobs/67d9bc14bdaf24.23726253"
