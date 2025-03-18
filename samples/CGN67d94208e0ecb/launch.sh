#!/bin/csh
cd /mnt/Web/CGN67d94208e0ecb

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_circular_all.pl /mnt/Web/CGN67d94208e0ecb/inputSequence.txt 10 5 -1 25000000 0 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 2 0 1

echo "## Execute end-of-work routines ##"
wget -q -O- "https://mmb.irbbarcelona.org/CGNAW/end/jobs/67d94208e19073.61173417"
