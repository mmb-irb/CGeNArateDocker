#!/bin/csh
cd /mnt/Web/CGN67d87509b8772

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_circular_all.pl /mnt/Web/CGN67d87509b8772/inputSequence.txt 10 5 -1 25000000 1 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 2 1 1

echo "## Execute end-of-work routines ##"
wget -q -O- "https://mmb.irbbarcelona.org/CGNAW/end/jobs/67d87509b8d2b2.67967905"
