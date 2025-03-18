#!/bin/csh
cd /mnt/Web/CGN67d94200afbbb

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_all_new.pl /mnt/Web/CGN67d94200afbbb/inputSequence.txt 100 1 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 1 1 1

echo "## Execute end-of-work routines ##"
wget -q -O- "https://mmb.irbbarcelona.org/CGNAW/end/jobs/67d94200b03b99.28165192"
