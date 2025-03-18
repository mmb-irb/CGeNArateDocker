#!/bin/csh
cd /mnt/Web/CGN67d86f36f0e92

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_all_new.pl /mnt/Web/CGN67d86f36f0e92/inputSequence.txt 100 0 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 1 0 1

echo "## Execute end-of-work routines ##"
wget -q -O- "https://mmb.irbbarcelona.org/CGNAW/end/jobs/67d86f37021ff2.50529559"
