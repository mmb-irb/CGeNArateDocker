#!/bin/csh
cd /mnt/Web/CGN67d86f56ec511

echo "## CGeNArate ##"
perl /app/Scripts/MCDNA/runMCDNA_Prots_all.pl /mnt/Web/CGN67d86f56ec511/inputSequence.txt 50 3 proteins.mcdna.in 1 1

echo "## Analysis ##"
perl /app/Scripts/MCDNA/runMCDNA_Analysis.pl 3 1 1

echo "## Execute end-of-work routines ##"
wget -q -O- "https://mmb.irbbarcelona.org/CGNAW/end/jobs/67d86f56eca641.18630477"
