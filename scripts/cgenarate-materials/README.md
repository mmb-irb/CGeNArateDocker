# CGeNArate Materials

The programs provided are the following:

- CGeNArateCode/CGeNArate: To simulate linear DNA from a starting pdb structure, with settings file simulationssettings.txt, and run by executing run.sh, which reads the relevant name files from an appropriate namesfile****.dat. Pre-compiled static binaries are provided for unix machines (CGeNArate_st.exe) and MAC-M1 (CGeNArate_macM1.exe).

- CGeNArateCode/CGeNArateCircular: To simulate circular DNA from a starting pdb structure, with settings file simulationssettingsCircular.txt, and run by executing Circularrun.sh, which reads the relevant name files from an appropriate namesfile****.dat.Pre-compiled static binaries are provided for unix machines (CGeNArateCircular_st.exe) and MAC-M1 (CGeNArateCircular_macM1.exe).

- GLIMPS/Rebuild_nmer.py: To rebuild a simulated Coarse-Grained DNA. The input and output files have
to be edited in the first few lines of the python script.

- fdhelix/fdhelix.c: Contains a simplified version of the original code by David A. Case. fdhelix/MitochondriafromLine.py is a script to convert linear DNA into circular DNA with a specified supercoiling.
