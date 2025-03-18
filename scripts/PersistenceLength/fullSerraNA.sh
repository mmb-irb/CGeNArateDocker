#!/bin/bash
##GPUCLUSTER
## "56merLpredictednoh"
name=$1

#2. If the structure is double-stranded ("2") or if it is single-stranded ("1") 
#3. If the structure is linear ("1") or closed ("2")

echo "c Input file for SerraNA                              
c Lines beginning with c are comment lines.
c Always leave one empty line before entering a parameter.
c Do not change the order of inputs.
c------------------------------------------------------------------------------
c Indicate if the structure to be analysed is double or single stranded
c Write 1 for single and 2 for double.

   2
c------------------------------------------------------------------------------
c Indicate if the structure is linear or cicular.
c Write 1 for linear and 2 for circular

   1
c -------------------------------------------
c Topology file. The file must be an amber topology file.                                       
c                                          

   input/${name}.prmtop
c -------------------------------------------
c Trajectory file. The file must be an amber trajectory file.                                         
c                                          

   input/${name}.mdcrd
c
c No more inputs to read
c end" > s_NA.in

#echo "export LD_LIBRARY_PATH=/orozco/projects/BioExcel/Test_Kim_Fede/PL_Agnes
#/orozco/projects/BioExcel/Test_Kim_Fede/PL_Agnes/SerraNA < s_NA.in"

#export LD_LIBRARY_PATH=/orozco/projects/BioExcel/Test_Kim_Fede/PL_Agnes
#/orozco/projects/BioExcel/Test_Kim_Fede/PL_Agnes/SerraNA < s_NA.in
/app/Scripts/PersistenceLength/PL_Agnes/SerraNA < s_NA.in

#echo 23
#exit 0

mv BPP.out output/BPP${name}.out
mv BSP.out output/BSP${name}.out
mv structural_parameters.out output/structural_parameters${name}.out
mv elastic_parameters.out output/elastic_parameters${name}.out

echo "Finished SerraNA"


echo "c Input file for Analysis programs                        
c Lines beginning with c are comment lines.
c Always leave one empty line before entering a parameter.
c Do not change the order of inputs.
c------------------------------------------------------------------------------
c Indicate ranges of Persistence lengths, Tilt and Roll.
c You should write 4 integers: a, b, c, d
c Default a=b=0 corresponds to considering the whole fragment.
c Default c=d=0 corresponds to default methodology (see manual).

   0, 0, 0, 0
c------------------------------------------------------------------------------
c Indicate ranges of Twist
c

   0, 0, 0, 0
c------------------------------------------------------------------------------
c Indicate ranges of Stretch
c

   0, 0, 0, 0
c -------------------------------------------
c Elastic parameters file (elastic_parameters.out). Note that if you analysed
c one single snapshot with SerraNA, you won't have this file...                             
c                                          

   output/elastic_parameters${name}.out
c -------------------------------------------
c Structural parameters file (structural_parameters.out)                                   
c                                          

   output/structural_parameters${name}.out
c
c No more inputs to read
c end" > ov_NA.in

#/orozco/projects/BioExcel/Test_Kim_Fede/PL_Agnes/Analysis < ov_NA.in > Analysis/Analysis_${name}.out
/app/Scripts/PersistenceLength/PL_Agnes/Analysis < ov_NA.in > Analysis/Analysis_${name}.out

echo "Finished Analysis"

exit 0

###############################################################


echo "c Input file for Extract program                        
c Lines beginning with c are comment lines.
c Always leave one empty line before entering a parameter.
c------------------------------------------------------------------------------
c Path to data file (BPP.out BSP.out, structural_parameters.out, elastic_parameters.out)
c

   output/structural_parameters${name}.out
c------------------------------------------------------------------------------
c Indicate if you want to extract a sublength [l] by typing     :    0
c Or if you want to get overalls from a particular region [a,b] :    1
c Any other values different to 0 or 1 will caused an error
c

  1
c------------------------------------------------------------------------------
c If you typed 0, then indicate the length [l] you want to process.
c This length must be:  0 < l < N, where N is the number of bp-steps.
c
c If you typed 1, then indicate the region (a,b) from which you want 
c to extract avg+-sd as a function of length. 
c If the structure is opened (linear DNA), then 0 < a < b < N.
c If the structure is closed (circular DNA), then both a < b or b < a, are valid
c since those ranges correspond to different regions along a closed structure.
c Default a=b=0 consider all possible lengths in the whole fragment,
c The only invalid input is a = b, if both a and b are different than 0.
c

   0,0
c -------------------------------------------
c No more inputs to read
c end" > ex_NA.in

exit 0

SerraNA/Extract < ex_NA.in 

mv structural_plot.out output/structural_plot${name}.out


echo "Finished Extraction"

: '
1.- Path to either BPP, BSP, structural or elastic parameters file. If you selected to extract BPP or BSP, then all other inputs will be ignored. 2.- Type "0" for extracting a sub-length or "1" for getting avg+-sd as a function of length. 3.- (i) If you type 0, then indicate the length (l) you want to process, which should be 0 < l < N, where N is the number of bp-steps. (ii) If you typed 1, then indicate the region (a,b) from which you want to extract avg+-sd as a function of length. If it is linear DNA, then 0 < a < b < N If it is circular DNA, then both a < b or b < a, are valid DEFAULT OPTION, a=b=0, consider all possible lengths in the whole fragment.

The program can create different types of outputs:

    BPP_plot.out, if BPP.out is processed
    BSP_plot.out, if BSP.out is processed
    structural_lmer.out, if structural_parameters.out is processed to extract a particular sub-length l
    structural_[a:b].out, if structural_parameters.out is processed to extract a length-dependence for a particular sub-fragment
    structural_plot.out, if a=b=0
    elastic_lmer.out, if elastic_parameters.out is processed to extract a particular sub-length l
    elastic_[a:b].out, if elastic_parameters.out is processed to extract a length-dependence for a particular sub-fragment
    elastic_plot.out, if a=b=0
'

