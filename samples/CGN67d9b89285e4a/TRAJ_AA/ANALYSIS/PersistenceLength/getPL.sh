#!/bin/bash
##GPUCLUSTER

#names="56merLpredictednoh"

names=$1

bash fullSerraNA.sh ${names}

file=Analysis/Analysis_${names}.out
var=$(grep "A \[d\]" ${file})


echo $file
echo "Persistence length", $names, ${var:14:16}

