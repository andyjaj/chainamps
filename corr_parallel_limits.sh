#!/bin/bash

clear

#EDIT these if necessary
#parameter input filename
INPUT=input
#chain length
R=10
#number of measurements from 0 to R/2 (actually there is an implicit extra measurement for 0)
resolution=10

MINFILENUM=3
MAXFILENUM=1000
INCFILENUM=3
PRENAME=Ortho_

pi=3.14159265358979324

echo "Generating SPARSEMATRIX and BASIS files"
$PWD/STORE_MODEL.bin $INPUT

echo "Measuring intra vertex correlations from x=0 to x=R/2"
echo "R is" $R
echo "Resolution is" $resolution
echo "PI is" $pi

declare -a FILENAMES

COUNTER=0
for ((i=MINFILENUM; i<=MAXFILENUM;i+=INCFILENUM)); do
    FILENAMES+=$PRENAME$i.UNITCELL" "
    let COUNTER=COUNTER+1
done

echo $FILENAMES

for ((i=0; i<=resolution;i++)); do

    float=$(bc <<< "scale=4; $pi*$i/$resolution")

    $PWD/UNITCELL_MEASURE.bin -O $INPUT-Spin_Operator.SPARSEMATRIX -O $INPUT-Spin_Operator.SPARSEMATRIX@0@$float $FILENAMES

done
