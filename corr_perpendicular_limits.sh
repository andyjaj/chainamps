#!/bin/bash

clear

#EDIT these if necessary
#parameter input filename
INPUT=input
#number of measurements from 0 to R/2 (actually there is an implicit extra measurement for 0)
MAXSEP=20

MINFILENUM=3
MAXFILENUM=1000
INCFILENUM=3
PRENAME=Ortho_

#uncomment if you want the script to also generate the SPARSEMATRIX files
#$PWD/STORE_MODEL.bin $INPUT

echo "Measuring inter vertex correlations" 
echo "Max Separation is" $MAXSEP

declare -a FILENAMES

COUNTER=0
for ((i=MINFILENUM; i<=MAXFILENUM;i+=INCFILENUM)); do
    FILENAMES+=$PRENAME$i.UNITCELL" "
    let COUNTER=COUNTER+1
done

echo $FILENAMES
    
for ((i=7; i<=MAXSEP;i++)); do

    SEP=$i

    $PWD/UNITCELL_MEASURE.bin -O $INPUT-Spin_Operator.SPARSEMATRIX -O $INPUT-Spin_Operator.SPARSEMATRIX -S $SEP $FILENAMES

done
