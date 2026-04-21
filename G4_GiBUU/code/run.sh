#!/bin/bash   


ifile=$1
ievent=$2
iseed=$3

echo "../output/B4_${ifile}.root"
cp test.mac "../output/B4_${ifile}.mac" 
echo "/random/setSeeds ${iseed} ${iseed}" >>  "../output/B4_${ifile}.mac" 
echo "/analysis/setFileName ../output/B4_${ifile}.root" >> "../output/B4_${ifile}.mac" 
echo "/run/beamOn ${ievent}" >> "../output/B4_${ifile}.mac" 
./exampleB4a -t 1 -m ../output/B4_${ifile}.mac > ../output/B4_${ifile}.log

