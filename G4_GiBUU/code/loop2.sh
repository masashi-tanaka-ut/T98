#!/bin/bash   


for ((ifile=${1}; ifile<=${2}; ifile++))  
do
    echo "../output/B4_${ifile}.root"
    cp test.mac "../output/B4_${ifile}.mac" 
    echo "/random/setSeeds ${ifile} ${ifile}" >>  "../output/B4_${ifile}.mac" 
    echo "/analysis/setFileName ../output/B4_${ifile}.root" >> "../output/B4_${ifile}.mac" 
    echo "/run/beamOn $3" >> "../output/B4_${ifile}.mac" 
    ./exampleB4a -t 1 -m ../output/B4_${ifile}.mac > ../output/B4_${ifile}.log

done 
