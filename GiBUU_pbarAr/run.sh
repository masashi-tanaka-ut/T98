#/bin/bash
# sh run.sh "directory" "Energy" "Seed"
mkdir $1
cd $1
sed "s/0.03/$2/g" ../001_antiP_Ar_seed1.job | sed "s/11111/$3/g" | sed "s/-999.9/0.0/g" > job.job
../../GiBUU.x < job.job > log.log
