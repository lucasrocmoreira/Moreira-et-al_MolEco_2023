#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=1000:00:00
#PBS -N fastsimcoal
#PBS -j oe
#PBS -m ae
#PBS -M lmoreira@amnh.org

# change to the working directory 
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

#####PARAMETERS#####
#par = prefix
#num - number of iterations

./run_fsc.sh -p $par -n $num -m > out.log
