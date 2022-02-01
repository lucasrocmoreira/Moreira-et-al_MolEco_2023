#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N Chromosemble
#PBS -j oe
#PBS -m ae
#PBS -k oe

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

Chromosemble -t /home/lmoreira/reference/Other_references/Taeniopygia_guttata.whole.genome+chicken_W-2.fa \
-q /home/lmoreira/reference/Picoides_pubescens_ref-sorted.fa \
-o /home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2 \
-pseudochr 1
