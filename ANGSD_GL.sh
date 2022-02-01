#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N ANGSD_GL
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

angsd -dovcf 1 -dopost 1 -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMaf 0.05 -SNP_pval 0.01 -minInd 50 -minMapQ 30 -minQ 20 -nThreads 16 -bam Picoides_pubescens.bamlist -out Picoides-pubescens2
