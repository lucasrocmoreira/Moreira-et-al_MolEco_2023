#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N iqtree2
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

# Convert to iqtree POMO format

FastaToCounts.py SNP-only.fasta SNP-only.cf --iupac

iqtree -s SNP-only.cf -m HKY+P -b 100 -nt 16
