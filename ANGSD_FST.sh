#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N ANGSD_FST
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

# This for loop calculates FST across all pairwise comparisons
for i in /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/SAF/*saf.idx;
do
	for j in /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/SAF/*saf.idx;
	do
	
	name1=`echo $i | cut -d '.' -f3`
	name2=`echo $j | cut -d '.' -f3`
	
	echo $name1 $name2

	echo "Calculate the 2dsfs prior"

	realSFS $i $j -P 16 > FST/$name1-$name2.sfs 

	echo "Prep for window analysis"

	realSFS fst index $i $j -sfs FST/$name1-$name2.sfs -fstout FST/$name1-$name2 -P 16

	echo "Get global estimate"

	realSFS fst stats FST/$name1-$name2.fst.idx -P 16

	echo "Sliding Windows"

	realSFS  fst stats2 FST/$name1-$name2.fst.idx -win 50000 -step 10000 -P 16 > FST/$name1-$name2.slidingwindow.fst
	
	done
done
