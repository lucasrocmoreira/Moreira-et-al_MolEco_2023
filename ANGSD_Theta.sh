#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=5000:00:00
#PBS -N ANGSD_Theta
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

date
time

ref=/home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2/pseudochromosomes-subset_sorted.fasta

for pop in AK NE MW SE SR NR NW E R;
do

	# Creating saf files
	angsd -GL 2 -dosaf 1 -anc $ref -minMapQ 30 -minQ 20 -nThreads 32 -rf regions.txt -bam Picoides_pubescens.$pop.bamlist -out theta_every_pop/Picoides-pubescens.$pop

	# Estimating sfs for each 10 Mb 
	realSFS theta_every_pop/Picoides-pubescens.$pop.saf.idx -nSites 10000000 -P 32 -fold 1 > theta_every_pop/Picoides_pubescens.$pop.sfs

	# Summing sfs across all windows
	perl -lane '$sum[$_] += $F[$_] for 0..$#F; END {print join $", @sum}' theta_every_pop/Picoides_pubescens.$pop.sfs > theta_every_pop/Picoides_pubescens.$pop.total.sfs

	# Estimating theta
	realSFS saf2theta theta_every_pop/Picoides-pubescens.$pop.saf.idx -sfs theta_every_pop/Picoides_pubescens.$pop.total.sfs -fold 1 -outname theta_every_pop/Picoides_pubescens.$pop
	thetaStat do_stat theta_every_pop/Picoides_pubescens.$pop.thetas.idx

	rm theta_every_pop/Picoides_pubescens.$pop.thetas.gz theta_every_pop/Picoides_pubescens.$pop.thetas.idx theta_every_pop/Picoides-pubescens.$pop.saf.*

done
