#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N Indel_Realigner
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

#Arguments:
## bam = clipped BAM file

ref=/home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2/pseudochromosomes-subset_sorted.fasta
name=`echo $bam | cut -d '.' -f1`

echo
echo "#######################"
echo $name
echo "#######################"

echo
echo
echo Indel Realigner
echo
echo

java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $ref \
-I $bam \
-targetIntervals /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/all_samples.intervals \
-o /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/BAM_files/$name.chrom.realigned.bam
