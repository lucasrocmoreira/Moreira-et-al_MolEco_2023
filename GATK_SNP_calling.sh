#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N GATK_SNP_calling
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
ref=/home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2/pseudochromosomes-subset_sorted.fasta
## ref=reference sequence
## bam=bam file

name=`echo $bam | cut -d '/' -f2 | cut -d '.' -f1`
	
echo
echo "#######################"
echo $name
echo "#######################"

echo
echo
echo Calling SNPs with GATK Haplotype Caller
echo
echo

java -Djava.io.tmpdir=./tmp -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $ref \
-I $bam \
-o vcf_files/$name.raw.g.vcf \
--emitRefConfidence GVCF \
-minPruning 1 \
-minDanglingBranchLength 1 \
-hets 0.05 \
-nct 16
