#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N GATK_SNP_select_variants
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
## vcf
vcf=combined.db.vcf

echo
echo "#######################"
echo $vcf
echo "#######################"

echo
echo
echo Selecting Variants
echo
echo

vcftools \
--vcf $vcf \
--remove-indels \
--remove-filtered-all \
--min-meanDP 2 \
--max-meanDP 100 \
--maf 0.05 \
--max-missing 0.25 \
--hwe 0.01 \
--min-alleles 2 \
--max-alleles 2 \
--recode \
--recode-INFO-all \
--out SNP-only
