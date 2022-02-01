#!/bin/bash
#PBS -l select=1:ncpus=16:mem=200Gb
#PBS -l walltime=5000:00:00
#PBS -N ReLERNN
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

SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
CPU="16"
MU="4.007e-9"
RTR="1"
DIR="./E.highMU"
VCF="/home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/SNP-only.E.maf002.recode.vcf"
GENOME="/home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ReLERNN/pseudochromosomes-subset_sorted.chrom.bed"
MASK="/home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ReLERNNy/pseudochromosomes-subset_sorted.chrom.Nmask.bed"

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --mask ${MASK} \
    --unphased \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --nCPU ${CPU}

# Train network
${TRAIN} \
    --projectDir ${DIR} \

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR}

# Parametric Bootstrapping
${BSCORRECT} \
    --projectDir ${DIR} \
    --nCPU ${CPU} \
    --nSlice 2 \
    --nReps 2
