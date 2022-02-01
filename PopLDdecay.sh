#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N PopLDdecay
#PBS -j oe
#PBS -m ae
#PBS -M lmoreira@amnh.org
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

vcf=/home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/SNP-only.recode.vcf

echo
echo "#######################"
echo $vcf
echo "#######################"

PopLDdecay -InVCF  $vcf  -OutStat AK -SubPop AK 
PopLDdecay -InVCF  $vcf  -OutStat NE -SubPop NE
PopLDdecay -InVCF  $vcf  -OutStat MW -SubPop MW
PopLDdecay -InVCF  $vcf  -OutStat SE -SubPop SE
PopLDdecay -InVCF  $vcf  -OutStat NR -SubPop NR
PopLDdecay -InVCF  $vcf  -OutStat SR -SubPop SR
PopLDdecay -InVCF  $vcf  -OutStat NW -SubPop NW
