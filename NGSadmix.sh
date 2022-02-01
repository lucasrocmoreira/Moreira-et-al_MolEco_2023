#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=5000:00:00
#PBS -N NGSadmix
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

echo
echo "#########################"
echo K = 1
echo "#########################"
echo

NGSadmix -likes /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/Picoides-pubescens_clean.beagle.gz -K 1 -P 32 -outfiles K1.Picoides_pubescens_NGSadmix -minMaf 0.05 -minInd 50

echo
echo "#########################"
echo K = 2
echo "#########################"
echo

NGSadmix -likes /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/Picoides-pubescens_clean.beagle.gz -K 2 -P 32 -outfiles K2.Picoides_pubescens_NGSadmix -minMaf 0.05 -minInd 50

echo
echo "#########################"
echo K = 3
echo "#########################"
echo

NGSadmix -likes /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/Picoides-pubescens_clean.beagle.gz -K 3 -P 32 -outfiles K3.Picoides_pubescens_NGSadmix.qopt -minMaf 0.05 -minInd 50

echo
echo "#########################"
echo K = 4
echo "#########################"
echo

NGSadmix -likes /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/ANGSD/Picoides-pubescens_clean.beagle.gz -K 4 -P 32 -outfiles K4.Picoides_pubescens_NGSadmix.qopt -minMaf 0.05 -minInd 50
