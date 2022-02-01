#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N GATK_SNP_GVCF
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
## ref=reference sequence

ref=/home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2/pseudochromosomes-subset_sorted.fasta

echo
echo
echo Combine GVCFs
echo
echo

time java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $ref \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-10.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-1.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-2.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-3.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-4.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-5.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-6.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-7.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-8.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-AK-9.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-11.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-18.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-19.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-20.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-21.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-22.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-2.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-3.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-4.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-MW-7.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-26.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-37.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-38.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-39.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-40.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-42.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-43.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-47.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-48.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NE-49.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-01.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-02.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-03.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-04.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-05.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-06.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-07.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-08.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-09.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NR-10.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-10.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-11.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-12.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-13.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-15.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-16.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-18.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-5.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-8.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-NW-9.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-01.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-02.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-08.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-09.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-10.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-12.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-14.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-15.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-16.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SE-18.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-09.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-10.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-12.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-13.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-15.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-16.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-17.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-18.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-19.raw.g.vcf \
--variant /home/lmoreira/Picoides_pubescens/SNP_calling_from_pseudochrom/vcf_files/PP-SR-21.raw.g.vcf \
--heterozygosity 0.05 \
-o allsamples.raw.vcf
