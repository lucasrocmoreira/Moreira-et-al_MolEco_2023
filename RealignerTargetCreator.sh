#!/bin/bash
#PBS -l select=1:ncpus=4
#PBS -l walltime=5000:00:00
#PBS -N RealignerTargetCreator
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
echo RealignerTargetCreator
echo
echo

time java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $ref \
-I ../BAM_files/PP-AK-10.chrom.dedup.bam \
-I ../BAM_files/PP-AK-1.chrom.dedup.bam \
-I ../BAM_files/PP-AK-2.chrom.dedup.bam \
-I ../BAM_files/PP-AK-3.chrom.dedup.bam \
-I ../BAM_files/PP-AK-4.chrom.dedup.bam \
-I ../BAM_files/PP-AK-5.chrom.dedup.bam \
-I ../BAM_files/PP-AK-6.chrom.dedup.bam \
-I ../BAM_files/PP-AK-7.chrom.dedup.bam \
-I ../BAM_files/PP-AK-8.chrom.dedup.bam \
-I ../BAM_files/PP-AK-9.chrom.dedup.bam \
-I ../BAM_files/PP-MW-11.chrom.dedup.bam \
-I ../BAM_files/PP-MW-18.chrom.dedup.bam \
-I ../BAM_files/PP-MW-19.chrom.dedup.bam \
-I ../BAM_files/PP-MW-20.chrom.dedup.bam \
-I ../BAM_files/PP-MW-21.chrom.dedup.bam \
-I ../BAM_files/PP-MW-22.chrom.dedup.bam \
-I ../BAM_files/PP-MW-2.chrom.dedup.bam \
-I ../BAM_files/PP-MW-3.chrom.dedup.bam \
-I ../BAM_files/PP-MW-4.chrom.dedup.bam \
-I ../BAM_files/PP-MW-7.chrom.dedup.bam \
-I ../BAM_files/PP-NE-26.chrom.dedup.bam \
-I ../BAM_files/PP-NE-37.chrom.dedup.bam \
-I ../BAM_files/PP-NE-38.chrom.dedup.bam \
-I ../BAM_files/PP-NE-39.chrom.dedup.bam \
-I ../BAM_files/PP-NE-40.chrom.dedup.bam \
-I ../BAM_files/PP-NE-42.chrom.dedup.bam \
-I ../BAM_files/PP-NE-43.chrom.dedup.bam \
-I ../BAM_files/PP-NE-47.chrom.dedup.bam \
-I ../BAM_files/PP-NE-48.chrom.dedup.bam \
-I ../BAM_files/PP-NE-49.chrom.dedup.bam \
-I ../BAM_files/PP-NR-01.chrom.dedup.bam \
-I ../BAM_files/PP-NR-02.chrom.dedup.bam \
-I ../BAM_files/PP-NR-03.chrom.dedup.bam \
-I ../BAM_files/PP-NR-04.chrom.dedup.bam \
-I ../BAM_files/PP-NR-05.chrom.dedup.bam \
-I ../BAM_files/PP-NR-06.chrom.dedup.bam \
-I ../BAM_files/PP-NR-07.chrom.dedup.bam \
-I ../BAM_files/PP-NR-08.chrom.dedup.bam \
-I ../BAM_files/PP-NR-09.chrom.dedup.bam \
-I ../BAM_files/PP-NR-10.chrom.dedup.bam \
-I ../BAM_files/PP-NW-10.chrom.dedup.bam \
-I ../BAM_files/PP-NW-11.chrom.dedup.bam \
-I ../BAM_files/PP-NW-12.chrom.dedup.bam \
-I ../BAM_files/PP-NW-13.chrom.dedup.bam \
-I ../BAM_files/PP-NW-15.chrom.dedup.bam \
-I ../BAM_files/PP-NW-16.chrom.dedup.bam \
-I ../BAM_files/PP-NW-18.chrom.dedup.bam \
-I ../BAM_files/PP-NW-5.chrom.dedup.bam \
-I ../BAM_files/PP-NW-8.chrom.dedup.bam \
-I ../BAM_files/PP-NW-9.chrom.dedup.bam \
-I ../BAM_files/PP-SE-01.chrom.dedup.bam \
-I ../BAM_files/PP-SE-02.chrom.dedup.bam \
-I ../BAM_files/PP-SE-08.chrom.dedup.bam \
-I ../BAM_files/PP-SE-09.chrom.dedup.bam \
-I ../BAM_files/PP-SE-10.chrom.dedup.bam \
-I ../BAM_files/PP-SE-12.chrom.dedup.bam \
-I ../BAM_files/PP-SE-14.chrom.dedup.bam \
-I ../BAM_files/PP-SE-15.chrom.dedup.bam \
-I ../BAM_files/PP-SE-16.chrom.dedup.bam \
-I ../BAM_files/PP-SE-18.chrom.dedup.bam \
-I ../BAM_files/PP-SR-09.chrom.dedup.bam \
-I ../BAM_files/PP-SR-10.chrom.dedup.bam \
-I ../BAM_files/PP-SR-12.chrom.dedup.bam \
-I ../BAM_files/PP-SR-13.chrom.dedup.bam \
-I ../BAM_files/PP-SR-15.chrom.dedup.bam \
-I ../BAM_files/PP-SR-16.chrom.dedup.bam \
-I ../BAM_files/PP-SR-17.chrom.dedup.bam \
-I ../BAM_files/PP-SR-18.chrom.dedup.bam \
-I ../BAM_files/PP-SR-19.chrom.dedup.bam \
-I ../BAM_files/PP-SR-21.chrom.dedup.bam \
-o all_samples.intervals
