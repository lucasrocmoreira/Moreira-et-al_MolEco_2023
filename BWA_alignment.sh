#!/bin/bash
#PBS -l select=1:ncpus=32
#PBS -l walltime=5000:00:00
#PBS -N BWA_alignment
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
## fastq=fastq file read 1

ref=/home/lmoreira/reference/Picoides_pubescens_ref-pseudochromosome.v2/pseudochromosomes-subset_sorted.fasta

module load bwa-0.7.15
module load fastqc-0.11.5
module load R-3.4.1
source activate stampy

read1=$fastq
name=`echo $read1 | cut -d '.' -f1`
read2ending=".reverse_paired.fq.gz"
read2=$name$read2ending
file1=`echo $read1 | cut -d '.' -f1`
file2=`echo $read2 | cut -d '.' -f1`
header=`zcat $read1 | head -1`
IFS=':' read -a header <<< "$header"
INSTRUMENT=${header[0]}
RUN_ID=${header[1]}
FLOWCELL_BARCODE=${header[2]}
LANE=${header[3]}
ID=$FLOWCELL_BARCODE.$LANE
PU=$FLOWCELL_BARCODE.$LANE.$name
SM=$name
PL=ILLUMINA
LB=TrueSeq

echo $name $read1 $read2 $file1 $file2 $ID $PU $SM $PL $LB

echo
echo "#######################"
echo $name
echo "#######################"

echo
echo
echo Running FASTQ SCREEN
echo
echo

fastq_screen --nohits --aligner bwa $read1
fastq_screen --nohits --aligner bwa $read2

rm *.tagged.fastq.gz

echo
echo
echo Separating unpaired
echo
echo

gunzip $name.forward_paired.tagged_filter.fastq.gz $name.reverse_paired.tagged_filter.fastq.gz 
fastq_pair $name.forward_paired.tagged_filter.fastq $name.reverse_paired.tagged_filter.fastq 

echo
echo
echo Trimming adapters
echo
echo

java -jar /home/lmoreira/programs/Trimmomatic-0.36/trimmomatic-0.36.jar \
PE $file1.tagged_filter.fastq.gz $file2.tagged_filter.fastq.gz \
$name.forward_paired.fq.gz \
$name.forward_unpaired.fq.gz \
$name.reverse_paired.fq.gz \
$name.reverse_unpaired.fq.gz \
ILLUMINACLIP:/home/lmoreira/programs/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:true

echo
echo	
echo Quality assessment
echo
echo

fastqc $name.forward_paired.fq.gz
fastqc $name.reverse_paired.fq.gz

echo
echo	
echo Running BWA alignment
echo
echo

if [ ! -d "BAM_files/" ]; then mkdir BAM_files; fi 
time bwa mem \
-M \
-t 32 \
$ref \
$name.forward_paired.tagged_filter.fastq.paired.fq \
$name.reverse_paired.tagged_filter.fastq.paired.fq \
> BAM_files/$name.sam

echo	
echo
echo Converting SAM to sorted BAM
echo
echo	

mkdir tmp
time java -jar ~/programs/picard/picard.jar SortSam \
INPUT=BAM_files/$name.sam \
OUTPUT=BAM_files/$name.bam \
SORT_ORDER=coordinate \
TMP_DIR=`pwd`/tmp
rm -r tmp/
rm BAM_files/$name.sam

echo
echo
echo Adding Read Group
echo
echo

time java -jar ~/programs/picard/picard.jar AddOrReplaceReadGroups \
I=BAM_files/$name.bam \
O=BAM_files/$name.groups_added.bam \
RGID=$ID \
RGLB=$LB \
RGPL=$PL \
RGPU=$PU \
RGSM=$SM
rm BAM_files/$name.bam

echo
echo
echo Deduplication
echo
echo

time java -jar ~/programs/picard/picard.jar MarkDuplicates \
TMP_DIR=tmp \
INPUT=BAM_files/$name.groups_added.bam \
OUTPUT=BAM_files/$name.dedup.bam \
METRICS_FILE=Stats/$name.dedup.metrics.txt \
REMOVE_DUPLICATES=false \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024 \
TAGGING_POLICY=All
rm -r tmp/
rm BAM_files/$name.groups_added.bam

echo
echo
echo Collecting alignment metrics
echo
echo

if [ ! -d "Stats/" ]; then mkdir Stats; fi

java -jar ~/programs/picard/picard.jar CollectAlignmentSummaryMetrics \
INPUT=BAM_files/$name.dedup.bam \
OUTPUT=Stats/$name.alignment_metrics.txt \
R=$ref \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \

echo
echo
echo Collect Insert Size Metrics
echo
echo

time java -jar ~/programs/picard/picard.jar CollectInsertSizeMetrics \
INPUT=BAM_files/$name.dedup.bam \
OUTPUT=Stats/$name.insert_metrics.txt \
HISTOGRAM_FILE=Stats/$name.insert_size_histogram.pdf

echo
echo
echo Collect Coverage Metrics
echo
echo

java -jar ~/programs/picard/picard.jar CollectRawWgsMetrics \
I=BAM_files/$name.dedup.bam \
O=Stats/$name.raw_wgs_metrics.txt \
R=$ref \
INCLUDE_BQ_HISTOGRAM=true

echo
echo
echo Indexing BAM
echo
echo

time java -jar ~/programs/picard/picard.jar BuildBamIndex \
I=BAM_files/$name.dedup.bam

echo
echo
echo Running QUALIMAP
echo
echo

qualimap bamqc -nt 8 -bam BAM_files/$name.dedup.bam -outdir Stats/$name.QUALIMAP --java-mem-size=8G 

echo
echo "#######################"
echo $name DONE!
echo "#######################"
echo
echo

