# Code from: Demography and linked selection interact to shape the genomic landscape of codistributed woodpeckers during the Ice Age

Code used for analyses in Moreira et al. 2022:

### Read alignment, variant calling, and filtering

* `Chromosemble.sh`: maps the scaffolds of the *Picoides pubescens* genome onto the Zebra finch chromosomes using [Satsuma Chromosemble](http://satsuma.sourceforge.net/), generating a pseudochromosome reference.

* `BWA_alignment.sh`: trims adaptors from raw reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), assesses quality of reads with fastqc, aligns them to the genome reference using [BWA](http://bio-bwa.sourceforge.net/), converts SAM file to sorted BAM format, adds read group information, marks PCR/optical duplicated reads, and collects alignmnet QC metrics using [Picard](https://broadinstitute.github.io/picard/) and [QualiMap](http://qualimap.conesalab.org/).

* `RealignerTargetCreator.sh`: creates list of target intervals for indel realignment using [GATK](https://gatk.broadinstitute.org/hc/en-us).

* `Indel_Realigner.sh`: performs local realignment around indels using [GATK](https://gatk.broadinstitute.org/hc/en-us).

* `GATK_SNP_calling.sh`: calls genotypes using [GATK Haplotype Caller](https://gatk.broadinstitute.org/hc/en-us).

* `GATK_SNP_GVCF.sh`: jointly calls genotypes across all samples using [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us).

* `SNP_select_variants-25missingdata.sh`: filters SNPs from vcf file using [VCFtools](http://vcftools.sourceforge.net/).

* `ANGSD_GL.sh`: estimates genotype likelihood from BAM files using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

### Population structure

* `ANGSD_FST.sh`: estimates pairwise global and window-based FST using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `SNPRelate.R`: performs principal component analysis (PCA using the R package [SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html).

### Demographic inference

### Genetic diversity, recombination rates, and linkage disequilibrium

* `ANGSD_Theta.sh`: estimates genetic diversity (theta) using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

### Genomic predictors of regional variation in nucleotide diversity

### Natural selection and genetic load
