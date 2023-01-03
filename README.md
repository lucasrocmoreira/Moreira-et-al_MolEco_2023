# Code from: Demography and linked selection interact to shape the genomic landscape of codistributed woodpeckers during the Ice Age

Code used for analyses in Moreira et al. 2023:

### Read alignment, variant calling, and filtering

* `Chromosemble.sh`: maps the scaffolds of the *Dryobates pubescens* genome onto the Zebra finch chromosomes using [Satsuma Chromosemble](http://satsuma.sourceforge.net/), generating a pseudochromosome reference.

* `BWA_alignment.sh`: trims adaptors from raw reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), assesses quality of reads with [FastaQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), aligns them to the genome reference using [BWA](http://bio-bwa.sourceforge.net/), converts SAM file to sorted BAM format, adds read group information, marks PCR/optical duplicated reads, and collects alignmnet QC metrics using [Picard](https://broadinstitute.github.io/picard/) and [QualiMap](http://qualimap.conesalab.org/).

* `RealignerTargetCreator.sh`: creates list of target intervals for indel realignment using [GATK](https://gatk.broadinstitute.org/hc/en-us).

* `Indel_Realigner.sh`: performs local realignment around indels using [GATK](https://gatk.broadinstitute.org/hc/en-us).

* `GATK_SNP_calling.sh`: calls genotypes using [GATK Haplotype Caller](https://gatk.broadinstitute.org/hc/en-us).

* `GATK_SNP_GVCF.sh`: jointly calls genotypes across all samples using [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us).

* `SNP_select_variants-25missingdata.sh`: filters SNPs from vcf file using [VCFtools](http://vcftools.sourceforge.net/).

* `ANGSD_GL.sh`: estimates genotype likelihood from BAM files using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

### Population structure

* `ANGSD_FST.sh`: estimates pairwise global and window-based FST using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `SNPRelate.R`: performs principal component analysis (PCA) using the R package [SNPRelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html).

* `IQ-TREE_POMO.sh`: builds a maximum likelihood tree based on the polymorphism-aware phylogenetic model (PoMo) implemented in [IQ-Tree 2](http://www.iqtree.org/).

* `NGSadmix.sh`: estimates admixture proportions with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix).

* `NGSadmix_plot.R`: plots admixture proportions.

* `EEMS.sh`: calculates estimated effective migration surface ([EEMS](https://github.com/dipetkov/eems)). **Note**: the folder `EEMS_par_file` contains the parameter file used.

### Demographic inference

* `Fastsimcoal.sh`: runs demographic models with [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal27/). **Note**: parameter files are prodived in folder `fastsimcoal2`.

* `Stairwayplot2.sh`: performs demographic inference using [Stairway Plot 2](https://github.com/xiaoming-liu/stairway-plot-v2). **Note**: blueprint files provided in folder `stairwayplot2_bluprints`.

* `Stairwayplot2_plot.R`: plots results from [Stairway Plot 2](https://github.com/xiaoming-liu/stairway-plot-v2).

### Genetic diversity, recombination rates, and linkage disequilibrium

* `ANGSD_Theta.sh`: estimates genetic diversity (theta) using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

* `ReLERNN.sh`: estimates per-base recombination rates using [ReLERNN](https://github.com/kr-colab/ReLERNN).

* `PopLDdecay.sh`: estimates linkage decay with [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay).

### Genomic predictors of regional variation in nucleotide diversity

* `GC_content.R`: computes GC content along the genome using the R package [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html).

* `Gene_density.R`: computes gene density along the genome.

* `Recombination-vs-theta.R`: computes weighted recombination rates along windows matching theta estimates.

* `Genomic_predictors.R`: statistical models for predictors of genomic diversity.

### Natural selection and genetic load

* `Genetic load.R`: polarizes SNPs in .geno format and computes genetic load metrics.

* `JustOrthlog_processing.R`: processes outputs from [JustOrthologs](https://github.com/ridgelab/JustOrthologs).

* `PAML_codeml_dnds.py`: estimates dN/dS ratio along the branches of a tree using [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html).
