## Made by Lucas R. Moreira
## Last updated 12 April 2020
## Usage: Calculates gene density along specified sliding windows

size_of_genome <- read.table("pseudochromosomes-subset_sorted.fasta.size_of_genome.txt",header=F)

bed <- read.table("Picoides_pubescens.pseudogenome.cds.bed",sep="\t",header=F)
colnames(bed) <- c("chromosome","start","end","source","match","strand","genome","feature","frame","id")

# Now let's fill the table with info

whole_data <- read.csv("whole_data.csv",header=T)
chromos <- c(1,10,11,12,13,14,15,17,18,19,"1A",2,20,21,22,23,24,25,26,27,28,3,4,"4A",5,6,7,8,9,"LGE22","Z")

# create table to store data

data_gene <- numeric()

for(i in 1:length(chromos)){
  chr <- chromos[i]
  size <- size_of_genome[size_of_genome$V1==chr,2]
  print(chr)
  
  sub_data <- bed[bed$chromosome==chr,] #select data for a single chromosome
  sub_data <- sub_data[sub_data$feature=="CDS",]
  
  if(nrow(sub_data)!=0){
    # create table to represent which bp are part of a gene
    whole_seq <- 1:size
    gene_present <- rep(0,size)
    where_are_genes <- data.frame(chr,whole_seq,gene_present)
    
    print(paste0("Starting locating CDS in chromosome ",chr))
    # all bp that are part of a gene get a 1, otherwise 0
    for(feature in 1:nrow(sub_data)){
      print(paste0(feature," of ",nrow(sub_data)," CDS"))
      start_ <- sub_data$start[feature]
      end_ <- sub_data$end[feature]
      where_are_genes$gene_present[start_:end_] <- 1
    }
  }

  chrom_data <- whole_data[whole_data$chrom==chr,]
  for(i in 1:nrow(chrom_data)){
    chromosome <- as.character(chrom_data$chrom[i])
    first <- chrom_data$start[i]
    last <- chrom_data$end[i]
    
    print(paste0("Calculating CDS density for chromosome ",chromosome," at ",first,":",last))
    coding_seq_length <- sum(where_are_genes$gene_present[first:last])
    coding_seq_percent <- coding_seq_length/(last-first)
    
    # count the number of CDS features that start within the window
    sub_data <- bed[bed$chromosome==chromosome,]
    sub_data <- sub_data[sub_data$feature=="CDS",]
    a <- subset(sub_data,start>=first&start<=last)
    count_cds <- nrow(a)
  
    table <- data.frame(chromosome,start=first,end=last,coding_seq_percent,count_cds)
    data_gene <- rbind(data_gene,table)
  }
}

new_data_table <- cbind(whole_data,coding_seq_percent=data_gene$coding_seq_percent,count_cds=data_gene$count_cds)
write.csv(new_data_table,"Total_data.csv",row.names=F)
