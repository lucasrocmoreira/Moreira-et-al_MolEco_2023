## Made by Lucas R. Moreira
## Last updated 12 April 2020
## Usage: Calculates GC content along specified sliding windows

#install.packages('seqinr')
library("seqinr")

seq <- read.fasta(file = "C:\\Users\\lexlu\\OneDrive\\Biology\\AMNH\\Projects\\Woodpeckers project\\Reference Genome\\pseudochromosomes-subset_sorted.fasta")
whole_data <- read.csv("Total_data.csv",header=T)

GCs <- numeric()

for(i in 1:nrow(whole_data)){
  chrom <- as.character(whole_data$chrom[i])
  first <- whole_data$start[i]
  last <- whole_data$end[i]
  
  print(paste0("Calculating GC content for chromosome ",chrom," at ",first,":",last))
  
  chromosome <- seq[[chrom]]
  chunk <- chromosome[first:last]
  chunkGC <- try(GC(chunk),TRUE)
  if(class(chunkGC)=="try-error"){
    chunkGC <- NA
  } else{
  chunkGC <- GC(chunk)
  
  data <- data.frame(chrom,start=first,end=last,chunkGC)
  GCs <- rbind(GCs,data)
  }
}

complete_data <- cbind(whole_data,GC_content=GCs$chunkGC)

write.csv(complete_data,"Predictors_of_theta.csv",row.names=F)
