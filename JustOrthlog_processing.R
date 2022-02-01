library(ips)
library(parallel)
library(MASS)
library(rBLAST)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:\\Program Files\\NCBI\\blast-2.11.0+\\bin", sep= .Platform$path.sep))

# Read the output from Justortholog (this assumes the file has the same number of rows across lines)
tabs <- read.table("othologs.complete.txt", sep = "\t", header = F, na.strings=c("","NA"))

#i=1
#z=2
for (i in 1:nrow(tabs)){      # looping across different ortholog groups
  print(paste(i,"out of",nrow(tabs),"ortholog groups"))
  a <- as.matrix(tabs[i,])
  
  f_tot=NULL
  for (z in 2:4) {            # looping species
    if (!is.na(a[z])){
      print(a[z])
      sp <- gsub(".fa.*", "", a[z])
      gene <- tail(unlist(strsplit(a[z],":")), n=1)
      orth <- gsub(" ", "_", a[1])
      orth <- gsub(":", "", orth)
      file <- paste0("cds//", sp,".sorted.cleaned.cds.tab.fa")
      w <- grep(gene,readLines(file))
      f <- readLines(file)[w]
      f <- strsplit(f, "\t|\\*")            # split each cds
      f1 <- cbind(sp, length(f[[1]]))
      f_tot <- as.data.frame(rbind(f_tot, f1))
      # p=2
      for (p in 2:(length(f[[1]])-1)){          # looping across cds
        n <- paste(sp,"_",gsub(";", "", f[[1]][1]), "_", p-1, "\t", f[[1]][p], sep="")
        n <- gsub("^", ">", n)
        n <- gsub("\t", "\n", n, perl = T)
        new_file <- paste("Orthologs_cds_splitted/", orth, ".fasta", sep = "")
        write(n, file=new_file ,append=TRUE)
      }
    }
  }
    
  f_tot <- f_tot[which.max(f_tot$V2),]
  print(f_tot[,1])
  
  #################### Converting fasta to tabular #########################
  
  input <- readLines(new_file)
  output_file <- gsub(".fasta", ".tab", new_file)
  output <- file(output_file,"w")
  
  currentSeq <- 0
  newLine <- 0
  
  for(i in 1:length(input)) {
    if(strtrim(input[i], 1) == ">") {
      if(currentSeq == 0) {
        writeLines(paste(input[i],"\t"), output, sep="")
        currentSeq <- currentSeq + 1
      } else {
        writeLines(paste("\n",input[i],"\t", sep =""), output, sep="")
      }
    } else {
      writeLines(paste(input[i]), output, sep="")
    }
  }
  
  close(output)
  
  #########################################################################

  # Creating database for BLAST
  where <- grep(f_tot[,1],readLines(paste("Orthologs_cds_splitted/", orth, ".tab", sep = ""),warn=F))
  sin <- readLines(paste("Orthologs_cds_splitted/", orth, ".tab", sep = ""),warn=F)[where]
  sin <- gsub("\t", "\n", sin)
  write(sin, file= paste("Orthologs_cds_splitted/", orth, ".db.fasta", sep = "") ,append=FALSE)
  
  #prepare for a BLAST query
  print(paste("Blasting",new_file,sep=" "))
  makeblastdb(paste("Orthologs_cds_splitted/", orth, ".db.fasta", sep = ""), dbtype = "nucl")
  blast1 <- blast(db=paste("Orthologs_cds_splitted/", orth, ".db.fasta", sep = ""),type="blastn")
  print(blast1,info=TRUE)
  #Run BLAST query
  query <- predict(blast1,readDNAStringSet(new_file, format='fasta'),BLAST_args="-max_target_seqs 1")
  
  print(paste("Separating cds for",new_file,sep=" "))
  for (cds in 1:length(unique(query$SubjectID))) {
    print(paste("CDS:",unique(query$SubjectID)[cds]))
    #cds=12
    for (tes in 1:nrow(query)) {
      #print(tes)
      if (identical(query$SubjectID[tes], unique(query$SubjectID)[cds])) {
        #print(paste("MATCH",query$SubjectID[tes],query$SubjectID)[cds])
        #print(as.character(blast$V1[[tes]]))
        ss <- paste0(query$QueryID[tes],"[ \t]")  # $ guarantees it's the end of the string for grep
        this <- grep(ss,readLines(output_file,warn=F))
        sequence <- readLines(output_file,warn=F)[this]
        write(sequence, file= paste("Orthologs_cds_splitted/Separated_cds/", orth, "_", cds, ".fasta", sep = "") ,append=T)
        print("I wrote a sequence")
        }
      }
    }
}
