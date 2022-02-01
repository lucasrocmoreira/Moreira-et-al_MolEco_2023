## Made by Lucas R. Moreira
## Last updated 12 April 2020
## Usage: Calculates weighted recombination rates for 100kb windows

###################
### Import theta ##
###################

theta <- read.table("theta.thetasWindow.gz.pestPG",header=F)
colnames(theta) <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
# We then have 5 different estimators of theta, these are: Watterson, pairwise, FuLi, fayH, L. And we have 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data in the window.

# Now simplify theta table from ANGSD

theta_simplified <- numeric()
for(i in 1:nrow(theta)){
  windows_start <- strsplit(as.character(theta[i,1]),"[^A-z0-9_]")[[1]][5]  # also, 2
  windows_end <- strsplit(as.character(theta[i,1]),"[^A-z0-9_]")[[1]][6]    # also, 3
  x <- cbind(windows_start,windows_end)
  theta_simplified <- rbind(theta_simplified,x)
}
chromosome <- theta$Chr
theta_values <- theta[,4:13]

theta_simplified <- cbind(theta_simplified,chromosome,theta_values)
write.csv(theta_simplified,"theta_simplified.csv",row.names = F)

#read.csv("theta_simplified.csv",header=T)

theta_simplified <- read.csv("theta_simplified.csv",header=T) # after producing this file according to above, remove first column of numbers. The code below only works for csv imported data.
attach(theta_simplified)
  
###########################
### Import recombination ##
###########################

rr <- read.table("SNP-only.E.maf002.recode.PREDICT.BSCORRECTED.txt",header=T)

chromosomes <- c(1,"1A",2,3,4,"4A",5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,"LGE22","Z")
for(chr in chromosomes){
  print(chr)
  # set chromosome you are analyzing (add 1,2,3 or "1A" "Z")
  chromo <- chr
  
  sub_rr <- rr[rr$chrom==chromo,]
  size_of_chromosome <- sub_rr$end[nrow(sub_rr)]
  recombination_per_bp <- data.frame("base_pair"=0:size_of_chromosome,"chrom"=chromo,"recombination_rate"=NA)
  
  for(i in 1:nrow(sub_rr)){
    print(paste0(i," of ",nrow(sub_rr)))
    start_ <- sub_rr$start[i]
    end_ <- sub_rr$end[i]
    recombination_per_bp$recombination_rate[start_:end_] <- sub_rr$recombRate[i]
  }
  
  # Calculate weighted average of recombination rates, according to 
  
  #as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}  # we need a function to convert factors to integers
  weighted_recR <- numeric()
  for(i in 1:nrow(theta_simplified[chromosome==chromo,])){  # only the correct chromosome
    #starting <- as.numeric.factor(theta_simplified[chromosome==chromo,][i,1])
    starting <- as.numeric(theta_simplified[chromosome==chromo,][i,1])
    #ending <- as.numeric.factor(theta_simplified[chromosome==chromo,][i,2])
    ending <- as.numeric(theta_simplified[chromosome==chromo,][i,2])
    weighted_mean <- sum(na.omit(recombination_per_bp[starting:ending,3]))/length(na.omit(recombination_per_bp[starting:ending,3]))
    weighted_recTable <- cbind(starting,ending,weighted_mean)
    weighted_recR <- rbind(weighted_recR,weighted_recTable)
  }
  
  # Correlation between theta and recomb rates
  pdf(paste0(chromo,".hist.pdf"),height=4,width=6)
  hist(weighted_recR[,3])
  dev.off()
  
  model <- lm(as.numeric(theta_simplified[chromosome==chromo,4])~weighted_recR[,3])
  sink(paste0(chromo,".tW.lm.txt"))
  print(summary(model))
  sink()

  model2 <- lm(as.numeric(theta_simplified[chromosome==chromo,5])~weighted_recR[,3])
  sink(paste0(chromo,".tP.lm.txt"))
  print(summary(model))
  sink()
  
  pdf(paste0("chr",chromo,".thetaW-recombR.pdf"),height=4,width=6)
  plot(weighted_recR[,3],as.numeric(theta_simplified[chromosome==chromo,4]),xlab="Recombination rate",ylab="Theta W")
  abline(lm(as.numeric(theta_simplified[chromosome==chromo,4])~weighted_recR[,3]))
  dev.off()

  pdf(paste0("chr",chromo,".thetaP-recombR.pdf"),height=4,width=6)
  plot(weighted_recR[,3],as.numeric(theta_simplified[chromosome==chromo,5]),xlab="Recombination rate",ylab="Theta Pi")
  abline(lm(as.numeric(theta_simplified[chromosome==chromo,5])~weighted_recR[,3]))
  dev.off()
  
  chr_data <- cbind(chrom=chromo,start=weighted_recR[,1],end=weighted_recR[,2],thetaW=as.numeric(theta_simplified[chromosome==chromo,4]),thetaP=as.numeric(theta_simplified[chromosome==chromo,5]),tajimasD=as.numeric(theta_simplified[chromosome==chromo,9]),recR=weighted_recR[,3])
  write.csv(chr_data,paste0("chr",chromo,".csv"))
}

whole_data <- numeric()
for(file in list.files("Individual chromosome data/",full.names = T)){
  print(file)
  chr_data <- read.csv(file,header=T)[,-1]
  whole_data <- rbind(whole_data,chr_data)
}
write.csv(whole_data,"whole_data.csv",row.names=F)


model <- lm(whole_data$thetaP~whole_data$recR)
summary(model)
plot(whole_data$recR,whole_data$thetaP,xlab="Recombination rate (c/bp)",ylab="Theta pi")
abline(lm(whole_data$thetaP~whole_data$recR),col="red")

mean(na.omit(whole_data$recR))

mean(whole_data$thetaP/100000)
