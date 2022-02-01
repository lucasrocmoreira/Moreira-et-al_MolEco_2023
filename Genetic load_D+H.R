#####################################################################
# Script to calculate population genetic load                       #
# Written by Lucas R. Moreira
# Last updated 01.28.21
####################################################################

# load geno file

A <- function(x){strsplit(as.character(x),"")}

Downy <- 1:12
Hairy <- 13:24

#################### LOW IMPACT (synonymous) ########################

low <- read.table("SNP-only.LOW.ANN_simplified.geno",header=F,colClasses = c("factor"))
low_character <- lapply(low[,1],A)
df_low <- data.frame(matrix(unlist(low_character), nrow=length(low_character), byrow=T))
colnames(df_low) <- c("PP-AK-10","PP-AK-4","PP-AK-7","PP-MW-18","PP-NE-48","PP-NR-02","PP-NR-04","PP-NW-10","PP-NW-16","PP-NW-9","PP-SE-12","PP-SR-17","PV-AK-10","PV-AK-6","PV-AK-9","PV-MW-7","PV-NE-28","PV-NR-03","PV-NR-07","PV-NW-12","PV-NW-18","PV-NW-21","PV-SE-8","PV-SR-19")
df_low[df_low==9] <- NA

#################### MODERATE IMPACT (non-synonymous) ####################

moderate <- read.table("SNP-only.MODERATE.ANN_simplified.geno",header=F,colClasses = c("factor"))
moderate_character <- lapply(moderate[,1],A)
df_moderate <- data.frame(matrix(unlist(moderate_character), nrow=length(moderate_character), byrow=T))
colnames(df_moderate) <- c("PP-AK-10","PP-AK-4","PP-AK-7","PP-MW-18","PP-NE-48","PP-NR-02","PP-NR-04","PP-NW-10","PP-NW-16","PP-NW-9","PP-SE-12","PP-SR-17","PV-AK-10","PV-AK-6","PV-AK-9","PV-MW-7","PV-NE-28","PV-NR-03","PV-NR-07","PV-NW-12","PV-NW-18","PV-NW-21","PV-SE-8","PV-SR-19")
df_moderate[df_moderate==9] <- NA

#################### HIGH IMPACT (loss of start and stop codons) ########################

high <- read.table("SNP-only.HIGH.ANN_simplified.geno",header=F,colClasses = c("factor"))
high_character <- lapply(high[,1],A)
df_high <- data.frame(matrix(unlist(high_character), nrow=length(high_character), byrow=T))
colnames(df_high) <- c("PP-AK-10","PP-AK-4","PP-AK-7","PP-MW-18","PP-NE-48","PP-NR-02","PP-NR-04","PP-NW-10","PP-NW-16","PP-NW-9","PP-SE-12","PP-SR-17","PV-AK-10","PV-AK-6","PV-AK-9","PV-MW-7","PV-NE-28","PV-NR-03","PV-NR-07","PV-NW-12","PV-NW-18","PV-NW-21","PV-SE-8","PV-SR-19")
df_high[df_high==9] <- NA

#####################################################
############### Calculate allele freq ###############
#####################################################

species <- list(Downy, Hairy)

af <- function(x,species){
  # x=df_high
  num <- 1:nrow(x)
  allele_frequency <- numeric()
  
  # Allele frequencies
  
  for(i in species){
    
    #i = Downy
    subset.data <- x[,i]
    
    # Frequency
    freq_per_locus <- numeric()
    for(i in 1:nrow(subset.data)){
      #i=12 # for testing
      print(i)
      SNP <- subset.data[i,]
      SNP.no.missing.data <- SNP[!is.na(SNP)]
      n <- 2*length(SNP.no.missing.data)
      freq_alternative <- sum(as.integer(SNP.no.missing.data))/n #each individual is diploid
      freq_per_locus <- c(freq_per_locus,freq_alternative)
    }
    
    allele_frequency <- cbind(allele_frequency,freq_per_locus)
  }
  
  allele_frequency <- cbind(num,allele_frequency)
  
  return(allele_frequency)
}

allele.freq.low <- af(x=df_low,species=species)
allele.freq.moderate <- af(x=df_moderate,species=species)
allele.freq.high <- af(x=df_high,species=species)

######### Select only fixed SNPs (where ancestral state is unambiguous) #########

############ HIGH ############
fixed_D.high <- which(allele.freq.high[,2]==0|allele.freq.high[,2]==1)
fixed_H.high <- which(allele.freq.high[,3]==0|allele.freq.high[,3]==1,)
any_fixed.high <- union(fixed_D.high,fixed_H.high)

a <- rowSums(sapply(df_high[fixed_H.high,Hairy],function(x) as.numeric(as.character(x))),na.rm=T)

change_ancestral_D.high <- which(rowSums(sapply(df_high[fixed_H.high,Hairy],function(x) as.numeric(as.character(x))),na.rm=T)>0)
change_ancestral_H.high <- which(rowSums(sapply(df_high[fixed_D.high,Downy],function(x) as.numeric(as.character(x))),na.rm=T)>0)

############ LOW ############
fixed_D.low <- which(allele.freq.low[,2]==0|allele.freq.low[,2]==1)
fixed_H.low <- which(allele.freq.low[,3]==0|allele.freq.low[,3]==1,)
any_fixed.low <- union(fixed_D.low,fixed_H.low)

change_ancestral_D.low <- which(rowSums(sapply(df_low[fixed_H.low,Hairy],function(x) as.numeric(as.character(x))),na.rm=T)>0)
change_ancestral_H.low <- which(rowSums(sapply(df_low[fixed_D.low,Downy],function(x) as.numeric(as.character(x))),na.rm=T)>0)

############ MODERATE ############
fixed_D.moderate <- which(allele.freq.moderate[,2]==0|allele.freq.moderate[,2]==1)
fixed_H.moderate <- which(allele.freq.moderate[,3]==0|allele.freq.moderate[,3]==1,)
any_fixed.moderate <- union(fixed_D.moderate,fixed_H.moderate)

change_ancestral_D.moderate <- which(rowSums(sapply(df_moderate[fixed_H.moderate,Hairy],function(x) as.numeric(as.character(x))),na.rm=T)>0)
change_ancestral_H.moderate <- which(rowSums(sapply(df_moderate[fixed_D.moderate,Downy],function(x) as.numeric(as.character(x))),na.rm=T)>0)

#################################
######### Polarize SNPs #########
#################################

# These are the data sets of SNPs that are eligible (they have ancestral state info)
SNPS_Downy.high <- df_high[any_fixed.high,Downy]
SNPS_Downy.moderate <- df_moderate[any_fixed.moderate,Downy]
SNPS_Downy.low <- df_low[any_fixed.low,Downy]
SNPS_Hairy.high <- df_high[any_fixed.high,Hairy]
SNPS_Hairy.moderate <- df_moderate[any_fixed.moderate,Hairy]
SNPS_Hairy.low <- df_low[any_fixed.low,Hairy]

polarize <- function(x,change){
  genetic.data <- x
  
  for(row in change){
    print(paste(row))
    genetic.data[row,][genetic.data[row,]==2] <- 44       #change to 44 and 55 first because if I change to 2 and 0 they get mixed up
    genetic.data[row,][genetic.data[row,]==0] <- 55
  }
  
  genetic.data[genetic.data==44] <- 0
  genetic.data[genetic.data==55] <- 2
  return(genetic.data)
}

PP.high.variants <- polarize(x=SNPS_Downy.high,change=change_ancestral_D.high)
PP.moderate.variants <- polarize(x=SNPS_Downy.moderate,change=change_ancestral_D.moderate)
PP.low.variants <- polarize(x=SNPS_Downy.low,change=change_ancestral_D.low)
PV.high.variants <- polarize(x=SNPS_Hairy.high,change=change_ancestral_H.high)
PV.moderate.variants <- polarize(x=SNPS_Hairy.moderate,change=change_ancestral_H.moderate)
PV.low.variants <- polarize(x=SNPS_Hairy.low,change=change_ancestral_H.low)
# write.csv(PP.high.variants,"PP.high.variants.csv",row.names=FALSE)
# write.csv(PP.moderate.variants,"PP.moderate.variants.csv",row.names=FALSE)
# write.csv(PP.low.variants,"PP.low.variants.csv",row.names=FALSE)
# write.csv(PV.high.variants,"PV.high.variants.csv",row.names=FALSE)
# write.csv(PV.moderate.variants,"PV.moderate.variants.csv",row.names=FALSE)
# write.csv(PV.low.variants,"PV.low.variants.csv",row.names=FALSE)
PP.high.variants <- read.csv("PP.high.variants.csv",header = T)
PP.moderate.variants <- read.csv("PP.moderate.variants.csv",header = T)
PP.low.variants <- read.csv("PP.low.variants.csv",header = T)
PV.high.variants <- read.csv("PV.high.variants.csv",header = T)
PV.moderate.variants <- read.csv("PV.moderate.variants.csv",header = T)
PV.low.variants <- read.csv("PV.low.variants.csv",header = T)

################# COUNTING #######################

# Missing data
PP.high_miss <- colSums(is.na(PP.high.variants))
PV.high_miss <- colSums(is.na(PV.high.variants))
PP.moderate_miss <- colSums(is.na(PP.moderate.variants))
PV.moderate_miss <- colSums(is.na(PV.moderate.variants))
PP.low_miss <- colSums(is.na(PP.low.variants))
PV.low_miss <- colSums(is.na(PV.low.variants))

# Additive model

PP.high <- do.call(cbind, lapply(PP.high.variants, as.numeric))
PP.high_add <- colSums(PP.high,na.rm=T)
PV.high <- do.call(cbind, lapply(PV.high.variants, as.numeric))
PV.high_add <- colSums(PV.high,na.rm=T)
High_additive <- data.frame(rbind(cbind(Count=PP.high_add,Species="Downy",
                                        Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                  Missing=PP.high_miss),
                                  cbind(Count=PV.high_add,Species="Hairy",
                                        Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                  Missing=PV.high_miss)))
High_additive$Count.ratio <- ((as.numeric(High_additive$Count))/(nrow(PP.high)-as.numeric(High_additive$Missing))*nrow(PP.high))


PP.moderate <- do.call(cbind, lapply(PP.moderate.variants, as.numeric))
PP.moderate_add <- colSums(PP.moderate,na.rm=T)
PV.moderate <- do.call(cbind, lapply(PV.moderate.variants, as.numeric))
PV.moderate_add <- colSums(PV.moderate,na.rm=T)
Moderate_additive <- data.frame(rbind(cbind(Count=PP.moderate_add,Species="Downy",
                                        Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                      Missing=PP.moderate_miss),
                                  cbind(Count=PV.moderate_add,Species="Hairy",
                                        Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                      Missing=PV.moderate_miss)))
Moderate_additive$Count.ratio <- ((as.numeric(Moderate_additive$Count))/(nrow(PP.moderate)-as.numeric(Moderate_additive$Missing))*nrow(PP.moderate))


PP.low <- do.call(cbind, lapply(PP.low.variants, as.numeric))
PP.low_add <- colSums(PP.low,na.rm=T)
PV.low <- do.call(cbind, lapply(PV.low.variants, as.numeric))
PV.low_add <- colSums(PV.low,na.rm=T)
Low_additive <- data.frame(rbind(cbind(Count=PP.low_add,Species="Downy",
                                            Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                     Missing=PP.low_miss),
                                      cbind(Count=PV.low_add,Species="Hairy",
                                            Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                      Missing=PV.low_miss)))
Low_additive$Count.ratio <- ((as.numeric(Low_additive$Count))/(nrow(PP.low)-as.numeric(Low_additive$Missing))*nrow(PP.low))


# Recessive model

PP.high_rec <- colSums(PP.high==2,na.rm=T)
PP.high_rec_miss <- colSums(is.na(PP.high))
PV.high_rec <- colSums(PV.high==2,na.rm=T)
PV.high_rec_miss <- colSums(is.na(PV.high))
High_recessive <- data.frame(rbind(cbind(Count=PP.high_rec,Species="Downy",
                                            Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                         Missing=PP.high_miss),
                                      cbind(Count=PV.high_rec,Species="Hairy",
                                            Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                            Missing=PV.high_miss)))
High_recessive$Count.ratio <- ((as.numeric(High_recessive$Count))/(nrow(PP.high)-as.numeric(High_recessive$Missing))*nrow(PP.high))

PP.moderate_rec <- colSums(PP.moderate==2,na.rm=T)
PV.moderate_rec <- colSums(PV.moderate==2,na.rm=T)
Moderate_recessive <- data.frame(rbind(cbind(Count=PP.moderate_rec,Species="Downy",
                                         Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                         Missing=PP.moderate_miss),
                                   cbind(Count=PV.moderate_rec,Species="Hairy",
                                         Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                         Missing=PV.moderate_miss)))
Moderate_recessive$Count.ratio <- ((as.numeric(Moderate_recessive$Count))/(nrow(PP.moderate)-as.numeric(Moderate_recessive$Missing))*nrow(PP.moderate))


PP.low_rec <- colSums(PP.low==2,na.rm=T)
PV.low_rec <- colSums(PV.low==2,na.rm=T)
Low_recessive <- data.frame(rbind(cbind(Count=PP.low_rec,Species="Downy",
                                         Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                        Missing=PP.low_miss),
                                   cbind(Count=PV.low_rec,Species="Hairy",
                                         Pop=c("AK","AK","AK","E","E","R","R","NW","NW","NW","E","R"),
                                         Missing=PV.low_miss)))
Low_recessive$Count.ratio <- ((as.numeric(Low_recessive$Count))/(nrow(PP.low)-as.numeric(Low_recessive$Missing))*nrow(PP.low))

##########################################################

write.csv(High_additive,"High_additive.csv",row.names=FALSE)
write.csv(Moderate_additive,"Moderate_additive.csv",row.names=FALSE)
write.csv(Low_additive,"Low_additive.csv",row.names=FALSE)
write.csv(High_recessive,"High_recessive.csv",row.names=FALSE)
write.csv(Moderate_recessive,"Moderate_recessive.csv",row.names=FALSE)
write.csv(Low_recessive,"Low_recessive.csv",row.names=FALSE)

High_additive <- read.csv("High_additive.csv",header = T)
Moderate_additive <- read.csv("Moderate_additive.csv",header = T)
Low_additive <- read.csv("Low_additive.csv",header = T)
High_recessive <- read.csv("High_recessive.csv",header = T)
Moderate_additive <- read.csv("Moderate_recessive.csv",header = T)
Low_recessive <- read.csv("Low_recessive.csv",header = T)

# Let's calculate the ration high/low

High.Low.RATIO <- High_recessive$Count.ratio/Low_recessive$Count.ratio
High.Low.RATIO <- data.frame(Species=High_recessive$Species,Pop=High_recessive$Pop,H.L.recessive.Ratio=High.Low.RATIO)
write.csv(High.Low.RATIO,"High.Low.RATIO.csv",row.names=FALSE)

pdf("High.over.Low.recessive.pdf",width=5,height =5)
p <- ggplot(data=High.Low.RATIO, aes(x=Pop, y=H.L.recessive.Ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Recessive load (High impact / Low impact)")
dev.off()

# Downy vs Hairy
kruskal.test(H.L.recessive.Ratio ~ Species, data = High.Low.RATIO)

kruskal.test(H.L.recessive.Ratio ~ Species, data = High.Low.RATIO[High.Low.RATIO$Pop=="AK",])
kruskal.test(H.L.recessive.Ratio ~ Species, data = High.Low.RATIO[High.Low.RATIO$Pop=="E",])
kruskal.test(H.L.recessive.Ratio ~ Species, data = High.Low.RATIO[High.Low.RATIO$Pop=="NW",])
kruskal.test(H.L.recessive.Ratio ~ Species, data = High.Low.RATIO[High.Low.RATIO$Pop=="R",])

boxplot(H.L.recessive.Ratio ~ Species, data = High.Low.RATIO)

High.Low.RATIO.add <- High_additive$Count.ratio/Low_additive$Count.ratio
High.Low.RATIO.add <- data.frame(Species=High_additive$Species,Pop=High_additive$Pop,H.L.additive.Ratio=High.Low.RATIO.add)
write.csv(High.Low.RATIO.add,"High.Low.RATIO.add.csv",row.names=FALSE)

pdf("High.over.Low.additive.pdf",width=5,height =5)
p <- ggplot(data=High.Low.RATIO.add, aes(x=Pop, y=H.L.additive.Ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Recessive load (High impact / Low impact)")
dev.off()

######################## PLOTTING ########################

library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

pdf("High_additive.pdf",width=5,height =5.5)
p <- ggplot(data=High_additive, aes(x=Pop, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Additive load")
dev.off()

pdf("Moderate_additive.pdf",width=5,height =6)
p <- ggplot(data=Moderate_additive, aes(x=Pop, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Additive load")
dev.off()

pdf("Low_additive.pdf",width=5,height =6)
p <- ggplot(data=Low_additive, aes(x=Pop, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Additive load")
dev.off()

pdf("High_recessive.pdf",width=5,height =5.5)
p <- ggplot(data=High_recessive, aes(x=Pop, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Recessive load")
dev.off()

pdf("Moderate_recessive.pdf",width=5,height =6)
p <- ggplot(data=Moderate_recessive, aes(x=Pop, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Recessive load")
dev.off()

pdf("Low_recessive.pdf",width=5,height =6)
p <- ggplot(data=Low_recessive, aes(x=Pop, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Recessive load")
dev.off()

# Species

pdf("High_additive-species.pdf",width=5,height =5.5)
p <- ggplot(data=High_additive, aes(x=Species, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Species") + labs(y = "Additive load")
dev.off()

pdf("Moderate_additive-species.pdf",width=5,height =6)
p <- ggplot(data=Moderate_additive, aes(x=Species, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Species") + labs(y = "Additive load")
dev.off()

pdf("Low_additive-species.pdf",width=5,height =6)
p <- ggplot(data=Low_additive, aes(x=Species, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Species") + labs(y = "Additive load")
dev.off()

pdf("High_recessive-species.pdf",width=5,height =5.5)
p <- ggplot(data=High_recessive, aes(x=Species, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Species") + labs(y = "Recessive load")
dev.off()

pdf("Moderate_recessive-species.pdf",width=5,height =6)
p <- ggplot(data=Moderate_recessive, aes(x=Species, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Species") + labs(y = "Recessive load")
dev.off()

pdf("Low_recessive-species.pdf",width=5,height =6)
p <- ggplot(data=Low_recessive, aes(x=Species, y=Count.ratio, fill=Species)) + 
  geom_boxplot()
p + labs(x = "Species") + labs(y = "Recessive load")
dev.off()

###################################################################################
###################################################################################
###################################################################################

PP.HIGH <- read.csv("PP.high.variants.csv",header=T)
PP.MODERATE <-read.csv("PP.moderate.variants.csv",header=T)
PP.LOW <-read.csv("PP.low.variants.csv",header=T)
PV.HIGH <-read.csv("PV.high.variants.csv",header=T)
PV.MODERATE <-read.csv("PV.moderate.variants.csv",header=T)
PV.LOW <-read.csv("PV.low.variants.csv",header=T)

af2 <- function(x){

  # Allele frequencies
    
  freq_per_locus <- numeric()
  for (i in 1:nrow(x)) {
    #i=1 # for testing
    print(i)
    SNP <- x[i, ]
    SNP.no.missing.data <- SNP[!is.na(SNP)]
    n <- 2 * length(SNP.no.missing.data)
    freq_alternative <-
      sum(as.integer(SNP.no.missing.data)) / n #each individual is diploid
    freq_per_locus <- c(freq_per_locus, freq_alternative)
  }
    
  allele_frequency <- data.frame(freq_per_locus=freq_per_locus)
  
  allele_frequency <- allele_frequency[!is.na(allele_frequency),]
  
  return(allele_frequency)
}

PP.af.low <- af2(x=PP.LOW)
PP.af.moderate <- af2(x=PP.MODERATE)
PP.af.high <- af2(x=PP.HIGH)
PV.af.low <- af2(x=PV.LOW)
PV.af.moderate <- af2(x=PV.MODERATE)
PV.af.high <- af2(x=PV.HIGH)

# Dataset for plotting

PP.af <- rbind(data.frame(freq_per_locus=PP.af.low,type="Low",species="Downy"),
               data.frame(freq_per_locus=PP.af.moderate,type="Moderate",species="Downy"),
               data.frame(freq_per_locus=PP.af.high,type="High",species="Downy"))
PV.af <- rbind(data.frame(freq_per_locus=PV.af.low,type="Low",species="Hairy"),
               data.frame(freq_per_locus=PV.af.moderate,type="Moderate",species="Hairy"),
               data.frame(freq_per_locus=PV.af.high,type="High",species="Hairy"))
allele_freq <- rbind(PP.af,PV.af)
allele_freq <- allele_freq[allele_freq$freq_per_locus!=0&allele_freq$freq_per_locus!=1,]
allele_freq$type <- factor(allele_freq$type,levels = c("High", "Moderate", "Low"))

write.csv(allele_freq,"Allele_freq_type.csv",row.names = F)
allele_freq <- read.csv("Allele_freq_type.csv",header = T)

pdf("Freq_Boxplot.pdf",width=5,height =5.5)
p <- ggplot(data=allele_freq, aes(x=type, y=freq_per_locus, fill=species)) + 
  geom_boxplot()
p + labs(x = "Population") + labs(y = "Additive load")
dev.off()

# Interleaved histograms
af.pp <- allele_freq[allele_freq$species=="Downy",]
af.pv <- allele_freq[allele_freq$species=="Hairy",]

pdf("PP.SFS_type.pdf",width=6,height=4)
g <- ggplot(af.pp, aes(x=freq_per_locus, color=type, fill=type)) +
  geom_histogram(position=position_dodge(0.08),bins = 10,aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                                                    ..count..[..group..==2]/sum(..count..[..group..==2]),
                                                    ..count..[..group..==3]/sum(..count..[..group..==3]))*100)) 
g + labs(x = "Frequency", y = "Percent") + theme(axis.title = element_text(size = 16),
                                                 axis.text = element_text(size = 14),
                                                 legend.text = element_text(size = 14),
                                                 legend.title = element_text(size = 14)) + ylim(0,90)
dev.off()

pdf("PV.SFS_type.pdf",width=6,height=4)
g <- ggplot(af.pv, aes(x=freq_per_locus, color=type, fill=type)) +
  geom_histogram(position=position_dodge(0.08),bins = 10,aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                                                                 ..count..[..group..==2]/sum(..count..[..group..==2]),
                                                                 ..count..[..group..==3]/sum(..count..[..group..==3]))*100)) 
g + labs(x = "Frequency", y = "Percent") + theme(axis.title = element_text(size = 16),
                                                 axis.text = element_text(size = 14),
                                                 legend.text = element_text(size = 14),
                                                 legend.title = element_text(size = 14)) + ylim(0,90)
dev.off()

# Are there differences in the freq of HIGH, MODERATE, and LOW SNPs in D & H?

# Low
kruskal.test(freq_per_locus ~ species, data = allele_freq[allele_freq$type=="Low",])

# Moderate
kruskal.test(freq_per_locus ~ species, data = allele_freq[allele_freq$type=="Moderate",])

# High
kruskal.test(freq_per_locus ~ species, data = allele_freq[allele_freq$type=="High",])


pdf("SFS.low_species.pdf",width=6,height=4)
g <- ggplot(allele_freq[allele_freq$type=="Low",], aes(x=freq_per_locus, color=species, fill=species)) +
  geom_histogram(position=position_dodge(0.08),bins = 10,aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                                                                 ..count..[..group..==2]/sum(..count..[..group..==2]))*100)) 
g + labs(x = "Frequency", y = "Percent") + theme(axis.title = element_text(size = 16),
                                                 axis.text = element_text(size = 14),
                                                 legend.text = element_text(size = 14),
                                                 legend.title = element_text(size = 14))
dev.off()

pdf("SFS.moderate_species.pdf",width=6,height=4)
g <- ggplot(allele_freq[allele_freq$type=="Moderate",], aes(x=freq_per_locus, color=species, fill=species)) +
  geom_histogram(position=position_dodge(0.08),bins = 10,aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                                                                 ..count..[..group..==2]/sum(..count..[..group..==2]))*100)) 
g + labs(x = "Frequency", y = "Percent") + theme(axis.title = element_text(size = 16),
                                                 axis.text = element_text(size = 14),
                                                 legend.text = element_text(size = 14),
                                                 legend.title = element_text(size = 14))
dev.off()

pdf("SFS.high_species.pdf",width=6,height=4)
g <- ggplot(allele_freq[allele_freq$type=="High",], aes(x=freq_per_locus, color=species, fill=species)) +
  geom_histogram(position=position_dodge(0.08),bins = 10,aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                                                                 ..count..[..group..==2]/sum(..count..[..group..==2]))*100)) 
g + labs(x = "Frequency", y = "Percent") + theme(axis.title = element_text(size = 16),
                                                 axis.text = element_text(size = 14),
                                                 legend.text = element_text(size = 14),
                                                 legend.title = element_text(size = 14))
dev.off()
