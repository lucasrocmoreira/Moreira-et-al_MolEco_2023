# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install(c("gdsfmt","SNPRelate"))

library(gdsfmt)
library(SNPRelate)

snpgdsVCF2GDS("SNP-only.vcf", "Picoides_pubescens-raw.gds", method="biallelic.only")

snpgdsSummary("Picoides_pubescens-raw.gds")

genofile <- snpgdsOpen("Picoides_pubescens-raw.gds")

samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

AK <- c("PP-AK-1","PP-AK-10","PP-AK-2","PP-AK-3","PP-AK-4","PP-AK-5","PP-AK-6","PP-AK-7","PP-AK-8","PP-AK-9")
MW <- c("PP-MW-11","PP-MW-18","PP-MW-19","PP-MW-2","PP-MW-20","PP-MW-21","PP-MW-22","PP-MW-3","PP-MW-4","PP-MW-7")
NE <- c("PP-NE-26","PP-NE-37","PP-NE-38","PP-NE-39","PP-NE-40","PP-NE-42","PP-NE-43","PP-NE-47","PP-NE-48","PP-NE-49")
NR <- c("PP-NR-01","PP-NR-02","PP-NR-03","PP-NR-04","PP-NR-05","PP-NR-06","PP-NR-07","PP-NR-08","PP-NR-09","PP-NR-10")
NW <- c("PP-NW-10","PP-NW-11","PP-NW-12","PP-NW-13","PP-NW-15","PP-NW-16","PP-NW-18","PP-NW-5","PP-NW-8","PP-NW-9")
SE <- c("PP-SE-01","PP-SE-02","PP-SE-08","PP-SE-09","PP-SE-10","PP-SE-12","PP-SE-14","PP-SE-15","PP-SE-16","PP-SE-18")
SR <- c("PP-SR-09","PP-SR-10","PP-SR-12","PP-SR-13","PP-SR-15","PP-SR-16","PP-SR-17","PP-SR-18","PP-SR-19","PP-SR-21")

pop <- numeric()
for(i in samples){
  if(i %in% AK){
    pop <- c(pop,"Alaska")
  }
  if(i %in% MW){
    pop <- c(pop,"Midwest")
  }
  if(i %in% NE){
    pop <- c(pop,"Northeast")
  }
  if(i %in% NR){
    pop <- c(pop,"Northern Rockies")
  }
  if(i %in% NW){
    pop <- c(pop,"Northwest")
  }
  if(i %in% SE){
    pop <- c(pop,"Southeast")
  }
  if(i %in% SR){
    pop <- c(pop,"Souther Rockies")
  }
}

color <- numeric()
for(i in samples){
  if(i %in% AK){
    color <- c(color,"blue")
  }
  if(i %in% MW){
    color <- c(color,"green")
  }
  if(i %in% NE){
    color <- c(color,"purple")
  }
  if(i %in% NR){
    color <- c(color,"orange")
  }
  if(i %in% NW){
    color <- c(color,"black")
  }
  if(i %in% SE){
    color <- c(color,"pink")
  }
  if(i %in% SR){
    color <- c(color,"red")
  }
}

#=============================================================================================================#
# LD-based SNP pruning
#   Recursively removes SNPs within a sliding window
#=============================================================================================================#

#SNP set per amount of missing data
snpset_25 <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold=0.2, maf = 0.03, missing.rate = 0.25) #missing.rate: to use the SNPs with "<= missing.rate"
snpset_20 <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold=0.2, maf = 0.03, missing.rate = 0.20) #missing.rate: to use the SNPs with "<= missing.rate"
snpset_15 <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold=0.2, maf = 0.03, missing.rate = 0.15) #missing.rate: to use the SNPs with "<= missing.rate"
snpset_10 <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold=0.2, maf = 0.03, missing.rate = 0.10) #missing.rate: to use the SNPs with "<= missing.rate"

# Get all selected snp id

#snpset.id <- unlist(snpset_10)
#saveRDS(snpset_25,"PP-snpset_25.rds")

snpset.id <- unlist(snpset_25)
snpset.id <- unlist(readRDS("PP-snpset_25.rds"))

#=============================================================================================================#
# Principal Component Analysis (PCA)
#   Recursively removes SNPs within a sliding window
#=============================================================================================================#

pca <- snpgdsPCA(genofile,snp.id=snpset.id,num.thread=4,autosome.only=FALSE)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

write.csv(pc.percent,"pc.percent.csv")

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  stringsAsFactors = FALSE)
head(tab)

# Draw
pdf('PC1vsPC2.pdf',width=8,height=5, bg='transparent')
plot(tab$EV1, tab$EV2, xlab="PC1", ylab="PC2", pch=16, col=color)
legend("bottomright",
       col= c("blue","black","green","purple","pink","red","orange"),
       bg="white", pch=c(20,20),yjust=0,
       legend = c("Alaska","Northwest","Midwest","Northeast","Southeast","Souther Rockies","Northern Rockies"),
       cex = 0.6)
dev.off()

#PC1 (as well as PC2) separates Alaska, Eastern+Northwest and Rockies

pdf('PC1vsPC3.pdf',width=8,height=5, bg='transparent')
plot(tab$EV1, tab$EV3, xlab="PC1", ylab="PC3", pch=16, col=color)
legend("bottomright",
       col= c("blue","black","green","purple","pink","red","orange"),
       bg="white", pch=c(20,20),yjust=0,
       legend = c("Alaska","Northwest","Midwest","Northeast","Southeast","Souther Rockies","Northern Rockies"),
       cex = 0.5)
dev.off()

#PC3 separates Northwest from the rest and samples within Southeast

pdf('PC2vsPC3.pdf',width=8,height=5, bg='transparent')
plot(tab$EV2, tab$EV3, xlab="PC2", ylab="PC3", pch=16, col=color)
legend("bottomleft",
       col= c("blue","black","green","purple","pink","red","orange"),
       bg="white", pch=c(20,20),yjust=0,
       legend = c("Alaska","Northwest","Midwest","Northeast","Southeast","Souther Rockies","Northern Rockies"),
       cex = 0.8)
dev.off()

#PC2 separates Alaska+Rockies from the rest

#Write eigenvalues to a file

tab2 <- data.frame(sample.id = pca$sample.id,
                   EV1 = pca$eigenvect[,1],    
                   EV2 = pca$eigenvect[,2],
                   EV3 = pca$eigenvect[,3],
                   EV4 = pca$eigenvect[,4],
                   EV5 = pca$eigenvect[,5],
                   EV6 = pca$eigenvect[,6],
                   EV7 = pca$eigenvect[,7],
                   EV8 = pca$eigenvect[,8],
                   EV9 = pca$eigenvect[,9],
                   EV10 = pca$eigenvect[,10],
                   EV11 = pca$eigenvect[,11],
                   EV12 = pca$eigenvect[,12],
                   EV13 = pca$eigenvect[,13],
                   EV14 = pca$eigenvect[,14],
                   EV15 = pca$eigenvect[,15],
                   EV16 = pca$eigenvect[,16],
                   EV17 = pca$eigenvect[,17],
                   EV18 = pca$eigenvect[,18],
                   EV19 = pca$eigenvect[,19],
                   EV20 = pca$eigenvect[,20],
                   stringsAsFactors = FALSE)

write.csv(tab2,"PCA.eign.csv")

#==============================================
# 3D PCA PLOT
#==============================================

## If you want to use scores from somewhere else
# tab2 <- read.table("Picoides_pubescens.cov")
#

#install.packages("scatterplot3d")
#install.packages("rgl")
library(scatterplot3d)
library(rgl)

pcs123<-cbind(tab2$EV1,tab2$EV2,tab2$EV3)
rownames(pcs123) <- pop
color_dots <- color

pdf('PP-PCA.pdf',width=8,height=5, bg='transparent')

p<-scatterplot3d(pcs123[,1],pcs123[,2]
                 ,pcs123[,3],color=color_dots,
                 pch=20,cex.symbols=2,xlab="PC1",ylab="PC2",zlab="PC3",type="h",angle = 80,
                 xlim=c(min(pcs123[,1]),max(pcs123[,1])),
                 ylim=c(min(pcs123[,2]),max(pcs123[,2])),
                 zlim=c(min(pcs123[,3]),max(pcs123[,3])))
legend(p$xyz.convert(-0.4,0.4,-0.2),
       col= c("blue","black","green","purple","pink","red","orange"),
       bg="white", pch=c(20,20),yjust=0,
       legend = c("Alaska","Northwest","Midwest","Northeast","Southeast","Southern Rockies","Northern Rockies"),
       cex = 0.8)

#text(p$xyz.convert(pcs123[,1:3]), labels = samples,
#     cex= 0.7, col = "steelblue")

dev.off()

# pdf('PP-PCA_smaller.pdf',width=6,height=5, bg='transparent')
# 
# p<-scatterplot3d(pcs123[,1],pcs123[,2]
#                  ,pcs123[,3],color=color_dots,
#                  pch=20,cex.symbols=2,xlab="PC1",ylab="PC2",zlab="PC3",type="h",angle = 128.5,
#                  xlim=c(min(pcs123[,1]),max(pcs123[,1])),
#                  ylim=c(min(pcs123[,2]),max(pcs123[,2])),
#                  zlim=c(min(pcs123[,3]),max(pcs123[,3])))
# legend(p$xyz.convert(-0.4,0.4,-0.2),
#        col= c("blue","black","green","purple","pink","red","orange"),
#        bg="white", pch=c(20,20),yjust=0,
#        legend = c("Alaska","Northwest","Midwest","Northeast","Southeast","Southern Rockies","Northern Rockies"),
#        cex = 0.8)
# 
# #text(p$xyz.convert(pcs123[,1:3]), labels = samples,
# #     cex= 0.7, col = "steelblue")
# 
# dev.off()

#Get multidimensional distance
DIST = dist(tab2)
write.csv(as.matrix(DIST),"PCA_dist.csv")

#=============================================================================================================#
# Estimating IBD Using PLINK method of moments (MoM)
#   Calculate three IBD coefficients for non-inbred individual pairs by PLINK method of moment (MoM)
#=============================================================================================================#


# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=NULL, snp.id=snpset.id, autosome.only=FALSE, num.thread=4)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

related <- ibd.coeff[ibd.coeff[,5]!=0,]

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="Kinship of samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

#=====================================================================================================================#
# Estimating IBD Using Maximum Likelihood Estimation (MLE)
#   Calculate the three IBD coefficients (k0, k1, k2) for non-inbred individual pairs by Maximum Likelihood Estimation
#=====================================================================================================================#

ibdML <- snpgdsIBDMLE(genofile, sample.id=NULL, snp.id=snpset.id, autosome.only=FALSE, num.thread=4)

# Make a data.frame
ibdML.coeff <- snpgdsIBDSelection(ibdML)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="Kinship of samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)

#=====================================================================================================================#
# Relationship inference Using KING method of moments
#   Calculate IBD coefficients by KING method of moment.
#=====================================================================================================================#

ibd.robust <- snpgdsIBDKING(genofile, snp.id=snpset.id, autosome.only=FALSE, num.thread=4, type="KING-robust")

dat <- snpgdsIBDSelection(ibd.robust)
head(dat)

plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
     ylab="Estimated Kinship Coefficient (KING-robust)")

#=======================================================================================================================================#
# Identity-By-State Analysis
#   For the nn individuals in a sample, snpgdsIBS() can be used to create a n×nn×n matrix of genome-wide average IBS pairwise identities
#=======================================================================================================================================#

ibs <- snpgdsIBS(genofile, snp.id=snpset.id, autosome.only=FALSE, num.thread=2)

### Generate heat map

# individulas in the same population are clustered together
pop_code <- c("AK","AK","AK","AK","AK","AK","AK","AK","AK","AK","MW","MW","MW","MW","MW","MW","MW","MW","MW","MW","NE","NE","NE","NE","NE","NE","NE","NE","NE","NE","NR","NR","NR","NR","NR","NR","NR","NR","NR","NR","NW","NW","NW","NW","NW","NW","NW","NW","NW","NW","SE","SE","SE","SE","SE","SE","SE","SE","SE","SE","SR","SR","SR","SR","SR","SR","SR","SR","SR","SR")
pop.idx <- order(pop_code)

pdf('heat_map_IBS.pdf',width=5,height=5, bg='transparent')
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(30))
dev.off()

# Multiscaling

loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)

plot(x, y, col=race, xlab = "", ylab = "",
     main = "Multidimensional Scaling Analysis (IBS)")
legend("bottomright", legend=levels(race), text.col=1:nlevels(race),cex = 0.5)

#=======================================================================================================================================#
# Dissimilarity matrix
#=======================================================================================================================================#

ibs.diss <- matrix(1,70,70)-(as.matrix(ibs$ibs))
write.csv(ibs.diss,"Picoides_pubescens_IBS_dissimilarity.csv")

AK <- c(ibs.diss[1:10,1],ibs.diss[1:10,2],ibs.diss[1:10,3],ibs.diss[1:10,4],ibs.diss[1:10,5],ibs.diss[1:10,6],ibs.diss[1:10,7],ibs.diss[1:10,8],ibs.diss[1:10,9],ibs.diss[1:10,10])
MW <- c(ibs.diss[11:20,11],ibs.diss[11:20,12],ibs.diss[11:20,13],ibs.diss[11:20,14],ibs.diss[11:20,15],ibs.diss[11:20,16],ibs.diss[11:20,17],ibs.diss[11:20,18],ibs.diss[11:20,19],ibs.diss[11:20,20])
NE <- c(ibs.diss[21:30,21],ibs.diss[21:30,22],ibs.diss[21:30,23],ibs.diss[21:30,24],ibs.diss[21:30,25],ibs.diss[21:30,26],ibs.diss[21:30,27],ibs.diss[21:30,28],ibs.diss[21:30,29],ibs.diss[21:30,30])
NR <- c(ibs.diss[31:40,31],ibs.diss[31:40,32],ibs.diss[31:40,33],ibs.diss[31:40,34],ibs.diss[31:40,35],ibs.diss[31:40,36],ibs.diss[31:40,37],ibs.diss[31:40,38],ibs.diss[31:40,39],ibs.diss[31:40,40])
NW <- c(ibs.diss[41:50,41],ibs.diss[41:50,42],ibs.diss[41:50,43],ibs.diss[41:50,44],ibs.diss[41:50,45],ibs.diss[41:50,46],ibs.diss[41:50,47],ibs.diss[41:50,48],ibs.diss[41:50,49],ibs.diss[41:50,50])
SE <- c(ibs.diss[51:60,51],ibs.diss[51:60,52],ibs.diss[51:60,53],ibs.diss[51:60,54],ibs.diss[51:60,55],ibs.diss[51:60,56],ibs.diss[51:60,57],ibs.diss[51:60,58],ibs.diss[51:60,59],ibs.diss[51:60,60])
SR <- c(ibs.diss[61:70,61],ibs.diss[61:70,62],ibs.diss[61:70,63],ibs.diss[61:70,64],ibs.diss[61:70,65],ibs.diss[61:70,66],ibs.diss[61:70,67],ibs.diss[61:70,68],ibs.diss[61:70,69],ibs.diss[61:70,70])

data_for_box_plot <- data.frame(AK=AK[which(AK!=0)],MW=MW[which(NW!=0)],NE=NE[which(NE!=0)],NR=NR[which(NR!=0)],NW=NW[which(NW!=0)],SE=SE[which(SE!=0)],SR=SR[which(SR!=0)])
boxplot(data_for_box_plot,main="Intrapopulation dissimilarity")

#or

diss <- snpgdsDiss(genofile, sample.id=NULL, snp.id=snpset.id, autosome.only=FALSE,
           remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
ibs.diss_2 <- diss$diss

mantel.randtest(as.dist(ibs.diss),as.dist(ibs.diss_2))

#=======================================================================================================================================#
# FST Estimation
#   Given two or more populations, FstFst can be estimated by the method of Weir & Cockerham (1984).
#=======================================================================================================================================#

#pop_code <- scan("pop.txt", what=character())
pop_code <- c("AK","AK","AK","AK","AK","AK","AK","AK","AK","AK","MW","MW","MW","MW","MW","MW","MW","MW","MW","MW","NE","NE","NE","NE","NE","NE","NE","NE","NE","NE","NR","NR","NR","NR","NR","NR","NR","NR","NR","NR","NW","NW","NW","NW","NW","NW","NW","NW","NW","NW","SE","SE","SE","SE","SE","SE","SE","SE","SE","SE","SR","SR","SR","SR","SR","SR","SR","SR","SR","SR")
fst <- snpgdsFst(genofile, population=as.factor(pop_code), autosome.only=FALSE, method="W&C84")

FST <- fst$Fst # Weighted Fst (weighted by sample size per SNP ?)
Mean_FST <- fst$MeanFst # Mean Fst across SNPs
FST_SNPS <- fst$FstSNP # Fst per SNP

#=======================================================================================================================================#
# Cluster Analysis
#=======================================================================================================================================#

#To perform cluster analysis on the n×nn×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score
set.seed(100)
ibs.hc <- snpgdsHCluster(ibs)

# Determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc)

pdf('dendrogram.pdf',width=5,height=5, bg='transparent')
plot(rv$dendrogram, leaflab="none", main="Cluster Analysis")
dev.off()

# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

pdf('Cluster_dendrogram.pdf',width=5,height=5, bg='transparent')
plot(rv2$dendrogram, leaflab="none", main="Cluster Analysis")
legend("bottom", legend=levels(as.factor(pop_code)), col=1:nlevels(as.factor(pop_code)), pch=19, ncol=5,xjust=0.5,yjust=0.5,cex=0.5)
dev.off()

#### With dissimilarity matrix

#To perform cluster analysis on the n×nn×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score
set.seed(100)
diss.hc <- snpgdsHCluster(diss)

# Determine groups of individuals automatically
rvdiss <- snpgdsCutTree(diss.hc)

pdf('dendrogram.pdf',width=5,height=5, bg='transparent')
plot(rvdiss$dendrogram, leaflab="none", main="Cluster Analysis")
dev.off()

# Determine groups of individuals by population information
rv2diss <- snpgdsCutTree(diss.hc, samp.group=as.factor(pop_code))

pdf('Cluster_dendrogram.pdf',width=5,height=5, bg='transparent')
plot(rv2diss$dendrogram, leaflab="none", main="Cluster Analysis")
legend("bottom", legend=levels(as.factor(pop_code)), col=1:nlevels(as.factor(pop_code)), pch=19, ncol=5,xjust=0.5,yjust=0.5,cex=0.5)
dev.off()
