## Made by Lucas R. Moreira
## Last updated 14 April 2020
## Usage: Predictors of genetic diversity along the genome of P. pubescens

data <- read.csv("Predictors_of_theta.csv",header=T)

# Univariate models #

hist(data$thetaP/100000)

max(data$thetaP/100000)
min(data$thetaP/100000)

mean(data$recR)
max(data$recR)
min(data$recR)

# theta pi is normally distributed, being therefore suited for linear model

# Effect of recombination rates

sink("lm_thetaP~recR.txt")
model1 <- lm(scale(data$thetaP/100000)~scale(data$recR))
summary(model1)
sink()

pdf("Theta~RecR.pdf",width = 8,height = 6)
plot(data$recR,data$thetaP/100000,ylab="Theta pi",xlab="Recombination Rate")
abline(lm(data$thetaP/100000~data$recR),col="red")
dev.off()


# Effect of Gene Density

sink("lm_thetaP~gene density.txt")
model2 <- lm(scale(data$thetaP/100000)~scale(data$coding_seq_percent))
summary(model2)
sink()

pdf("Theta~Gene density.pdf",width = 8,height = 6)
plot(data$coding_seq_percent,data$thetaP/100000,ylab="Theta pi",xlab="Gene Density (% of coding sequence)")
abline(lm(data$thetaP/100000~data$coding_seq_percent),col="red")
dev.off()

model3 <- lm(scale(data$thetaP/100000)~scale(data$count_cds))
summary(model3)

pdf("Theta~Gene density_count.pdf",width = 8,height = 6)
plot(data$count_cds,data$thetaP/100000,ylab="Theta pi",xlab="Gene Density (number of genes)")
abline(lm(data$thetaP/100000~data$count_cds),col="red")
dev.off()

# Effect of GC content

sink("lm_thetaP~GC_content.txt")
model4 <- lm(scale(data$thetaP/100000)~scale(data$GC_content))
summary(model4)
sink()

pdf("Theta~GC_content.pdf",width = 8,height = 6)
plot(data$GC_content,data$thetaP/100000,ylab="Theta pi",xlab="GC content")
abline(lm(data$thetaP/100000~data$GC_content),col="red")
dev.off()

# Multivariate model #

sink("lm_thetaP~all_variables.txt")
model5 <- lm(scale(data$thetaP/100000)~scale(data$recR)+scale(data$coding_seq_percent)+scale(data$GC_content)+scale(data$recR)*scale(data$coding_seq_percent)+scale(data$recR)*scale(data$GC_content)+scale(data$coding_seq_percent)*scale(data$GC_content)+scale(data$recR)*scale(data$coding_seq_percent)*scale(data$recR))
summary(model5)
sink()

# Let's see how correlated the predictor variables are

cor(data[,4:10])
heatmap(cor(data[complete.cases(data),4:10]))

sink("corr.txt")
cor.test(data$recR,data$GC_content, method = "pearson") #not significant
cor.test(data$recR,data$coding_seq_percent, method = "pearson")
cor.test(data$coding_seq_percent,data$GC_content, method = "pearson")
sink()

data2 <- na.exclude(data)

library(corrplot)
M <- cor(data2[,c(5,7,8,11)])
colnames(M) <- c("Theta pi","Recombination rate","Gene density","GC content")
rownames(M) <- c("Theta pi","Recombination rate","Gene density","GC content")
corrplot(M, method = "number", type = "upper")

pdf("Correlations.pdf")
p.mat <- cor.mtest(data2[,c(5,7,8,11)])$p
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)
dev.off()

# Principal component regression

library(pls)
library("factoextra")

sink("PCR.txt")
pcr_model <- pcr(thetaP/100000~recR+coding_seq_percent+GC_content, data = data, scale = TRUE, validation = "CV")
summary(pcr_model)
res.pca <- prcomp(data[,c(7,8,11)], scale = TRUE, center = TRUE)
res.var <- get_pca_var(res.pca)
res.var$contrib
summary(lm(data$thetaP/100000~pcr_model$scores[,1]))
summary(lm(data$thetaP/100000~pcr_model$scores[,2]))
summary(lm(data$thetaP/100000~pcr_model$scores[,3]))
m.model <- lm(data$thetaP/100000~pcr_model$scores[,1]+pcr_model$scores[,2]+pcr_model$scores[,3])

af <- anova(m.model)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
sink()

###################################################
######### PLOTTING DENSITIES WITH GGPLOT  #########
###################################################

library(MASS)
library(ggplot2)
library(viridis)
#> Loading required package: viridisLite
theme_set(theme_bw(base_size = 16))

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

data2$density <- get_density(data2$recR, data2$thetaP/100000, n = 100)
pdf("Theta~RecR_densities.pdf",width = 8,height = 6)
ggplot(data2, aes(data2$recR, data2$thetaP/100000)) + geom_point(aes(recR, thetaP/100000, color = density)) + scale_color_viridis() + ylab("Theta pi") + xlab("Recombination Rate") + labs(color="Density") + geom_smooth(method = "lm",colour="red") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

data2$density <- get_density(data2$coding_seq_percent, data2$thetaP/100000, n = 100)
pdf("Theta~Gene_density_densities.pdf",width = 8,height = 6)
ggplot(data2, aes(data2$coding_seq_percent, data2$thetaP/100000)) + geom_point(aes(coding_seq_percent, thetaP/100000, color = density)) + scale_color_viridis() + ylab("Theta pi") + xlab("Gene density") + labs(color="Density") + geom_smooth(method = "lm",colour="red") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

data2$density <- get_density(data2$GC_content, data2$thetaP/100000, n = 100)
pdf("Theta~GC_content_densities.pdf",width = 8,height = 6)
ggplot(data2, aes(data2$GC_content, data2$thetaP/100000)) + geom_point(aes(GC_content, thetaP/100000, color = density)) + scale_color_viridis() + ylab("Theta pi") + xlab("GC content") + labs(color="Density") + geom_smooth(method = "lm",colour="red") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


# Are the top 25% of rec rates still correlated with diversity?

top_data <- data[with(data,order(-recR)),]
top_data1 <- top_data[1:115,]

model_top <- lm(scale(top_data1$thetaP/100000)~scale(top_data1$recR))
summary(model_top)

plot(top_data1$recR,top_data1$thetaP/100000)


# LOESS model
library(caret)

dt <- data
train.control <- trainControl(method = "cv", p = 0.80, number = 20)
lm1 <- train(thetaP~recR, data = dt, method = "gamLoess", na.action = na.omit, trControl = train.control)
model_results1 <- lm1$results
sink("recR_vs_thetaP-LOESS.txt")
model_results1
sink()

lm2 <- train(thetaP~coding_seq_percent, data = dt, method = "gamLoess", na.action = na.omit, trControl = train.control)
model_results2 <- lm2$results
sink("coding_seq_vs_thetaP-LOESS.txt")
model_results2
sink()

library("ggpubr")
library(gghighlight)
library(ggplot2)
cor_plot <- ggscatter(dt, y = "thetaP", x = "recR", color = "coding_seq_percent", size = 0.7, alpha=1.0,
                      ylab = "Nucleotide Diversity (Pi)", xlab = "Recombination rate c/bp")  + #sc1 + theme_bw()# +
  scale_color_gradient(low = "blue", high = "green", name= "Gene density\nquantiles") + theme_minimal(base_size = 18) +
  stat_smooth(geom = "smooth", position = "identity", method = "loess", se = TRUE, n = 80,
              span = 0.75, fullrange = FALSE, level = 0.95, method.args = list(), na.rm = T, inherit.aes = TRUE) 

pdf("recR_vs_thetaP-LOESS.pdf",width=8,height=5)
cor_plot
dev.off()

library("ggpubr")
library(gghighlight)
library(ggplot2)
cor_plot <- ggscatter(dt, y = "thetaP", x = "coding_seq_percent", color = "recR", size = 0.7, alpha=1.0,
                      ylab = "Nucleotide Diversity (Pi)", xlab = "Gene Density (%)")  + #sc1 + theme_bw()# +
  scale_color_gradient(low = "blue", high = "green", name= "Recombination rate\nquantiles") + theme_minimal(base_size = 18) +
  stat_smooth(geom = "smooth", position = "identity", method = "loess", se = TRUE, n = 80,
              span = 0.75, fullrange = FALSE, level = 0.95, method.args = list(), na.rm = T, inherit.aes = TRUE) 

pdf("coding_seq_vs_thetaP-LOESS.pdf",width=8,height=5)
cor_plot
dev.off()
