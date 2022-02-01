#Picoides pubescens

popPP<-read.table("PPpop.info.txt",as.is=T)

qpotPP <- list.files(pattern ="Picoides_pubescens_NGSadmix.qopt")

pdf("PP-K2.pdf",8,4)
admixPP<-t(as.matrix(read.table(qpotPP[1])))
admixPP<-admixPP[,order(popPP[,1])]
popPP<-popPP[order(popPP[,1]),]
h<-barplot(admixPP,col=c("blue","orange"),space=c(0.3,rep(0.1,9)),border=NA,xlab="Individuals",ylab="Admixture")
text(c(6,17.2,28.4,39.6,50.8,62,73.2),-0.05,unique(popPP[,1]),xpd=T)
dev.off()

pdf("PP-K3.pdf",8,4)
admixPP<-t(as.matrix(read.table(qpotPP[2])))
admixPP<-admixPP[,order(popPP[,1])]
popPP<-popPP[order(popPP[,1]),]
h<-barplot(admixPP,col=c("blue","orange","green"),space=c(0.3,rep(0.1,9)),border=NA,xlab="Individuals",ylab="Admixture")
text(c(6,17.2,28.4,39.6,50.8,62,73.2),-0.05,unique(popPP[,1]),xpd=T)
dev.off()

pdf("PP-K4.pdf",8,4)
admixPP<-t(as.matrix(read.table(qpotPP[3])))
admixPP<-admixPP[,order(popPP[,1])]
popPP<-popPP[order(popPP[,1]),]
h<-barplot(admixPP,col=c("red","orange","blue","green"),space=c(0.3,rep(0.1,9)),border=NA,xlab="Individuals",ylab="Admixture")
text(c(6,17.2,28.4,39.6,50.8,62,73.2),-0.05,unique(popPP[,1]),xpd=T)
dev.off()

#Picoides villosus

popPV<-read.table("PVpop.info.txt",as.is=T)
qpotPV <- list.files(pattern ="Picoides_villosus_NGSadmix.qopt")

pdf("PV-K2.pdf",8,4)
admixPV<-t(as.matrix(read.table(qpotPV[1])))
admixPV<-admixPV[,order(popPV[,1])]
popPV<-popPV[order(popPV[,1]),]
h<-barplot(admixPV,col=c("blue","orange"),space=c(0.3,rep(0.1,9)),border=NA,xlab="Individuals",ylab="Admixture")
text(c(6,17.2,28.4,39.6,50.8,62,73.2),-0.05,unique(popPV[,1]),xpd=T)
dev.off()

pdf("PV-K3.pdf",8,4)
admixPV<-t(as.matrix(read.table(qpotPV[2])))
admixPV<-admixPV[,order(popPV[,1])]
popPV<-popPV[order(popPV[,1]),]
h<-barplot(admixPV,col=c("green","blue","orange"),space=c(0.3,rep(0.1,9)),border=NA,xlab="Individuals",ylab="Admixture")
text(c(6,17.2,28.4,39.6,50.8,62,73.2),-0.05,unique(popPV[,1]),xpd=T)
dev.off()

pdf("PV-K4.pdf",8,4)
admixPV<-t(as.matrix(read.table(qpotPV[3])))
admixPV<-admixPV[,order(popPV[,1])]
popPV<-popPV[order(popPV[,1]),]
h<-barplot(admixPV,col=c("green","orange","blue","red"),space=c(0.3,rep(0.1,9)),border=NA,xlab="Individuals",ylab="Admixture")
text(c(6,17.2,28.4,39.6,50.8,62,73.2),-0.05,unique(popPV[,1]),xpd=T)
dev.off()
