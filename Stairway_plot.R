Alaska <- read.table("Alaska.final.summary",header=T)
Eastern <- read.table("Eastern.final.summary",header=T)
Northwest <- read.table("Northwest.final.summary",header=T)
Rockies <- read.table("Rockies.final.summary",header=T)

pdf("PP-stairway2.plot.pdf",height = 6.25,width = 10)

plot(Alaska$year/1000,Alaska$Ne_median/1000, log=c("xy"), type="n", xlab="Time (1k years ago)", ylab="Effective Population Size (1k individuals)",xlim=c(4,1000),ylim=c(20,30000))

#LGM
abline(v=22,col="grey")

#LGP
abline(v=115,col="grey")

lines(Alaska$year/1000,Alaska$Ne_median/1000,type="s",col="orange",lwd = 5)
lines(Alaska$year/1000,Alaska$Ne_2.5./1000,type="s",col="orange",lty=3)
lines(Alaska$year/1000,Alaska$Ne_97.5./1000,type="s",col="orange",lty=3)

lines(Eastern$year/1000,Eastern$Ne_median/1000,type="s",col="green",lwd = 5)
lines(Eastern$year/1000,Eastern$Ne_2.5./1000,type="s",col="green",lty=3)
lines(Eastern$year/1000,Eastern$Ne_97.5./1000,type="s",col="green",lty=3)

lines(Northwest$year/1000,Northwest$Ne_median/1000,type="s",col="red",lwd = 5)
lines(Northwest$year/1000,Northwest$Ne_2.5./1000,type="s",col="red",lty=3)
lines(Northwest$year/1000,Northwest$Ne_97.5./1000,type="s",col="red",lty=3)

lines(Rockies$year/1000,Rockies$Ne_median/1000,type="s",col="blue",lwd = 5)
lines(Rockies$year/1000,Rockies$Ne_2.5./1000,type="s",col="blue",lty=3)
lines(Rockies$year/1000,Rockies$Ne_97.5./1000,type="s",col="blue",lty=3)

legend(400,28000,legend = c("Alaska","Eastern","Northwest","Rockies","95% CI"),col=c("orange","green","red","blue","black"),lty=c(1,1,1,1,3),lwd=c(5,5,5,5,1),cex=0.8)

dev.off()

######################## OLNY RECENT ########################

pdf("PP-stairway2_onlyrecent.plot.pdf",height = 6.25,width = 10)

plot(Alaska$year/1000,Alaska$Ne_median/1000, log=c("xy"), type="n", xlab="Time (1k years ago)", ylab="Effective Population Size (1k individuals)",xlim=c(4,200),ylim=c(20,30000))

#LGM
abline(v=22,col="grey")

#LGP
abline(v=115,col="grey")

lines(Alaska$year/1000,Alaska$Ne_median/1000,type="s",col="orange",lwd = 5)
lines(Alaska$year/1000,Alaska$Ne_2.5./1000,type="s",col="orange",lty=3)
lines(Alaska$year/1000,Alaska$Ne_97.5./1000,type="s",col="orange",lty=3)

lines(Eastern$year/1000,Eastern$Ne_median/1000,type="s",col="green",lwd = 5)
lines(Eastern$year/1000,Eastern$Ne_2.5./1000,type="s",col="green",lty=3)
lines(Eastern$year/1000,Eastern$Ne_97.5./1000,type="s",col="green",lty=3)

lines(Northwest$year/1000,Northwest$Ne_median/1000,type="s",col="red",lwd = 5)
lines(Northwest$year/1000,Northwest$Ne_2.5./1000,type="s",col="red",lty=3)
lines(Northwest$year/1000,Northwest$Ne_97.5./1000,type="s",col="red",lty=3)

lines(Rockies$year/1000,Rockies$Ne_median/1000,type="s",col="blue",lwd = 5)
lines(Rockies$year/1000,Rockies$Ne_2.5./1000,type="s",col="blue",lty=3)
lines(Rockies$year/1000,Rockies$Ne_97.5./1000,type="s",col="blue",lty=3)

legend(4,200,legend = c("Alaska","Eastern","Northwest","Rockies","95% CI"),col=c("orange","green","red","blue","black"),lty=c(1,1,1,1,3),lwd=c(5,5,5,5,1),cex=0.8)

dev.off()

################# SAME SCALE AS HAIRY #################

pdf("PP-stairway2_same-scale.plot.pdf",height = 6.25,width = 10)

plot(Alaska$year/1000,Alaska$Ne_median/1000, log=c("xy"), type="n", xlab="Time (1k years ago)", ylab="Effective Population Size (1k individuals)",xlim=c(2,2000),ylim=c(2,50000))

#LGM
abline(v=22,col="grey")

#LGP
abline(v=115,col="grey")

lines(Alaska$year/1000,Alaska$Ne_median/1000,type="s",col="orange",lwd = 5)
lines(Alaska$year/1000,Alaska$Ne_2.5./1000,type="s",col="orange",lty=3)
lines(Alaska$year/1000,Alaska$Ne_97.5./1000,type="s",col="orange",lty=3)

lines(Eastern$year/1000,Eastern$Ne_median/1000,type="s",col="green",lwd = 5)
lines(Eastern$year/1000,Eastern$Ne_2.5./1000,type="s",col="green",lty=3)
lines(Eastern$year/1000,Eastern$Ne_97.5./1000,type="s",col="green",lty=3)

lines(Northwest$year/1000,Northwest$Ne_median/1000,type="s",col="red",lwd = 5)
lines(Northwest$year/1000,Northwest$Ne_2.5./1000,type="s",col="red",lty=3)
lines(Northwest$year/1000,Northwest$Ne_97.5./1000,type="s",col="red",lty=3)

lines(Rockies$year/1000,Rockies$Ne_median/1000,type="s",col="blue",lwd = 5)
lines(Rockies$year/1000,Rockies$Ne_2.5./1000,type="s",col="blue",lty=3)
lines(Rockies$year/1000,Rockies$Ne_97.5./1000,type="s",col="blue",lty=3)

legend(800,48000,legend = c("Alaska","Eastern","Northwest","Rockies","95% CI"),col=c("orange","green","red","blue","black"),lty=c(1,1,1,1,3),lwd=c(5,5,5,5,1),cex=0.8)

dev.off()

