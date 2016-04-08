library(BEST)
source("~/scripts/R/dna.R")

mobio<-read.csv('mobio.csv')
psp<-read.csv('psp.csv')

#mobioNegs<-na.omit(unlist(mobio[,c('OR.Air.Swab','Sterile.Swab','Water','Extraction.Blank')]))
mobioNegs<-na.omit(unlist(mobio[,c('OR.Air.Swab')]))
mobioPlacenta<-na.omit(unlist(mobio[,c('Placenta..MS.','Placenta..FS.')]))
mobioSaliva<-na.omit(unlist(mobio[,c('Saliva')]))
mobioOtherNegs<-na.omit(unlist(mobio[,c('Sterile.Swab','Water','Extraction.Blank')]))
mobioSwab<-na.omit(unlist(mobio[,c('Sterile.Swab')]))
mobioBlank<-na.omit(unlist(mobio[,c('Extraction.Blank')]))
mobioWater<-na.omit(unlist(mobio[,c('Water')])) #only 1 samples

#pspNegs<-na.omit(unlist(psp[,c('OR.Air.Swab','Sterile.Swab','Water','Extraction.Blank')]))
pspNegs<-na.omit(unlist(psp[,c('OR.Air.Swab')]))
pspPlacenta<-na.omit(unlist(psp[,c('Placenta..MS.','Placenta..FS.')]))
pspSaliva<-na.omit(unlist(psp[,c('Saliva')]))
pspOtherNegs<-na.omit(unlist(psp[,c('Sterile.Swab','Water','Extraction.Blank')]))

pspQuick<-BESTmcmc(pspNegs,pspPlacenta)
mobioQuick<-BESTmcmc(mobioNegs,mobioPlacenta)
mobioBayes<-cacheOperation('mobioBayes.Rdat',BESTmcmc,mobioNegs,mobioPlacenta,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
pspBayes<-cacheOperation('pspBayes.Rdat',BESTmcmc,pspNegs,pspPlacenta,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
pspDiff<-pspBayes$mu1-pspBayes$mu2
mobioDiff<-mobioBayes$mu1-mobioBayes$mu2
pspDiffQuick<-pspQuick$mu1-pspQuick$mu2
mobioDiffQuick<-mobioQuick$mu1-mobioQuick$mu2

mobioSaliva<-cacheOperation('salivaMobioBayes.Rdat',BESTmcmc,mobioNegs,mobioSaliva,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
pspSaliva<-cacheOperation('salivapspBayes.Rdat',BESTmcmc,pspNegs,pspSaliva,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
pspSalivaDiff<-pspSaliva$mu1-pspSaliva$mu2
mobioSalivaDiff<-mobioSaliva$mu1-mobioSaliva$mu2

mobioOther<-cacheOperation('otherMobioBayes.Rdat',BESTmcmc,mobioNegs,mobioOtherNegs,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
mobioSwabBest<-cacheOperation('swabMobioBayes.Rdat',BESTmcmc,mobioNegs,mobioSwab,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
mobioBlankBest<-cacheOperation('blankMobioBayes.Rdat',BESTmcmc,mobioNegs,mobioBlank,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
#mobioWaterBest<-cacheOperation('waterMobioBayes.Rdat',BESTmcmc,mobioNegs,mobioWater,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
pspOther<-cacheOperation('otherPspBayes.Rdat',BESTmcmc,pspNegs,pspOtherNegs,burnInSteps=2e5,thinSteps=10,numSavedSteps=2e6,parallel=TRUE)
pspOtherDiff<-pspOther$mu1-pspOther$mu2
mobioOtherDiff<-mobioOther$mu1-mobioOther$mu2
mobioSwabDiff<-mobioSwabBest$mu1-mobioSwabBest$mu2
mobioBlankDiff<-mobioBlankBest$mu1-mobioBlankBest$mu2
#mobioWaterDiff<-mobioWaterBest$mu1-mobioWaterBest$mu2




xlim<-c(-5,25)
#prettyLabels<-pretty(c(-5,25))
prettyLabels<-pretty(log10(2^c(-5,25)))
bins<-seq(xlim[1],xlim[2],.1)
pspTab<-table(cut(pspDiff,bins))/length(pspDiff)
mobioTab<-table(cut(mobioDiff,bins))/length(mobioDiff)
pspSalivaTab<-table(cut(pspSalivaDiff,bins))/length(pspSalivaDiff)
mobioSalivaTab<-table(cut(mobioSalivaDiff,bins))/length(mobioSalivaDiff)
pspOtherTab<-table(cut(pspOtherDiff,bins))/length(pspOtherDiff)
mobioOtherTab<-table(cut(mobioOtherDiff,bins))/length(mobioOtherDiff)
mobioSwabTab<-table(cut(mobioSwabDiff,bins))/length(mobioSwabDiff)
mobioBlankTab<-table(cut(mobioBlankDiff,bins))/length(mobioBlankDiff)

pdf('mobioPsp.pdf')
par(mfrow=c(2,1),las=1,mar=c(4,3.6,1.1,.1))
plot(1,1,type='n',xlim=xlim,ylim=range(pspTab,pspSalivaTab,pspOtherTab),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),main='PSP')
#axis(1,prettyLabels,sapply(prettyLabels,function(x)as.expression(bquote(2^.(x)))),las=1)
axis(1,log2(10^prettyLabels),ifelse(prettyLabels==0,1,sapply(prettyLabels,function(x)as.expression(bquote(10^.(x))))),las=1)
title(xlab='Fold increase in DNA over air swabs',mgp=c(2,1,0))
polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(pspTab),0,0),col='#0000FF44')
polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(pspSalivaTab),0,0),col='#FF000044')
polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(pspOtherTab),0,0),col='#00FF0044')
abline(v=0,lty=2)
legend('topright',c('Placenta','Saliva','Negatives'),fill=c('#0000FF44','#FF000044','#00FF0044'),bty='n')
plot(1,1,type='n',xlim=xlim,ylim=range(mobioTab,mobioSalivaTab,mobioOtherTab),xlab='',xaxt='n',ylab='Posterior probability',mgp=c(2.7,.8,0),main='MoBio')
#axis(1,prettyLabels,sapply(prettyLabels,function(x)as.expression(bquote(2^.(x)))),las=1)
axis(1,log2(10^prettyLabels),ifelse(prettyLabels==0,1,sapply(prettyLabels,function(x)as.expression(bquote(10^.(x))))),las=1)
title(xlab='Fold increase in DNA over air swabs',mgp=c(2,1,0))
polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(mobioTab),0,0),col='#0000FF44')
polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(mobioSalivaTab),0,0),col='#FF000044')
polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(mobioOtherTab),0,0),col='#00FF0044')
#polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(mobioSwabTab),0,0),col='#00FF0044')
#polygon(c(-5,bins[-length(bins)]+.05,25,-5),c(0,unlist(mobioBlankTab),0,0),col='#00FF0044')
abline(v=0,lty=2)
legend('topright',c('Placenta','Saliva','Negatives'),fill=c('#0000FF44','#FF000044','#00FF0044'),bty='n')
dev.off()

