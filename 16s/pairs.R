###
pairedSampleTypes<-apply(table(samples$id,samples$SampleType)>1,2,any)
pairedSampleTypes<-names(pairedSampleTypes)[pairedSampleTypes]

pairs<-lapply(pairedSampleTypes,function(type){
  thisPsp<-samples[samples$SampleType==type&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit',]
  thisPsp<-thisPsp[order(thisPsp$id),]
  thisMobio<-samples[samples$SampleType==type&samples$ExtractionKit=='Mobio Powersoil',]
  thisMobio<-thisMobio[order(thisMobio$id),]
  if(any(thisPsp$id!=thisMobio$id))stop(simpleError('IDs do not match between PSP and Mobio'))
  paired<-mapply(function(pspId,mobioId)otus[,c(pspId,mobioId)],thisPsp$SampleID,thisMobio$SampleID,SIMPLIFY=FALSE)
  paired<-lapply(paired,function(x)apply(x,2,function(y)y/sum(y)))
  return(paired)
})
names(pairs)<-pairedSampleTypes

pspProp<-apply(pspControl,2,function(x)x/sum(x))
mobioProp<-apply(mobioControl,2,function(x)x/sum(x))
pairMin<-lapply(pairs,function(x)do.call(cbind,lapply(x,function(y)apply(y,1,min))))
pairMinMax<-do.call(cbind,lapply(pairMin,apply,1,max))

propCuts<-10^seq(log10(min(pairMinMax[pairMinMax>0])),0,.005)
propVsN<-do.call(rbind,dnar::cacheOperation('work/propVsN.Rdat',mclapply,propCuts,function(propCut){
  cat('.')
  pairMatch<-pairMinMax>propCut
  anyPsp<-apply(pspProp>propCut,1,sum)>0
  anyMobio<-apply(mobioProp>propCut,1,sum)>0
  noControl<-apply(pairMatch,2,function(x)x&!anyPsp&!anyMobio)
  return(apply(noControl,2,sum))
},mc.cores=4))

cols<-rev(rainbow(ncol(propVsN),alpha=.7))
names(cols)<-colnames(propVsN)
pdf('out/propVsN.pdf',width=4,height=4)
  par(mar=c(2.6,4.2,.1,.2),las=1,lheight=.7,mgp=c(2.5,.8,0))
  plot(1,1,type='n',xlab='',ylab='Number of OTUs shared by both\nextractions but not control',ylim=range(propVsN+1),xlim=range(propCuts),log='xy',xaxt='n',yaxt='n')
  title(xlab='Presence cutoff (proportion of sample)',mgp=c(1.5,1,0))
  prettyX<-pretty(log10(propCuts))
  axis(1,10^prettyX,sapply(prettyX,function(x)ifelse(x<0,as.expression(bquote(10^.(x))),10^x)),las=1,mgp=c(3,.5,0))
  prettyY<-c(0,1,5,10,50,100,500,1000)
  axis(2,prettyY+1,as.character(prettyY))
  ticks<-c(1:10,seq(20,100,10),seq(200,1000,100),seq(2000,10000,1000))+1
  axis(2,ticks,rep('',length(ticks)),tcl=-.25)
  sapply(colnames(propVsN),function(x)lines(propCuts,propVsN[,x]+1,lwd=2,col=cols[x]))
  legend('topright',names(cols),col=cols,bty='n',lty=1,lwd=2)
dev.off()

