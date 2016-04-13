###
pairedSampleTypes<-apply(table(ifelse(samples$id=='',NA,samples$id),samples$type)>1,2,any)
pairedSampleTypes<-names(pairedSampleTypes)[pairedSampleTypes]

secondMax<-function(x)tail(sort(x),2)[1]
singles<-lapply(unique(samples$airSplit),function(type){
  thisPsp<-samples[samples$airSplit==type&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit',]
  thisMobio<-samples[samples$airSplit==type&samples$ExtractionKit=='Mobio Powersoil',]
  if(nrow(thisPsp)>0)pspProps<-apply(otuProp[,thisPsp$SampleID],1,secondMax)
  else pspProps<-NULL
  if(nrow(thisMobio)>0)mobioProps<-apply(otuProp[,thisMobio$SampleID],1,secondMax)
  else mobioProps<-NULL
  out<-list( 'PSP'=pspProps, 'MoBio'=mobioProps)
  return(out)
})
names(singles)<-unique(samples$airSplit)
singleTreat<-lapply(c('MoBio','PSP'),function(x)do.call(cbind,lapply(singles,'[[',x)))
names(singleTreat)<-c('MoBio','PSP')

pairs<-lapply(pairedSampleTypes,function(type){
  thisPsp<-samples[samples$type==type&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit',]
  thisPsp<-thisPsp[order(thisPsp$id),]
  thisMobio<-samples[samples$type==type&samples$ExtractionKit=='Mobio Powersoil',]
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
pairMax<-lapply(pairs,function(x)do.call(cbind,lapply(x,function(y)apply(y,1,max))))
pairMaxMax<-do.call(cbind,lapply(pairMax,apply,1,max))

anyPlacenta<-apply(pairMinMax[,c('Placenta (FS)','Placenta (MS)')],1,max)>0
out<-data.frame('ID'=otus[anyPlacenta,'OTU.ID'],'Taxa'=otus[anyPlacenta,'Consensus.Lineage'],pairMinMax[anyPlacenta,],'PSP controls'=apply(pspProp[anyPlacenta,],1,max),'MoBio controls'=apply(mobioProp[anyPlacenta,],1,max),check.names=FALSE)
write.csv(out,'out/otusInPlacenta.csv',row.names=FALSE)


propCuts<-10^seq(log10(min(pairMinMax[pairMinMax>0])),0,.005)
propVsN<-do.call(rbind,dnar::cacheOperation('work/propVsN.Rdat',mclapply,propCuts,function(propCut){
  cat('.')
  pairMatch<-pairMinMax>propCut
  anyPsp<-apply(pspProp>propCut,1,sum)>0
  anyMobio<-apply(mobioProp>propCut,1,sum)>0
  noControl<-apply(pairMatch,2,function(x)x&!anyPsp&!anyMobio)
  return(apply(noControl,2,sum))
},mc.cores=4))

totalN<-do.call(rbind,dnar::cacheOperation('work/totalN.Rdat',mclapply,propCuts,function(propCut){
  cat('.')
  pairMatch<-pairMinMax>propCut
  return(apply(pairMatch,2,sum))
},mc.cores=4))

totalUnpairedN<-do.call(rbind,dnar::cacheOperation('work/totalUnpairedN.Rdat',mclapply,propCuts,function(propCut){
  cat('.')
  pairMatch<-pairMaxMax>propCut
  return(apply(pairMatch,2,sum))
},mc.cores=4))


singlePropVsN<-dnar::cacheOperation('work/singlePropVsN.Rdat',mclapply,propCuts,function(propCut,...){
  cat('.')
  thisPsp<-singleTreat[['PSP']][,colnames(singleTreat[['PSP']])!='Air Swab']
  thisMobio<-singleTreat[['MoBio']][,colnames(singleTreat[['MoBio']])!='Air Swab']
  thisPspControl<-singleTreat[['PSP']][,'Air Swab']>propCut
  thisMobioControl<-singleTreat[['MoBio']][,'Air Swab']>propCut
  pspSums<-apply(thisPsp>propCut,2,function(x)sum(x&!thisPspControl))
  mobioSums<-apply(thisMobio>propCut,2,function(x)sum(x&!thisMobioControl))
  pspTotal<-apply(thisPsp>propCut,2,function(x)sum(x))
  mobioTotal<-apply(thisMobio>propCut,2,function(x)sum(x))
  return(list('notControl'=list('PSP'=pspSums,'MoBio'=mobioSums),'total'=list('PSP'=pspTotal,'MoBio'=mobioTotal)))
},mc.cores=4)
singleSums<-lapply(c('MoBio','PSP'),function(x)do.call(rbind,lapply(singlePropVsN,function(y)y[['notControl']][[x]])))
singleTotals<-lapply(c('MoBio','PSP'),function(x)do.call(rbind,lapply(singlePropVsN,function(y)y[['total']][[x]])))
names(singleSums)<-names(singleTotals)<-c('MoBio','PSP')



cols<-rainbow.lab(ncol(propVsN),alpha=.7,lightMultiple=.65,lightScale=0)
cols2<-rainbow.lab(ncol(propVsN),alpha=.4,lightMultiple=.9,lightScale=0)
names(cols)<-names(cols2)<-colnames(propVsN)
for(showTotal in c('total','totalUnpaired','noTotal')){
pdf(sprintf('out/propVsN_%s.pdf',showTotal),width=4,height=4)
  par(mar=c(2.6,4.2,.1,.2),las=1,lheight=.7,mgp=c(2.5,.8,0))
  plot(1,1,type='n',xlab='',ylab='Number of OTUs shared by both\nextractions but not control',ylim=range(propVsN+1),xlim=range(propCuts),log='xy',xaxt='n',yaxt='n')
  title(xlab='Presence cutoff (proportion of sample)',mgp=c(1.5,1,0))
  prettyX<-pretty(log10(propCuts))
  axis(1,10^prettyX,sapply(prettyX,function(x)ifelse(x<0,as.expression(bquote(10^.(x))),10^x)),las=1,mgp=c(3,.5,0))
  prettyY<-c(0,1,5,10,50,100,500,1000)
  axis(2,prettyY+1,as.character(prettyY))
  ticks<-c(1:10,seq(20,100,10),seq(200,1000,100),seq(2000,10000,1000))+1
  axis(2,ticks,rep('',length(ticks)),tcl=-.25)
  if(showTotal=='total')sapply(colnames(propVsN),function(x)polygon(c(propCuts,rev(propCuts)),c(propVsN[,x],rev(totalN[,x]))+1,col=cols2[x],border=NA))#,border=cols[x]))
  if(showTotal=='totalUnpaired')sapply(colnames(propVsN),function(x)polygon(c(propCuts,rev(propCuts)),c(propVsN[,x],rev(totalUnpairedN[,x]))+1,col=cols2[x],border=NA))#,border=cols[x]))
  sapply(colnames(propVsN),function(x)lines(propCuts,propVsN[,x]+1,lwd=2,col=cols[x]))
  legend('topright',names(cols),col=cols,bty='n',lty=1,lwd=2)
  #legend('topright',names(cols),col=cols,bty='n',lty=1,lwd=2)
dev.off()
}

cols<-rainbow.lab(length(singles[names(singles)!='Air Swab']),alpha=.7,lightMultiple=.7,lightScale=0)
cols2<-rainbow.lab(length(singles[names(singles)!='Air Swab']),alpha=.3,lightMultiple=.9,lightScale=0)
names(cols)<-names(cols2)<-names(singles[names(singles)!='Air Swab'])
pdf('out/propVsN_single.pdf',width=4,height=4) 
  for(ii in c('MoBio','PSP')){
    thisSums<-singleSums[[ii]]
    thisTotals<-singleTotals[[ii]]
    par(mar=c(2.6,4.2,1.1,.2),las=1,lheight=.7,mgp=c(2.5,.8,0))
    plot(1,1,type='n',xlab='',ylab='Number of OTUs present in 2 or more samples\nand not in air swabs',ylim=range(thisSums+1),xlim=range(propCuts),log='xy',xaxt='n',yaxt='n',main=ii)
    title(xlab='Presence cutoff (proportion of sample)',mgp=c(1.5,1,0))
    prettyX<-pretty(log10(propCuts))
    axis(1,10^prettyX,sapply(prettyX,function(x)ifelse(x<0,as.expression(bquote(10^.(x))),10^x)),las=1,mgp=c(3,.5,0))
    prettyY<-c(0,1,5,10,50,100,500,1000)
    axis(2,prettyY+1,as.character(prettyY))
    ticks<-c(1:10,seq(20,100,10),seq(200,1000,100),seq(2000,10000,1000))+1
    axis(2,ticks,rep('',length(ticks)),tcl=-.25)
    sapply(colnames(thisSums),function(x)polygon(c(propCuts,rev(propCuts)),c(thisSums[,x],rev(thisTotals[,x]))+1,col=cols2[x],border=NA))#,border=cols[x]))
    sapply(colnames(thisSums),function(x)lines(propCuts,thisSums[,x]+1,lwd=2,col=cols[x]))
    legend('topright',names(cols),col=cols,bty='n',lty=1,lwd=2)
  }
dev.off()

propSelector<-propCuts<=.01
pdf('out/portion_single.pdf',width=8,height=4) 
  par(mfrow=c(1,2))
  for(ii in c('MoBio','PSP')){
    thisProps<-singleSums[[ii]]/singleTotals[[ii]]
    thisProps<-thisProps[propSelector,]
    par(mar=c(2.6,4.2,1.1,.2),las=1,lheight=.7,mgp=c(2.5,.8,0))
    plot(1,1,type='n',xlab='',ylab='Proportion of OTUs not shared\nwith negative controls',ylim=c(0,1),xlim=range(propCuts[propSelector]),log='x',xaxt='n',main=ii)
    title(xlab='Presence cutoff (proportion of sample)',mgp=c(1.5,1,0))
    prettyX<-pretty(log10(propCuts[propSelector]))
    axis(1,10^prettyX,sapply(prettyX,function(x)ifelse(x<0,as.expression(bquote(10^.(x))),10^x)),las=1,mgp=c(3,.5,0))
    sapply(colnames(thisProps),function(x){
      lines(propCuts[propSelector],thisProps[,x],lwd=2,col=cols[x])
      ci<-mapply(function(x,y)binom.test(x,y)$conf.int,singleSums[[ii]][propSelector,x],singleTotals[[ii]][propSelector,x])
      polygon(c(propCuts[propSelector],rev(propCuts[propSelector])),c(ci[1,],rev(ci[2,])),col=cols2[x],border=NA)
    })
    legend('topright',names(cols),col=cols,bty='n',lty=1,lwd=2)
  }
dev.off()


