source("~/scripts/R/dna.R")
if(!exists('otus')){
  raw<-readLines('data/classic_otu_table.txt')
  raw[2]<-sub('^#','',raw[2])
  otus<-read.table(textConnection(raw),comment='#',sep='\t',header=TRUE,stringsAsFactors=FALSE)

  raw<-readLines('data/combined_mapping_qiime_keepers.txt')
  raw[1]<-sub('^#','',raw[1])
  samples<-read.table(textConnection(raw),sep='\t',header=TRUE,stringsAsFactors=FALSE)

  if(any(!samples$SampleID %in% colnames(otus)))stop(simpleError('Missing sample'))
}

otus$taxa<-sapply(strsplit(otus$Consensus.Lineage,'; '),function(x){x<-x[!grepl('^[a-z]__$',x)];ifelse(grepl('^s__',tail(x,1)),paste(tail(x,2),collapse=' '),tail(x,1))})
rownames(otus)<-ave(otus$taxa,otus$taxa,FUN=function(x)sprintf('%s #%d',x,1:length(x)))

otuN<-apply(otus[,samples$SampleID],1,sum)
mobioFs<-otus[,samples[samples$SampleType=='Placenta (FS)'&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
mobioMs<-otus[,samples[samples$SampleType=='Placenta (MS)'&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
pspFs<-otus[,samples[samples$SampleType=='Placenta (FS)'&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]
pspMs<-otus[,samples[samples$SampleType=='Placenta (MS)'&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]
mobioControl<-otus[,samples[samples$SampleType %in% c('OR Air Swab','Extraction Blank','Sterile Swab')&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
pspControl<-otus[,samples[samples$SampleType %in% c('OR Air Swab','Extraction Blank','Sterile Swab')&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]

mobioSaliva<-otus[,samples[samples$SampleType=='Saliva'&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
pspSaliva<-otus[,samples[samples$SampleType=='Saliva'&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]

otuProp<-apply(otus[,samples$SampleID],2,function(x)x/sum(x))
mobioControlProp<-apply(mobioControl,2,function(x)x/sum(x))
pspControlProp<-apply(pspControl,2,function(x)x/sum(x))
mobioFsProp<-apply(mobioFs,2,function(x)x/sum(x))
mobioMsProp<-apply(mobioMs,2,function(x)x/sum(x))
pspFsProp<-apply(pspFs,2,function(x)x/sum(x))
pspMsProp<-apply(pspMs,2,function(x)x/sum(x))
mobioSalivaProp<-apply(mobioSaliva,2,function(x)x/sum(x))
pspSalivaProp<-apply(pspSaliva,2,function(x)x/sum(x))

controlReadCut<-.0001
controlPropCut<-.2
inControlMobio<-apply(mobioControlProp,1,function(x)sum(x>controlReadCut))>ncol(mobioControl)*controlPropCut
inControlPsp<-apply(pspControlProp,1,function(x)sum(x>controlReadCut))>ncol(pspControl)*controlPropCut

sampleReadCut<-.0001
samplePropCut<-.3
inBothFs<-apply(mobioFsProp,1,function(x)sum(x>sampleReadCut))>ncol(mobioFsProp)*samplePropCut & apply(pspFsProp,1,function(x)sum(x>sampleReadCut))>ncol(pspFsProp)*samplePropCut & !inControlMobio & !inControlPsp
print(summary(inBothFs))
inBothMs<-apply(mobioMsProp,1,function(x)sum(x>sampleReadCut))>ncol(mobioMsProp)*samplePropCut & apply(pspMsProp,1,function(x)sum(x>sampleReadCut))>ncol(pspMsProp)*samplePropCut & !inControlMobio & !inControlPsp
print(summary(inBothMs))
inBothSaliva<-apply(mobioSalivaProp,1,function(x)sum(x>sampleReadCut))>ncol(mobioSalivaProp)*samplePropCut & apply(pspSalivaProp,1,function(x)sum(x>sampleReadCut))>ncol(pspSalivaProp)*samplePropCut & !inControlMobio & !inControlPsp
print(summary(inBothSaliva))

otuMeanScale<-t(apply(otuProp,1,function(x)x/max(x)))
samples$sample<-ifelse(samples$SampleType %in% c('Extraction Blank','OR Air Swab','Sterile Swab'),'Negative control',samples$SampleType)
sampleCols<-rainbow.lab(length(unique(samples$sample)),lightScale=0,lightMultiple=.7)
names(sampleCols)<-unique(samples$sample)
cols<-c('white',tail(rev(heat.colors(130)),99))
breaks<-c(-1,seq(1e-9,1+1e-9,length.out=100))
pdf('out/heatmap.pdf',height=8,width=8)
  #heatmap(t(as.matrix(otus[otuN>20,samples$SampleID])),col=rev(heat.colors(100)),margins=c(5,10),scale='column')
  #heatmap(t(otuMeanScale[otuN>100,]),col=c('white',tail(rev(heat.colors(130)),99)),margins=c(1,10),scale='none',breaks=breaks,labCol=rep('',nrow(otuMeanScale)))
  heatmap(t(otuMeanScale[apply(otuProp,1,max)>.001,]),col=cols,margins=c(.1,14),scale='none',breaks=breaks,labCol=rep('',nrow(otuMeanScale)),RowSideColors=sampleCols[samples$sample])#,Colv=hclust(dist(otuProp[apply(otuProp,1,max)>.001,]>1)))
  insetPos<-c(grconvertX(0.02,'nfc','user'),grconvertY(0.85,'nfc','user'),grconvertX(0.21,'nfc','user'),grconvertY(0.87,'nfc','user'))
  breakPos<-breaks*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-length(breakPos)],insetPos[2],breakPos[-1]+1e-3,insetPos[4],col=cols,xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyLabs<-pretty(0:1)
  prettyPos<-prettyLabs*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,prettyLabs,xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.1,"Proportion of max",xpd=NA,adj=c(.5,0))
dev.off()

sampleOrder<-samples$SampleID[order(samples$SampleType,sub('([FM])(etal|aternal).Side.Biopsy','\\1SB.',samples$SampleID))]

addInset<-function(cols,breaks,x1=.7,x2=.99,y1=.02,y2=.04){
  insetPos<-c(grconvertX(x1,'nfc','user'),grconvertY(y1,'nfc','user'),grconvertX(x2,'nfc','user'),grconvertY(y2,'nfc','user'))
  breakPos<-breaks*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-length(breakPos)],insetPos[2],breakPos[-1]+1e-3,insetPos[4],col=cols,xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyLabs<-pretty(0:1)
  prettyPos<-prettyLabs*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,prettyLabs,xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.1,"Proportion of max",xpd=NA,adj=c(.5,0))
}


pdf('out/inBothNotControl.pdf',height=10,width=10)
  #heatmap(t(as.matrix(otus[inBothMs,samples$SampleID])>1)*1,col=rev(heat.colors(100)),margins=c(20,10),scale='none')
  heatmap(t(otuMeanScale[inBothFs,sampleOrder]),col=cols,breaks=breaks,margins=c(14,14),scale='none',main='Fetal side',cexCol=1,Rowv=NA,Colv=NA)
  addInset(cols,breaks)
  heatmap(t(otuMeanScale[inBothMs,sampleOrder]),col=cols,breaks=breaks,margins=c(14,14),scale='none',main='Maternal side',cexCol=1,Rowv=NA,Colv=NA)
  addInset(cols,breaks)
  heatmap(t(otuMeanScale[inBothSaliva,sampleOrder]),col=cols,breaks=breaks,margins=c(14,14),scale='none',main='Saliva',cexCol=.5,Rowv=NA,Colv=NA)
  addInset(cols,breaks)
dev.off()
