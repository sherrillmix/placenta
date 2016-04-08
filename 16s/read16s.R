if(!exists('otus')){
  raw<-readLines('data/classic_otu_table.txt')
  raw[2]<-sub('^#','',raw[2])
  otus<-read.table(textConnection(raw),comment='#',sep='\t',header=TRUE,stringsAsFactors=FALSE)

  raw<-readLines('data/combined_mapping_qiime_keepers.txt')
  raw[1]<-sub('^#','',raw[1])
  samples<-read.table(textConnection(raw),sep='\t',header=TRUE,stringsAsFactors=FALSE)

  if(any(!samples$SampleID %in% colnames(otus)))stop(simpleError('Missing sample'))
}

otuN<-apply(otus[,samples$SampleID],1,sum)
mobioFs<-otus[,samples[samples$SampleType=='Placenta (FS)'&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
mobioMs<-otus[,samples[samples$SampleType=='Placenta (MS)'&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
pspFs<-otus[,samples[samples$SampleType=='Placenta (FS)'&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]
pspMs<-otus[,samples[samples$SampleType=='Placenta (MS)'&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]
mobioControl<-otus[,samples[samples$SampleType %in% c('OR Air Swab','Extraction Blank','Sterile Swab')&samples$ExtractionKit=='Mobio Powersoil','SampleID']]
pspControl<-otus[,samples[samples$SampleType %in% c('OR Air Swab','Extraction Blank','Sterile Swab')&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]
otus$taxa<-sapply(strsplit(otus$Consensus.Lineage,'; '),function(x){x<-x[!grepl('^[a-z]__$',x)];ifelse(grepl('^s__',tail(x,1)),paste(tail(x,2),collapse=' '),tail(x,1))})
rownames(otus)<-ave(otus$taxa,otus$taxa,FUN=function(x)sprintf('%s #%d',x,1:length(x)))

controlReadCut<-10
controlPropCut<-.3
inControlMobio<-apply(mobioControl,1,function(x)sum(x>controlReadCut))>ncol(mobioControl)*controlPropCut
inControlPsp<-apply(pspControl,1,function(x)sum(x>controlReadCut))>ncol(pspControl)*controlPropCut

sampleReadCut<-5
samplePropCut<-.3
inBothFs<-apply(mobioFs,1,function(x)sum(x>sampleReadCut))>ncol(mobioFs)*samplePropCut & apply(pspFs,1,function(x)sum(x>sampleReadCut))>ncol(pspFs)*samplePropCut & !inControlMobio & !inControlPsp
print(summary(inBothFs))
inBothMs<-apply(mobioMs,1,function(x)sum(x>sampleReadCut))>ncol(mobioMs)*samplePropCut & apply(pspMs,1,function(x)sum(x>sampleReadCut))>ncol(pspMs)*samplePropCut & !inControlMobio & !inControlPsp
print(summary(inBothMs))

pdf('test.pdf',height=12,width=12)
heatmap(t(as.matrix(otus[otuN>20,samples$SampleID])),col=rev(heat.colors(100)),margins=c(5,10),scale='column')
dev.off()

pdf('test.pdf',height=12,width=12)
  #heatmap(t(as.matrix(otus[inBothMs,samples$SampleID])>1)*1,col=rev(heat.colors(100)),margins=c(20,10),scale='none')
  heatmap(log(t(as.matrix(otus[inBothMs|inBothFs,samples$SampleID]))+1),col=rev(heat.colors(100)),margins=c(20,10),scale='column')
dev.off()
