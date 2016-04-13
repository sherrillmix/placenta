library(dnar)
library(parallel)
if(!exists('otus')){
  raw<-readLines('data/classic_otu_table.txt')
  raw[2]<-sub('^#','',raw[2])
  otus<-read.table(textConnection(raw),comment='#',sep='\t',header=TRUE,stringsAsFactors=FALSE)
  otus$taxa<-sapply(strsplit(otus$Consensus.Lineage,'; '),function(x){x<-x[!grepl('^[a-z]__$',x)];ifelse(grepl('^s__',tail(x,1)),paste(tail(x,2),collapse=' '),tail(x,1))})
  rownames(otus)<-ave(otus$taxa,otus$taxa,FUN=function(x)sprintf('%s #%d',x,1:length(x)))

  raw<-readLines('data/combined_mapping_qiime_keepers.txt')
  raw[1]<-sub('^#','',raw[1])
  samples<-read.table(textConnection(raw),sep='\t',header=TRUE,stringsAsFactors=FALSE)
  #samples$id<-ave(samples$SampleID,paste(samples$ExtractionKit,samples$SampleType),FUN=function(x)1:length(x))
  samples$id<-''
  samples[samples$SampleType!='OR Air Swab'&grepl('[567][0-9]',samples$SampleID),'id']<-sub('.*([567][0-9]).*','\\1',samples[samples$SampleType!='OR Air Swab'&grepl('[567][0-9]',samples$SampleID),'SampleID'])
  if(any(samples$id!=''&!grepl('^[567][0-9]$',samples$id)))stop(simpleError("Problem finding patient id"))
  samples$extract<-sub(' .*$','',samples$ExtractionKit)
  samples$type<-sub('OR ','',samples$SampleType)
  samples$name<-sprintf('%s %s %s%s',samples$extract,samples$type,ifelse(samples$id=='','','#'),samples$id)
  samples$sample<-ifelse(samples$SampleType %in% c('Extraction Blank','OR Air Swab','Sterile Swab'),'Negative control',samples$SampleType)
  samples$airSplit<-ifelse(samples$type=='Air Swab',samples$type,samples$sample)
  rownames(samples)<-samples$SampleID

  if(any(!samples$SampleID %in% colnames(otus)))stop(simpleError('Missing sample'))
}


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



