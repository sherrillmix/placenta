bothExtract<-apply(table(samples$extract,samples$type)>1,2,all)
bothExtractTypes<-names(bothExtract)[bothExtract]

pairs<-lapply(bothExtractTypes,function(type){
#rownames(counts)<-names(pairs)
  thisPsp<-apply(otuProp[,samples[samples$type==type&samples$ExtractionKit=='PSP Spin Stool DNA Plus Kit','SampleID']]>0.0001,1,any)
  thisMobio<-apply(otuProp[,samples[samples$type==type&samples$ExtractionKit=='Mobio Powersoil','SampleID']]>0.0001,1,any)
  return(data.frame(psp=thisPsp,mobio=thisMobio))
})
names(pairs)<-bothExtractTypes

counts<-do.call(rbind,lapply(pairs,function(x){
  c(sum(x$psp&!x$mobio),sum(x$mobio&!x$psp),sum(x$psp&x$mobio))
}))
colnames(counts)<-c('PSP','MoBio','Both')


