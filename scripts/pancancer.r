library(tools)
library(ggfortify)
library(reshape)
library(grid)
library(rms)
#library(ggbio)
#library(GenomicRanges)
library(RColorBrewer)
library(colorspace)
library(gridExtra)
library(reshape2)
library(fpc)
library(ggplot2)
require(scales)
library(VennDiagram)
library(fitdistrplus)
library(stats)
library(LaplacesDemon)
library(survival)
library(survcomp)
library(dplyr)
library(OIsurv) # Aumatically loads KMsurv
library(ranger)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}


geometric_mean<-function(x) {
  n = length(x)
  return(prod(x)^(1/n))
}

geometric_mean_c<-function(x) {
  n = length(x)
  return(((exp(1)^(sum(log(1+x))))^(1/n))-1)
}

geometric_sd <- function(x, na.rm = FALSE, ...) {
  exp(sd(log(x, ...), na.rm = na.rm, ...))
}

heat<-function() {
  cnv<-read.csv('/pedigree2/projects/namphuon/programs/pancancer/analyses/tumor/TCGA-CESC/min_cnv1/TCGA-C5-A0TN-01A-21_WGS167271038740_CNV14267_mincnv1_20170905-133947/TCGA-C5-A0TN.cnv.bed',sep='\t')
    pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/cnv_dist.pdf",sep=""),width=3,height=4)
      p1 <- ggplot(cnv, aes(x=log(cnv$End-cnv$Start)/log(10)))+          
          geom_histogram()+
          geom_vline(xintercept=log(10000)/log(10))+
          geom_vline(xintercept=log(1000)/log(10))+
          geom_vline(xintercept=log(100)/log(10))+
          ylab("Counts")
      print(p1)
      dev.off()  

  
  
  meta<-read.csv('/pedigree2/projects/namphuon/programs/pancancer/data/heatmaps/meta.csv',sep=',',header=FALSE)
  names(meta)<-c('Sample','SampleID','Disease','Type')
  meta<-meta[order(meta$Type,meta$Disease),]
  meta$SampleID<-paste(meta$SampleID)
  
  aa_data<-read.csv('/pedigree2/projects/namphuon/programs/pancancer/data/heatmaps/aa_matrix.csv',sep=',',header=TRUE)
  aa_data$Sample<-gsub("-", "_", aa_data$Sample)
  rownames(aa_data)<-aa_data$Sample
  colnames(aa_data)<-names(aa_data)
  cnv_data<-read.csv('/pedigree2/projects/namphuon/programs/pancancer/data/heatmaps/cnv_matrix.csv',sep=',',header=TRUE)
  cnv_data$Sample<-gsub("-", "_", cnv_data$Sample)
  rownames(cnv_data)<-cnv_data$Sample
  colnames(cnv_data)<-names(cnv_data)
  
  aa_matrix<-data.matrix(aa_data[meta$SampleID,meta$SampleID])
  aa_matrix[aa_matrix==-1]<-0
  aa_min_max<-aa_matrix-min(aa_matrix)/(max(aa_matrix)-min(aa_matrix))
  cnv_matrix<-data.matrix(cnv_data[meta$SampleID,meta$SampleID])
  cnv_matrix[cnv_matrix==-1]<-0
  
  

data['GBM2','GBM1']
GBM1=TCGA-02-2483-01A-01_WGS185258413619_CNV17405_mincnv1_20170905-133947
GBM2=TCGA-02-2485-01A-01_WGS187013547986_CNV17480_mincnv1_20170905-133947
}
  stad_data<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/metadata/TCGA-STAD.study.csv',sep=''),header=TRUE)
data<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/ecDNA.csv',sep=''),header=TRUE)
names(stad_data)[1]<-'Sample'
merged<-merge(data[data$Disease=='TCGA-STAD',],stad_data,by=c('Sample'))
table(merged$CIMP.Category)
table(merged[merged$Tumor=='ecDNA',]$CIMP.Category)
#68% other, 82% other for ecDNA

(table(merged$Methylation.Cluster))/sum(table(merged$Methylation.Cluster))
(table(merged[merged$Tumor=='ecDNA',]$Methylation.Cluster))/sum(table(merged[merged$Tumor=='ecDNA',]$Methylation.Cluster))

#54% high, 82% high
(table(merged$Copy.Number.Cluster))/sum(table(merged$Copy.Number.Cluster))
(table(merged[merged$Tumor=='ecDNA',]$Copy.Number.Cluster))/sum(table(merged[merged$Tumor=='ecDNA',]$Copy.Number.Cluster))

#C3/C4 gene expression category
(table(merged$Gene.Expression.Cluster))/sum(table(merged$Gene.Expression.Cluster))
(table(merged[merged$Tumor=='ecDNA',]$Gene.Expression.Cluster))/sum(table(merged[merged$Tumor=='ecDNA',]$Gene.Expression.Cluster))

(table(merged$WHO.Class))/sum(table(merged$WHO.Class))
(table(merged[merged$Tumor=='ecDNA',]$WHO.Class))/sum(table(merged[merged$Tumor=='ecDNA',]$WHO.Class))


mean(merged[!is.na(merged$Total.Mutation.Rate),]$Total.Mutation.Rate)
mean(merged[merged$Tumor=='ecDNA' & !is.na(merged$Total.Mutation.Rate),]$Total.Mutation.Rate)
}


historgram_ecDNA<-function() {

  data<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/ecDNA.csv',sep=''),header=TRUE)
  
  melted<-melt(data[,c('Sample','Tlen','Nlen')],id=c('Sample'))    
  melted$Name<-'Tumor'
  melted[melted$variable=='Nlen',]$Name<-'Normal'
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/cycle_length.pdf",sep=""),width=3,height=4)
      p1 <- ggplot(melted[melted$value != 0,], aes(x=value))+
          geom_histogram()+
          xlab('Cycle length')+
          scale_x_log10(limits=c(1e0, 1e8), breaks=10^seq(0,8))+
          facet_wrap(~Name,nrow=2)+
          geom_vline(xintercept = 5000, size = 1, colour = "green", linetype = "dashed") +          
          geom_vline(xintercept = 10000, size = 1, colour = "blue", linetype = "dashed") +
          geom_vline(xintercept = 50000, size = 1, colour = "red", linetype = "dashed") +
          theme_bw()+     
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),  axis.line = element_line(colour = "black"))+                 
          ylab("Counts")
      print(p1)
      dev.off()  

  melted<-melt(data[,c('Sample','Disease','Tumor','Normal')],id=c('Sample','Disease'))
  melted$Count<-1
  m<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)  
  orders<-m[order(-m$ecDNA/(m$ecDNA+m$'No Amplified'+m$'No cyclic')),]$Disease
  for (len in c(5000,10000,25000,50000)) {
    subdata<-data
    subdata[subdata$Tlen < len & subdata$Tumor=='ecDNA',]$Tumor<-'No cyclic'
    if (dim(subdata[subdata$Nlen < len & subdata$Normal=='ecDNA',])[1] > 0) {
      subdata[subdata$Nlen < len & subdata$Normal=='ecDNA',]$Normal<-'No cyclic'
    }
    subdata<-subdata[,c('Sample','Disease','Tumor','Normal')]  
    melted<-melt(subdata,id=c('Sample','Disease'))
    melted$Count<-1
    m<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)
    melted$Disease<-factor(melted$Disease,levels=m[order(-m$ecDNA/(m$ecDNA+m$'No Amplified'+m$'No cyclic')),]$Disease)
    melted$value<-factor(melted$value,levels=c('ecDNA','No cyclic','No Amplified'))
    
    t<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)
    t$Disease<-factor(t$Disease,levels=t[order(-t$ecDNA/(t$ecDNA+m$'No Amplified'+t$'No cyclic')),]$Disease)
    t[,c('ecDNA','No Amplified','No cyclic')]<-t[,c('ecDNA','No Amplified','No cyclic')]/rowSums(t[,c('ecDNA','No Amplified','No cyclic')])
    t$variable<-'Tumor'
    t$variable<-factor(t$variable,levels=c('Tumor','Normal'))
  
    n<-cast(melted[melted$variable=='Normal',c('Disease','Count','value')],Disease~value)
    n[,c('ecDNA','No Amplified','No cyclic')]<-n[,c('ecDNA','No Amplified','No cyclic')]/rowSums(n[,c('ecDNA','No Amplified','No cyclic')])
    n$variable<-'Normal'
    combined<-rbind(t,n)  
    names(combined)<-c('Disease','ecDNA','No cyclic','No Amplified','Type')
    combined$Type<-factor(combined$Type,levels=c('Normal','Tumor'))
    melted<-melt(data.frame(combined),id.vars=c('Disease','Type'))
    melted$Disease<-factor(melted$Disease,levels=orders)
    melted$variable<-factor(melted$variable, levels=c('No.Amplified','No.cyclic','ecDNA'))
    pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/ecDNA_dist_norm.",len,".pdf",sep=""),width=5.5,height=5.5)
        p1 <- ggplot(melted, aes(x=Disease,fill=variable,y=value))+
            geom_bar(stat='identity')+
            scale_fill_brewer(palette='RdYlBu',name='',labels=c("Cyclic","No cyclic","No Amplified"),breaks=c("ecDNA","No.cyclic","No.Amplified"))+   
            #scale_fill_brewer(palette='RdYlBu',name='',labels=c("No Amplified","No cyclic","Cyclic"),breaks=c("No.Amplified","No.cyclic","ecDNA"))+            
            facet_wrap(~Type,nrow=2)+
            theme_bw()+     
            theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
            ylab("Frequency")
        print(p1)
        dev.off()  
    
  }
  
  types = c('Colon','Haemat','Prostate','Breast','Ovarian','Melan.','Renal','Lung','GBM')
  t1<-data[data$Disease=='TCGA-COAD',]
  t1$Name='Colon'
  t2<-data[data$Disease %in% c('TCGA-LAML'),]
  t2$Name='Haemat'
  t3<-data[data$Disease %in% c('TCGA-PRAD'),]
  t3$Name='Prostate'
  t4<-data[data$Disease %in% c('TCGA-BRCA'),]
  t4$Name='Breast'
  t5<-data[data$Disease %in% c('TCGA-OV'),]
  t5$Name='Ovarian'
  t6<-data[data$Disease %in% c('TCGA-SKCM'),]
  t6$Name='Melan.'
  t7<-data[data$Disease %in% c('TCGA-KIRP','TCGA-KIRC','TCGA-KICH'),]
  t7$Name='Renal'
  t8<-data[data$Disease %in% c('TCGA-LUSC','TCGA-LUAD'),]
  t8$Name='Lung'
  t9<-data[data$Disease %in% c('TCGA-GBM'),]
  t9$Name='GBM'
  res<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9)
  len=10000
  res[res$Tlen < len & res$Tumor=='ecDNA',]$Tumor<-'No cyclic'
  if (dim(res[res$Nlen < len & res$Normal=='ecDNA',])[1] > 0) {
    res[res$Nlen < len & res$Normal=='ecDNA',]$Normal<-'No cyclic'
  }
  res<-res[,c('Sample','Name','Tumor','Normal')]  
  names(res)[2]<-'Disease'
  res[res$Tumor != 'ecDNA',]$Tumor<-'No cyclic'
  res[res$Normal != 'ecDNA',]$Normal<-'No cyclic'
  melted<-melt(res,id=c('Sample','Disease'))
  melted$Count<-1
  m<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)
  melted$Disease<-factor(melted$Disease,levels=types)
  melted$value<-factor(melted$value,levels=c('ecDNA','No cyclic','No Amplified'))
  
  t<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)
  t$Disease<-factor(t$Disease,levels=types)
  t[,c('ecDNA','No cyclic')]<-t[,c('ecDNA','No cyclic')]/rowSums(t[,c('ecDNA','No cyclic')])
  t$variable<-'Tumor'
  t$variable<-factor(t$variable,levels=c('Tumor','Normal'))

  n<-cast(melted[melted$variable=='Normal',c('Disease','Count','value')],Disease~value)
  n[,c('ecDNA','No cyclic')]<-n[,c('ecDNA','No cyclic')]/rowSums(n[,c('ecDNA','No cyclic')])
  n$variable<-'Normal'
  combined<-rbind(t,n)  
  names(combined)<-c('Disease','ecDNA','No Amplified','Type')
  combined$Type<-factor(combined$Type,levels=c('Normal','Tumor'))
  melted<-melt(data.frame(combined),id.vars=c('Disease','Type'))
  melted$Disease<-factor(melted$Disease,levels=types)
  melted$variable<-factor(melted$variable,levels=c('No.Amplified','ecDNA'))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/ecDNA_dist_norm.turner.",len,".pdf",sep=""),width=5,height=3)
      p1 <- ggplot(melted[melted$Type=='Tumor',], aes(x=Disease,fill=variable,y=value))+
          geom_bar(stat='identity')+
          scale_fill_manual(values = c("blue1","red1"),name='', labels=c('No cyclic', 'Cyclic'), breaks=c('No.Amplified','ecDNA'))+
#          scale_fill_brewer(palette='RdYlBu')+            
#          facet_wrap(~Type,nrow=2)+
          theme_bw()+     
          theme(legend.position='top',legend.direction='horizontal', panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          ylab("Frequency")
      print(p1)
      dev.off()  
  
  
  
    
  melted<-melt(data,id=c('Sample','Disease'))
  melted$Count<-1
  m<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)
  melted$Disease<-factor(melted$Disease,levels=m[order(-m$ecDNA/(m$ecDNA+m$'No Amplified'+m$'No cyclic')),]$Disease)
  melted$value<-factor(melted$value,levels=c('ecDNA','No cyclic','No Amplified'))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/ecDNA_dist.pdf",sep=""),width=6,height=8)
      p1 <- ggplot(melted, aes(x=Disease))+
          geom_bar(aes(fill=value))+
          facet_wrap(~variable,nrow=2)+
          theme_bw()+     
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          ylab("Counts")
      print(p1)
      dev.off()  

  t<-cast(melted[melted$variable=='Tumor',c('Disease','Count','value')],Disease~value)
  t$Disease<-factor(t$Disease,levels=t[order(-t$ecDNA/(t$ecDNA+m$'No Amplified'+t$'No cyclic')),]$Disease)
  t[,c('ecDNA','No Amplified','No cyclic')]<-t[,c('ecDNA','No Amplified','No cyclic')]/rowSums(t[,c('ecDNA','No Amplified','No cyclic')])
  t$variable<-'Tumor'
  t$variable<-factor(t$variable,levels=c('Tumor','Normal'))
  
  n<-cast(melted[melted$variable=='Normal',c('Disease','Count','value')],Disease~value)
  n[,c('ecDNA','No Amplified','No cyclic')]<-n[,c('ecDNA','No Amplified','No cyclic')]/rowSums(n[,c('ecDNA','No Amplified','No cyclic')])
  n$variable<-'Normal'
  combined<-rbind(t,n)  
  names(combined)<-c('Disease','ecDNA','No cyclic','No Amplified','Type')
  combined$Type<-factor(combined$Type,levels=c('Tumor','Normal'))
  melted<-melt(data.frame(combined),id.vars=c('Disease','Type'))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/ecDNA_dist_norm.pdf",sep=""),width=6,height=8)
      p1 <- ggplot(melted, aes(x=Disease,fill=variable,y=value))+
          geom_bar(stat='identity')+
          facet_wrap(~Type,nrow=2)+
          theme_bw()+     
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          ylab("Frequency")
      print(p1)
      dev.off()  
  
}

intervals<-function() {
    intervals <-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/intervals.csv',sep=''),header=TRUE)
    intervals$NFPKM<-intervals$FPKM/intervals$CNV
    tmp<-intervals[intervals$IsAmplicon=='True',]
    tmp$amplicon_name<-paste(tmp$Patient,tmp$Amplicon)
    subint<-merge(intervals,tmp[,c('Interval','amplicon_name')])
    casted<-cast(melt(subint[,c('Patient','amplicon_name','NFPKM')]),fun=sum)
    tmp<-cast(melt(data.frame(casted[,c('amplicon_name','NFPKM')])),fun=mean)
    avint<-tmp[tmp$NFPKM != 0,]
    names(avint)<-c('amplicon_name','AvFPKM')
    subint<-subint[subint$amplicon_name %in% avint$amplicon_name & subint$IsAmplicon=='True',]   
    subint$num_cycles<-0
    subint$num_amplicons<-0
    subint[subint$IsCycle=='True',]$num_cycles<-1    
    subint[subint$PartAmplicon=='True' | subint$IsAmplicon == 'True',]$num_amplicons<-1    
    merged<-merge(subint,avint,by=c('amplicon_name'))
    casted<-cast(melt(merged[,c('amplicon_name', 'num_cycles','num_amplicons','Disease','NFPKM','AvFPKM')]),fun=sum)
    casted$fold_change<-(casted$NFPKM+1)/(casted$AvFPKM+1)
    casted$IsCycle<-'False'
    casted[casted$num_cycles>0,]$IsCycle<-'True'
    
    for (disease in unique(casted$Disease)){ 
      if (dim(casted[casted$Disease==disease & casted$IsCycle=='True',])[1] > 10)  {        
        res = ks.test(casted[casted$Disease==disease & casted$IsCycle=='True',]$fold_change,casted[casted$Disease==disease &  casted$IsCycle =='False',]$fold_change,alternative='less')
        print(paste(disease,res$p.value))
      }
    }

pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/amplicon.cycle.pdf",sep=""),width=8,height=4)
      p1 <- ggplot(casted, aes(x=Disease, y=fold_change, fill=IsCycle))+
          geom_boxplot()+
          scale_y_continuous(trans="log2", labels = scientific)+ 
          geom_hline(yintercept=1)+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")+
          ylab("FPKM Fold-change")
      print(p1)
      dev.off()  
      
      t<-melt(data.frame(casted[,c('Disease','IsCycle')]))
      t$Count<-1
      t<-cast(t,Disease~IsCycle)
      t$percent<-t$True/(t$True+t$False)      pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/amplicon.bar.pdf",sep=""),width=8,height=4)
      p1 <- ggplot(casted, aes(x=Disease, fill=IsCycle))+
          geom_bar()+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")+
          ylab("Counts")
      print(p1)
      dev.off()  
      
      
  
    tmp<-cast(melt(intervals[intervals$IsAmplicon == 'False' & intervals$PartAmplicon == 'False', c('Interval','NFPKM')]),fun=mean)
    names(tmp)<-c('Interval','AvFPKM')    
    keep = tmp[tmp$AvFPKM != 0,]            
    subint<-merge(intervals,keep)
    subint$Fold<-(subint$NFPKM+1)/(subint$AvFPKM+1)
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/intervals.pdf",sep=""),width=8,height=4)
      p1 <- ggplot(subint[subint$IsAmplicon=='True',], aes(x=Disease, y=Fold, fill=IsCycle))+
          geom_boxplot()+
          scale_y_continuous(trans="log2", labels = scientific)+ 
          geom_hline(yintercept=1)+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")+
          ylab("FPKM Fold-change")
      print(p1)
      dev.off()  

    for (disease in unique(casted$Disease)){ 
      if (dim(casted[casted$Disease==disease & casted$IsCycle=='True',])[1] > 5)  {        
        res = ks.test(subint[subint$IsAmplicon=='True' & subint$Disease==disease & subint$IsCycle=='True',]$Fold,subint[subint$IsAmplicon=='True' & subint$Disease==disease &  subint$IsCycle =='False',]$Fold,alternative='less')
        print(paste(disease,res$p.value))
      } else {
        print(paste("Can't test",disease))
      }
    }      
}

segments<-fuction() {
  graph<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/aa_overview.csv',sep=''))
  graph$amplicon_url<-NULL
  
  subgraph<-graph[graph$cyclic_cycles > 0 & !is.na(graph$cyclic_cycles) & graph$average_amplified_copy_count > 5,]
  head(subgraph[order(subgraph$total_amplified_amplicon_size)])
  
  subgraph<-graph[graph$cyclic_cycles > 0 & !is.na(graph$cyclic_cycles),]
  #Counts total amplicon per disease type
  counts_disease<-table(graph$Disease)
  counts_subgraph<-table(subgraph$Disease)
  cylic_amplicons_per_sample<-counts_subgraph/counts_disease  
  graph$cylic<-'False'
  graph[graph$cyclic_cycles > 0 & !is.na(graph$cyclic_cycles),]$cylic<-'True'
  graph$amplicon_url<-NULL
  graph<-graph[graph$total_amplified_amplicon_size > 5000 ,]
  
  
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/distribution_cyclic_amplified_amplicon_size.pdf",sep=""),width=8,height=8)
      p1 <- ggplot(graph, aes(x=cylic, y=total_amplified_amplicon_size))+
          geom_boxplot()+
          scale_y_continuous(trans="log10", labels = scientific)+ 
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")+
          ylab("Total Amplified Amplicon Size")
      print(p1)
      dev.off()  
  
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/distribution_cyclic_amplicon_size.pdf",sep=""),width=8,height=8)
      p1 <- ggplot(graph, aes(x=cylic, y=total_amplicon_size))+
          geom_boxplot()+
          scale_y_continuous(trans="log10", labels = scientific)+ 
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")+
          ylab("Total Amplicon Size")
      print(p1)
      dev.off()  
  
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/distribution_cyclic_breakpoint.pdf",sep=""),width=8,height=8)
      p1 <- ggplot(graph, aes(x=Disease, y=total_coverage_shifts_with_breakpoint,fill=cylic))+
          geom_boxplot()+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")+
          ylab("Total Coverage Shifts With Breakpoint")
      print(p1)
      dev.off()  
  
    pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/recovered_intervals.pdf",sep=""),width=16,height=16)
      p1 <- ggplot(graph, aes(x=Size, y=Amplification))+
          geom_point()+
          scale_x_continuous(trans="log10", labels = scientific)+  
          facet_wrap(~Reconstructed)+
          scale_colour_gradient(low = "blue", high = "red")+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Size of CNV")+
          ylab("Amplification")
      print(p1)
      dev.off()  

}

ecdna_graph<-function() {
  graph<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/graph.csv',sep=''))
  
  filtered_short<-graph[graph$total_length < 50000,]
  filtered_long<-graph[graph$total_length >= 50000,]
  
  melted<-melt(filtered_long[,c('Sample', 'total_copy_count','total_segment_count','total_length','total_chrom','total_breakpoints','size_of_heaviest_path','size_of_longest_path')])
  casted<-melted
  casted$Sample<-NULL
  casted<-merge(melt(cast(casted,fun=mean)),melt(cast(casted,fun=sd)),by=c('variable'))
  for (var in c('Sample', 'total_copy_count','total_segment_count','total_length','total_chrom','total_breakpoints','size_of_heaviest_path','size_of_longest_path')) {
    high = .75
    low = 1-high
    iqr=quantile(filtered_long[,var], high) - quantile(filtered_long[,var], low)
    keep <- filtered_long[filtered_long[,var] >= quantile(filtered_long[,var], low) - 1.5*iqr & 
        filtered_long[,var] <= quantile(filtered_long[,var], high) + 1.5*iqr,]
    removed <- filtered_long[!(filtered_long[,var] >= quantile(filtered_long[,var], low) - 1.5*iqr & 
        filtered_long[,var] <= quantile(filtered_long[,var], high) + 1.5*iqr),]
    removed[order(removed[,var]),c('Sample','Amplicon_id','total_breakpoints',var)]    
  }
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/graph_distributions.pdf",sep=""),width=16,height=16)
      p1 <- ggplot(melted, aes(x=value+1))+
          geom_histogram()+
          scale_x_continuous(trans="log10", labels = scientific)+  
          facet_wrap(~variable,scales='free')+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          ylab("Counts")
      print(p1)
      dev.off()  
  
  
  
  subgraph<-graph[graph$total_length > 100000,]
  subgraph[subgraph$total_breakpoints == 0 & subgraph$total_circular_cycles >= 1,]$total_breakpoints = 1
  subgraph$type<-paste(subgraph$type)
    
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/graph_breakpoint.pdf",sep=""),width=16,height=16)
      p1 <- ggplot(subgraph[subgraph$total_breakpoints > 1 | subgraph$total_circular_cycles >= 1,], aes(x=total_length, y=total_cycle_count,size=log(total_breakpoints+1),shape=oncogenes>0))+
          geom_point(aes(colour = total_chrom))+
          scale_x_continuous(trans="log10", labels = scientific)+  
          scale_y_continuous(trans="log10", labels = scientific)+   
          facet_wrap(~Disease)+
          scale_colour_gradient(low = "blue", high = "red")+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Total Size of Amplicon")+
          ylab("Total Cycle Counts")
      print(p1)
      dev.off()  
      
  subgraph$has_cycle<-'False'
  subgraph[subgraph$total_circular_cycles > 0,]$has_cycle<-'True'
  subgraph[subgraph$oncogenes > 0,]$oncogenes<-1
  
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/graph_cycles.pdf",sep=""),width=8,height=72)
      p1 <- ggplot(subgraph, aes(x=total_length, y=total_cycle_count,size=log(total_breakpoints+1),shape=type))+
          geom_point(aes(colour = total_chrom))+
          scale_x_continuous(trans="log10", labels = scientific)+  
          scale_y_continuous(trans="log10", labels = scientific)+   
          facet_grid(Disease~oncogenes)+
          scale_colour_gradient(low = "blue", high = "red")+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Total Size of Amplicon")+
          ylab("Total Cycle Counts")
      print(p1)
      dev.off()  
      
      
  smaller<-graph[graph$total_length < 100000,]
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/graph_breakpoint.smaller.pdf",sep=""),width=16,height=16)
      p1 <- ggplot(smaller, aes(x=total_length, y=total_cycle_count+1,size=log(total_breakpoints+1),shape=oncogenes>0))+
          geom_point(aes(colour = total_chrom))+
          scale_x_continuous(trans="log10", labels = scientific)+  
          scale_y_continuous(trans="log10", labels = scientific)+  
          facet_wrap(~Disease)+
          scale_colour_gradient(low = "blue", high = "red")+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Total Size of Amplicon")+
          ylab("Total Cycle Counts")
      print(p1)
      dev.off()    
  
}

pancancer_two<-function() {

  fusion<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/fusion.csv',sep=''))
  fusion$Count<-1
  fused<-cast(melt(fusion[,c('Disease','Fusion','Count')]))
  fused[fused$Fusion=='False']
  
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/fusion_amplicon.pdf",sep=""),width=8,height=4)
      p1 <- ggplot(fused, aes(x=Disease, y=Count,fill=Fusion))+
          geom_col()+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")
      print(p1)
      dev.off()  
  
  fused<-merge(fused, cast(melt(data.frame(fused[,c('Disease','Count')])),fun=sum), by=c('Disease'))
  names(fused)<-c('Disease','Fusion','Count','Total')
  fused$Percent<-fused$Count/fused$Total*100
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/fusion_amplicon.percent.pdf",sep=""),width=8,height=4)
      p1 <- ggplot(fused, aes(x=Disease, y=Percent,fill=Fusion))+
          geom_col()+
          theme_bw() +             
          theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
          xlab("Disease")
      print(p1)
      dev.off()    
  
  disease = 'ALL'
  tpm<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/gene_amp.csv',sep=''))
  tpm$wFPKM<-tpm$FPKM/tpm$CNV    
  tpm$keeper<-paste(tpm$Disease,tpm$Gene)
  temp<-cast(melt(tpm[tpm$Structure=='None',c('Disease', 'wFPKM','Gene'),]),fun=mean)
  names(temp)<-c('Disease', 'Gene','AvwFPKM')
  merged<-merge(tpm,temp,by=c('Disease','Gene'))
  
  temp<-cast(melt(tpm[tpm$Structure=='None',c('Disease', 'FPKM','Gene'),]),fun=mean)
  names(temp)<-c('Disease', 'Gene','AvFPKM')
  merged<-merge(merged,temp,by=c('Disease','Gene'))
  
  
  #Remove genes that have no expression in more than 50% of samples  
  keeper<-merged[,c('keeper','FPKM')]
  keeper[keeper$FPKM > 0,]$FPKM<-1
  keeper$Total<-1
  t<-cast(melt(data.frame(keeper)),fun=sum)
  subt<-t[t$FPKM/t$Total> 0.50,]
  subtpm<-merged[merged$keeper %in% subt$keeper,]
  subtpm$fold_change<-(1+subtpm$FPKM)/(1+subtpm$AvFPKM)
  casted<-cast(melt(subtpm[,c('Disease','Gene','Structure','Oncogene','fold_change')]),fun=geometric_mean)
  casted[casted$fold_change==0,]$fold_change<- 1e-10

  
  #Take only genes that have all three
  keeper<-unique(subtpm[,c('keeper','Structure')])
  keeper$Total<-1
  keeper<-cast(melt(keeper[,c('keeper','Total')]))
  only_three<-keeper[keeper$Total>=3,]
  complete<-subtpm[subtpm$keeper %in% keeper[keeper$Total>=2,]$keeper,]
  complete$wfold_change<-(1+complete$wFPKM)/(1+complete$AvwFPKM)
  complete$fold_change<-(1+complete$FPKM)/(1+complete$AvFPKM)
  complete[complete$Structure=='Cycle',]$Structure='Amplicon'  
    casted<-cast(melt(complete[,c('Disease','Gene','Structure','Oncogene','fold_change','wfold_change')]),fun=geometric_mean)
    
  data_annotation = c()
  for (disease in unique(casted$Disease)) {
    for (oncogene in c('True','False')) {
      d1<-casted[casted$Disease==disease & casted$Oncogene==oncogene & casted$Structure=='None',c('fold_change','Gene')]
      d2<-casted[casted$Disease==disease &  casted$Oncogene==oncogene  & casted$Structure=='Amplicon',c('fold_change','Gene')]
      if (dim(d2)[1] < 10) {
        data_annotation<-rbind(data_annotation,c(disease,oncogene,as.character( 1),10))              
        next
      }
      names(d2)<-c('fold_change_amp','Gene')
      big<-max(max(d1$fold_change_amp,d2$fold_change))      
      m1<-merge(d1,d2,by=c('Gene'))
      r<-wilcox.test(m1$fold_change_amp,m1$fold_change,alternative='greater', paired=TRUE)
      data_annotation<-rbind(data_annotation,c(disease,oncogene,as.character( r$p.value),big))      
    }
  } 
  
  ks_annotation = c()
  for (disease in unique(casted$Disease)) {
    d1<-casted[casted$Disease==disease & casted$Structure=="Amplicon" & casted$Oncogene=="True" ,c('fold_change','Gene')]
    d2<-casted[casted$Disease==disease & casted$Structure=="Amplicon" & casted$Oncogene=="False" ,c('fold_change','Gene')]
    if (dim(d2)[1] < 10 | dim(d1)[1] < 10) {
      ks_annotation<-rbind(ks_annotation,c(disease,oncogene,as.character( 1),10))              
      next
    }
    names(d2)<-c('fold_change_amp','Gene')
    big<-max(max(d1$fold_change_amp,d2$fold_change))      
    r<-ks.test(d1$fold_change,d2$fold_change)
    ks_annotation<-rbind(ks_annotation,c(disease,oncogene,as.character( r$p.value),big))      
  } 
  
  data_annotation<-data.frame(data_annotation) 
  names(data_annotation)<-c('Disease','Oncogene','Sign','wfold_change')
  data_annotation$Structure<-'Amplicon'
  data_annotation$fold_change<-as.numeric(paste(data_annotation$wfold_change))
  data_annotation$Sign<-as.numeric(paste(data_annotation$Sign))  
  data_annotation$Label<-'ns'
  data_annotation$fold_change = 1
  data_annotation[data_annotation$Sign == 10,]$Label<-'N/A'
  data_annotation[data_annotation$Sign <= 0.05,]$Label<-'*'
  data_annotation[data_annotation$Sign <= 0.01,]$Label<-'**'
  data_annotation[data_annotation$Sign <= 0.001,]$Label<-'***'
  data_annotation[data_annotation$Sign <= 0.0001,]$Label<-'****'
  casted<-casted[casted$fold_change!=0,]
    
 pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_complete_wbox.pdf",sep=""),width=14,height=8)
    p1 <- ggplot(casted[casted$Disease != 'TCGA-DLBC',], aes(x=Oncogene, y=wfold_change,fill=Structure))+
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        geom_boxplot(outlier.shape = NA)+
        geom_text(data=data_annotation[data_annotation$Disease != 'TCGA-DLBC',], aes(label=Label,size=4,x=Oncogene),vjust=-0.75, position = position_dodge(0.9))+
        facet_wrap(~Disease,scales='free_y')+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="log2", labels = scientific)+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Oncogene")
    print(p1)
    dev.off()
    
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_complete_box.disease.pdf",sep=""),width=8,height=4)
    p1 <- ggplot(casted[casted$Disease %in% c('TCGA-GBM','TCGA-LGG','TCGA-BRCA', 'TCGA-LUSC','TCGA-LUAD'),], aes(x=Oncogene, y=fold_change,fill=Structure))+
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        geom_boxplot(outlier.shape = NA)+
        geom_text(data=data_annotation[data_annotation$Disease %in% c('TCGA-GBM','TCGA-LGG','TCGA-BRCA', 'TCGA-LUSC','TCGA-LUAD'),], aes(label=Label,size=4,x=Oncogene),vjust=-0.75, position = position_dodge(0.9))+
        facet_wrap(~Disease,scales='free_y')+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="log2", labels = scientific)+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Oncogene")
    print(p1)
    dev.off()    
    
  collapsed<-cast(melt(complete[,c('Disease','Gene','Structure','Oncogene','fold_change')]),fun=geometric_mean)
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_complete_box.pdf",sep=""),width=12,height=8)
    p1 <- ggplot(casted, aes(x=Oncogene, y=fold_change,fill=Structure))+
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        geom_boxplot(draw_quantiles = c(0.25, 0.5, 0.75))+
        facet_wrap(~Disease)+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="log2", labels = scientific)+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()
    
unique(tpm[tpm$Disease=='TCGA-GBM']$Gene    )
write.csv(unique(tpm[tpm$Disease=='TCGA-GBM',]$Gene,file='gene.csv'))

ecdna<-tpm[,c('Disease','Gene','Structure')]
ecdna$Count<-1
ecdna<-cast(melt(ecdna),fun=sum)
tecdna<-cast(melt(data.frame(ecdna[,c('Disease', 'Gene','Count')])),fun=sum)
tecdna$keeper<-paste(tecdna$Disease,tecdna$Gene)
names(tecdna)<-c('Disease','Gene','Total','keeper')
ecdna$keeper<-paste(ecdna$Disease,ecdna$Gene)
ecdna<-merge(ecdna,tecdna,by=c('keeper','Gene','Disease'))
ecdna$Percent<-ecdna$Count/ecdna$Total
ecdna[ecdna$Structure == 'None' & ecdna$Percent < .25,]$keeper

gbm<-subtpm[subtpm$Gene_name %in% c('EGFR','MYC','ERBB2','CCND1','CCND3','TERT'),]
gbm[gbm$Structure!='None',]$Structure<-'Amplicon'
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_oncogenes.pdf",sep=""),width=12,height=8)
    p1 <- ggplot(gbm, aes(x=Gene_name, y=fold_change,fill=Structure))+
        geom_boxplot(draw_quantiles = c(0.25, 0.5, 0.75))+
        facet_wrap(~Disease,scales='free')+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="log2", labels = scientific)+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()

genes<-c('VSTM2A','SEC61G', 'EGFR','LANCL2','VOPP1','FKBP9L','SEPT14','ZNF713', 'GBAS', 'PSPH', 'SUMF2')
gbm<-subtpm[subtpm$Disease=='TCGA-GBM' & subtpm$Gene_name %in% genes,]
gbm$Gene_name<-factor(gbm$Gene_name, levels=genes)
gbm[gbm$Structure!='None',]$Structure<-'Amplicon'
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_gbm.pdf",sep=""),width=12,height=8)
    p1 <- ggplot(gbm, aes(x=Gene_name, y=fold_change,group=Sample,color=Structure,shape=Structure))+
        #geom_boxplot()+
        geom_line()+
        geom_point()+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+        
        scale_y_continuous(trans="log2", labels = scientific)+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()

pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_gbm_cnv.pdf",sep=""),width=12,height=8)
    p1 <- ggplot(gbm, aes(x=CNV, y=fold_change,color=Structure,shape=Structure))+
        #geom_boxplot()+
        geom_point()+
        ylab("Fold-change")+      
        facet_wrap(~Gene_name,scales='free')+   
        geom_hline(yintercept=1)+        
        scale_y_continuous(trans="log2", labels = scientific)+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("CNV")
    print(p1)
    dev.off()
    

}


sample_analysis<-function() {
  disease = 'ALL'
  data<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/egf_survival.',disease, '.csv',sep=''))
  data$Amp<-'False'
  data[data$Total_Size>10000,]$Amp<-'True'
  data$status<-0
  data[data$Status=='Dead',]$status<-1
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/egfr_curve.",disease,".pdf",sep=""),width=8,height=8)
  p <- ggplot(data, aes(x=Death, y=EGFR, color=Disease, shape=Oncogene, ))+
    geom_point()+
    xlab('Days till death')+
    ylab('EGFR amplification')
  print(p);
  dev.off()

  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/amplicon_size_curve.",disease,".pdf",sep=""),width=5,height=5)
  p <- ggplot(data, aes(x=Death, y=Total_Size+1, color=Disease, shape=Amp))+
    scale_y_log10()+  
    geom_point()+    
    xlab('Days till death')+
    ylab('Amplicon Size')
  print(p);
  dev.off()

  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/breakpoint_curve.",disease,".pdf",sep=""),width=5,height=5)
  p <- ggplot(data, aes(x=Death, y=Breakpoint_Edges, color=Disease, shape=Amplicon))+
    geom_point()+    
    xlab('Days till death')+
    ylab('Breakpoint edges')
  print(p);
  dev.off()
  
  data$SurvObj <- with(data, Surv(Death, status == 2))
  km.by.amp <- survfit(SurvObj ~ Amp, data = data, conf.type = "log-log")
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/amp_curve.",disease,".pdf",sep=""),width=5,height=5)
  p<-autoplot(km.by.amp)
  print(p)
  dev.off()

  km.by.disease <- survfit(SurvObj ~ Disease, data = data, conf.type = "log-log")
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/gbm_curve.pdf",sep=""),width=5,height=5)
  p<-autoplot(km.by.disease)
  print(p)
  dev.off()

  km.by.all <- survfit(SurvObj ~ Disease + Amp, data = data, conf.type = "log-log")
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/gbm_amp_curve.pdf",sep=""),width=5,height=5)
  p<-autoplot(aareg(km.by.all))
  print(p)
  dev.off()

  data$SurvObj <- with(data, Surv(Death, status))  
  form <- formula(data$SurvObj ~ Disease+Amp)   
  cox_bmt <- coxph(form,data = data)
  cox_fit_bmt <- survfit(cox_bmt)  
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/cox_curve.pdf",sep=""),width=5,height=5)  
  p<-plot(cox_fit_bmt)
  print(p)
  dev.off()
  
  data<-data[data$Disease != 'TCGA-LUSC',]
  data<-data[data$Disease != 'TCGA-OV',]
  data[data$Total_Size>10000 & data$Breakpoint_Edges>0,]$Amp<-'True'  
  for (disease in unique(data$Disease)) {
    subdata<-data[data$Disease==disease,]
    fit <- survfit(Surv(Death, status) ~ Amp, data = subdata)
    pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/amp_curve.", disease,".pdf",sep=""),width=5,height=5)
    p<-ggplot2::autoplot(fit)
    print(p)
    dev.off()    
  }

  for (disease in c('TCGA-GBM','TCGA-LGG')) {
    subdata<-data[data$Disease==disease,]
    fit <- survfit(Surv(Death, status) ~ Amp, data = subdata)
    pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/amp_curve.", disease,".pdf",sep=""),width=5,height=5)
    p<-ggplot2::autoplot(fit)
    print(p)
    dev.off()    
  }
  
  lfit <- aareg(Surv(Death, status) ~ Amp + Disease, data=data,
                     nmin=1)

  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/regression_curve.", disease,".pdf",sep=""),width=20,height=20)
  p<-ggplot2::autoplot(aareg(Surv(Death, status) ~ Amp + Disease, data = data))
  print(p)
  dev.off()    
  
  
  disease = 'ALL'
  tpm<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/gene_amp.csv',sep=''))
  #Uncorrected
  temp<-cast(melt(tpm[tpm$Structure=='None',c('Disease', 'FPKM','Gene'),]),fun=mean)
  names(temp)<-c('Disease', 'Gene','AvFPKM')
  merged<-merge(tpm,temp,by=c('Disease','Gene'))
  merged$fold_change<-(1+merged$FPKM)/(1+merged$AvFPKM)
  casted<-cast(melt(merged[,c('Disease','Gene','Structure','Oncogene','fold_change')]),fun=geometric_mean)
  temp<-casted[casted$Structure=='None',]
  casted$Gene<-factor(casted$Gene,levels=unique(temp[order(temp$fold_change),]$Gene))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_gene.pdf",sep=""),width=16,height=4)
    p1 <- ggplot(casted, aes(x=Gene,y=fold_change, shape=Structure, color=Structure))+
        facet_wrap(~Oncogene,scales='free_x')+
        geom_point()+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="asinh")+   
        theme_bw() +             
        theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()
  
  keeper<-tpm[,c('Disease','Gene','FPKM')]
  keeper[keeper$FPKM > 0,]$FPKM<-1
  keeper$Total<-1
  t<-cast(melt(data.frame(keeper)),fun=sum)
  subt<-t[t$FPKM/t$Total> 0.50,]
  subt$keeper<-paste(subt$Disease,subt$Gene)
  tpm$keeper<-paste(tpm$Disease,tpm$Gene)
  subtpm<-tpm[tpm$keeper %in% subt$keeper,]
  temp<-cast(melt(subtpm[subtpm$Structure=='None',c('Disease', 'FPKM','Gene'),]),fun=mean)
  names(temp)<-c('Disease', 'Gene','AvFPKM')
  merged<-merge(subtpm,temp,by=c('Disease','Gene'))
  merged$fold_change<-(1+merged$FPKM)/(1+merged$AvFPKM)
  casted<-cast(melt(merged[,c('Disease','fold_change','Structure','Gene')]),fun=geometric_mean_c)
  keeper<-casted[,c('Disease','Gene')]
  keeper$count<-1
  keeps<-cast(melt(data.frame(keeper)))
  keeps<-keeps[keeps$count==3,]
  keeps$keeps<-paste(keeps$Disease,keeps$Gene)
  casted$keeps<-paste(casted$Disease,casted$Gene)
  oncos<-unique(tpm[,c('Oncogene','Gene')])
  casted<-merge(casted,oncos)
  casted$Structure<-factor(casted$Structure,levels=c('Cycle','Amplicon','None'))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_violin.pdf",sep=""),width=12,height=8)
    p1 <- ggplot(casted[casted$keeps %in% keeps$keeps & casted$fold_change != 0,], aes(x=Oncogene, y=fold_change,fill=Structure))+
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        #geom_boxplot(draw_quantiles = c(0.25, 0.5, 0.75))+
        facet_wrap(~Disease)+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="log2")+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()
  
  subc$FPKM<-subtpm$FPKM/subtpm$CNV
  temp<-cast(melt(subc[subc$Structure=='None',c('Disease', 'FPKM','Gene'),]),fun=mean)
  merged<-merge(subc,temp,by=c('Disease','Gene'))
  merged$fold_change<-(1+merged$FPKM)/(1+merged$AvFPKM)
  casted<-cast(melt(merged[,c('Disease','fold_change','Structure','Gene')]),fun=geometric_mean_c)
  keeper<-casted[,c('Disease','Gene')]
  keeper$count<-1
  keeps<-cast(melt(data.frame(keeper)))
  keeps<-keeps[keeps$count==3,]
  keeps$keeps<-paste(keeps$Disease,keeps$Gene)
  casted$keeps<-paste(casted$Disease,casted$Gene)
  oncos<-unique(tpm[,c('Oncogene','Gene')])
  casted<-merge(casted,oncos)
  casted$Structure<-factor(casted$Structure,levels=c('Cycle','Amplicon','None'))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_wbox.pdf",sep=""),width=12,height=8)
    p1 <- ggplot(casted[casted$keeps %in% keeps$keeps & casted$fold_change != 0,], aes(x=Oncogene, y=fold_change,fill=Structure))+
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
        geom_boxplot(draw_quantiles = c(0.25, 0.5, 0.75))+
        facet_wrap(~Disease)+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="log2")+   
        theme_bw() +             
        theme(panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()
    

sub<-merged[merged$Disease=='TCGA-BLCA' & merged$Gene=='Y_RNA',]
      

  merged$Amplicon<-paste(merged$Structure)
  merged[merged$Structure != 'None',]$Amplicon<-'Amplicon'
  casted<-cast(melt(merged[,c('Gene','Amplicon','Oncogene','fold_change')]),fun=geometric_mean)
  temp<-casted[casted$Amplicon=='None',]
  casted$Gene<-factor(casted$Gene,levels=unique(temp[order(temp$fold_change),]$Gene))
  pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_amplicon_gene.pdf",sep=""),width=16,height=4)
    p1 <- ggplot(casted, aes(x=Gene,y=fold_change, shape=Amplicon, color=Amplicon))+
        #facet_wrap(~Oncogene,scales='free_x')+
        geom_point()+
        ylab("Fold-change")+ 
        geom_hline(yintercept=1)+
        scale_y_continuous(trans="asinh")+   
        theme_bw() +             
        theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
        xlab("Gene")
    print(p1)
    dev.off()    
    
  
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_gene.pdf",sep=""),width=26,height=26)
  p1 <- ggplot(merged, aes(x=Gene,y=fold_change, shape=Structure, color=Structure))+
      facet_wrap(~Oncogene,scales='free_x')+
      geom_jitter(height=0)+
      ylab("Fold-change")+ 
      geom_hline(yintercept=1)+
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Gene")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()
  
  disease = 'ALL'
  tpm<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/cyclic_amplicon_tpm.csv',sep=''))
  subtpm<-tpm[tpm$Size>=50000,]
  temp<-cast(melt(subtpm[subtpm$Recurrent=='False',c('FPKM','Amplicon_id'),]),fun=mean)
  names(temp)<-c('Amplicon_id','av_fpkm')
  merged<-merge(subtpm,temp,by=c('Amplicon_id'))
  merged$fold_change<-(1+merged$FPKM)/(1+merged$av_fpkm)
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_oncogene.pdf",sep=""),width=26,height=26)
  p1 <- ggplot(merged, aes(x=Amplicon_id,y=fold_change, color=Recurrent))+
      facet_wrap(Disease~Oncogene,scales='free_x')+
      geom_point()+
      scale_color_discrete(guide='none')+ 
      ylab("Fold-change")+ 
      geom_hline(yintercept=1)+
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()



  
  
  disease = 'ALL'
  tpm<-read.csv(paste('/pedigree2/projects/namphuon/programs/pancancer/analyses/amplicon_tpm.csv',sep=''))
  tpm<-tpm[tpm$Size>500000 & tpm$Breakpoints > 1,]
  temp<-cast(melt(tpm[tpm$Is_amplicon=='False',c('Reads','Amplicon','Amplification')]),fun=mean)
  names(temp)<-c('Amplicon','AvReads','AvAmp')
  merged<-merge(tpm,temp)
  merged$Classification<-paste(merged$Is_amplicon)
  merged[merged$Is_amplicon=='True',]$Classification<-"Shared overlap"  
  merged[substr(merged$Amplicon,1,12)==merged$Sample,]$Classification<-"Actual"
  merged$fold_change<-(0.01+merged$Reads)/(0.01+merged$AvReads)
  melted<-cast(melt(merged[merged$Is_amplicon=='True',c('Amplicon','fold_change')]),fun=mean)
  merged$Amplicon<-factor(merged$Amplicon,levels=melted[order(melted$fold_change),]$Amplicon)
  
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_oncogene.pdf",sep=""),width=24,height=24)
  p1 <- ggplot(merged[merged$Is_amplicon=='False',], aes(x=Amplicon,y=fold_change,))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),position='dodge')+
      geom_point(position=position_dodge(width=0.9), data=merged[merged$Is_amplicon=='True',],aes(x = Amplicon, y = fold_change,color=Classification))+
      facet_wrap(Disease~Oncogene,scales='free')+
      scale_color_discrete(guide='none')+ 
      ylab("Fold-change")+ 
      geom_hline(yintercept=1)+
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()
    
    
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_recurrent.pdf",sep=""),width=24,height=24)
  p1 <- ggplot(merged[merged$Is_amplicon=='False',], aes(x=Amplicon,y=fold_change,))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),position='dodge')+
      geom_point(position=position_dodge(width=0.9), data=merged[merged$Is_amplicon=='True',],aes(x = Amplicon, y = fold_change,color=Classification))+
      facet_wrap(Disease~Recurrent,scales='free')+
      scale_color_discrete(guide='none')+ 
      ylab("Fold-change")+ 
      geom_hline(yintercept=1)+
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()    
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_recurrent.pdf",sep=""),width=24,height=24)
  p1 <- ggplot(tpm[tpm$Is_amplicon=='False',], aes(x=Amplicon,y=fold_change,))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),position='dodge')+
      geom_point(position=position_dodge(width=0.9), data=tpm[tpm$Is_amplicon=='True',],aes(x = Amplicon, y = fold_change,color='Red'))+
      facet_wrap(Disease~Recurrent,scales='free')+
      scale_color_discrete(guide='none')+ 
      ylab("Fold-change")+ 
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()
            
  #subtpm=tpm[tpm$Disease %in% c('TCGA-GBM','TCGA-BRCA','TCGA-LGG','TCGA-CESC'),]
  subtpm=tpm[tpm$Disease %in% c('TCGA-GBM'),]
  subtpm$type<-paste(subtpm$Disease,subtpm$Oncogene,sep=',')
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_oncogene.gbm.pdf",sep=""),width=10,height=4)
  p1 <- ggplot(subtpm[subtpm$Is_amplicon=='False',], aes(x=Amplicon,y=fold_change))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),position='dodge')+
      geom_point(position=position_dodge(width=0.9), data=subtpm[subtpm$Is_amplicon=='True',],aes(x = Amplicon, y = fold_change,color='Red'))+
      facet_wrap(~type,scales='free_x',ncol=2)+
      scale_color_discrete(guide='none')+ 
      geom_hline(yintercept=1)+
      ylab("Fold-change")+ 
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()    
  
  subtpm$type<-paste(subtpm$Disease,subtpm$Recurrent,sep=',')  
pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_recurrent.gbm.pdf",sep=""),width=10,height=4)
  p1 <- ggplot(subtpm[subtpm$Is_amplicon=='False',], aes(x=Amplicon,y=fold_change))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),position='dodge')+
      geom_point(position=position_dodge(width=0.9), data=subtpm[subtpm$Is_amplicon=='True',],aes(x = Amplicon, y = fold_change,color='Red'))+
      facet_wrap(~type,scales='free_x',ncol=2)+
      scale_color_discrete(guide='none')+ 
      geom_hline(yintercept=1)+
      ylab("Fold-change")+ 
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()    
  
 
 pdf(paste("/pedigree2/projects/namphuon/programs/pancancer/pdf/tpm_disease_recurrent.pdf",sep=""),width=24,height=24)
  p1 <- ggplot(tpm[tpm$Is_amplicon=='False',], aes(x=Amplicon,y=fold_change))+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),position='dodge')+
      geom_point(position=position_dodge(width=0.9), data=tpm[tpm$Is_amplicon=='True',],aes(x = Amplicon, y = fold_change,color='Red'))+
      facet_grid(Disease~Recurrent,scales='free')+
      scale_color_discrete(guide='none')+ 
      ylab("Fold-change")+ 
      scale_y_continuous(trans="log2")+   
      theme_bw() +             
      theme(legend.position=c(0.9,1.), axis.text.x = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), axis.line = element_line(colour = "black"))+                 
      xlab("Amplicon")
  p1$layout$clip[p1$layout$name == "panel"] <- "off"
  print(p1)
  dev.off()  
}

stats<-function() {
  x = 1:39
  hits = 0
  total=1000000
  for (c in 1:total) {
    res = sample(x, 5)
    hits = hits + length(res[res==1])
  }  
  hits/total
}

histogram<-function() {
  data<-read.csv('/pedigree2/projects/namphuon/results/hbv/analyses/tree_analyses/stats/monophyly.r.csv')
 pdf(paste("/pedigree2/projects/namphuon/results/hbv/analyses/pdf/clades.histogram.pdf",sep=""),width=8,height=8)
  p <- ggplot(data, aes(x=Clade, y=Percent*100,fill=Method))+
    geom_bar(width = 0.5,stat='identity',position=position_dodge())+    
    xlab('Clade')+
    scale_fill_discrete(breaks = c('genomes_muscle_best_boot', 'genomes_manual_best_boot', 'genomes_pasta_linearized_best_boot', 'genomes_manual_sregion_best', 'Total_manual','Total_upp_manual','Total_upp_pasta'),
  labels = c('Genomes(Muscle)', 'Genomes(Manual)', 'Genomes(Pasta)', 'Genomes(S-Region)', 'Total(Manual)','Total(UPP-Manual)','Total(UPP-PASTA)'))+  
    ylab('Percent clade accuracy')+
    facet_wrap(~Support,ncol=1)
  print(p);
  dev.off()
    
data$Species<-'Genomes only'
data[grep('Total',data$Method),]$Species<-'Total alignment' 
pdf(paste("/pedigree2/projects/namphuon/results/hbv/analyses/pdf/clades.counts.histogram.pdf",sep=""),width=8,height=8)
  p <- ggplot(data, aes(x=Clade, y=Total/Percent-Total,fill=Method))+
    geom_bar(width = 0.5,stat='identity',position=position_dodge())+    
    xlab('Clade')+
    scale_fill_discrete(breaks = c('genomes_muscle_best_boot', 'genomes_manual_best_boot', 'genomes_pasta_linearized_best_boot', 'genomes_manual_sregion_best', 'Total_manual','Total_upp_manual','Total_upp_pasta'),
  labels = c('Genomes(Muscle)', 'Genomes(Manual)', 'Genomes(Pasta)', 'Genomes(S-Region)', 'Total(Manual)','Total(UPP-Manual)','Total(UPP-PASTA)'))+  
    ylab('Total number of incorrect subtypes ')+
    facet_grid(Species~Support,scales='free_y')
  print(p);
  dev.off()    
}

gbm<-function() {
  paste("Working\n")
  return(0)
  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/turner2017/rnaseq//egf.csv')
  bins<-c()
  bin_size=floor((max(data$Position)-min(data$Position))/10000)
  samples = unique(data$Patient)
  types = unique(data$Type)
  for (bin in (1:10000)) {
    print(bin)
    start=(bin-1)*bin_size+min(data$Position)
    end=start+bin_size
    for (sample in samples) {
      for (type in types) {
        tmp<-data[data$Patient == sample & (data$Position <= end & start <= data$Position) & data$Type==type,]
        av<-sum(tmp$Coverage)/(end-start)
        bins<-rbind(bins, c(sample,type,start,end,av))
      }
    }
  }
  bins<-data.frame(bins)
  names(bins)<-c('Sample','Type','Start','End','Coverage')
  bins$Start<-as.numeric(paste(bins$Start))
  bins$End<-as.numeric(paste(bins$End))  
  bins$Coverage<-as.numeric(paste(bins$Coverage))  
  save(bins,file='/pedigree2/projects/namphuon/tmp/results/turner2017/rnaseq//bins.dat')
pdf(paste("/pedigree2/projects/namphuon/tmp/results/turner2017/rnaseq//pdfs/coverage.pdf",sep=""),width=6,height=6)
  p <- ggplot(bins, aes(x = (Start+End)/2, y = log(Coverage+1), color=Type))+
    geom_segment(x=55086725,xend=55224644,y=0,yend=0,color='black')+
    geom_vline(xintercept = 54817993,, size = 1, colour = "red", linetype = "dashed") +
    geom_vline(xintercept = 55307917,, size = 1, colour = "red", linetype = "dashed") +
    geom_segment(x=54819940,xend=54826939,y=0,yend=0,color='black')+
    annotate('text',x=(55086725+55224644)/2,y=-1,label='EGF')+
    annotate('text',x=(54819940+54826939)/2,y=-1,label='SEC61G')+  
    geom_line()    
  print(p);
  dev.off()
  return(0)

  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/turner2017/rnaseq//segments.csv')
  data$Mid<-(data$End+data$Start)/2
  bins<-c()
  bin_size=floor((max(data$End)-min(data$Start))/1000)
  samples = unique(data$Sample)
  for (bin in (1:1000)) {
    print(bin)
    start=(bin-1)*bin_size+min(data$Start)
    end=start+bin_size
    for (sample in samples) {
      tmp<-data[data$Sample == sample & (data$Start <= end & start <= data$End),]
      av<-sum((tmp$End-tmp$Start)*tmp$Coverage)/sum(tmp$End-tmp$Start)
      bins<-rbind(bins, c(sample,start,end,av))
    }
  }
  bins<-data.frame(bins)
  names(bins)<-c('Sample','Start','End','Coverage')
  bins$Start<-as.numeric(paste(bins$Start))
  bins$End<-as.numeric(paste(bins$End))  
  bins$Coverage<-as.numeric(paste(bins$Coverage))  
  bins[is.na(bins$Coverage),]$Coverage<-0
  #save(bins,file='/pedigree2/projects/namphuon/tmp/results/turner2017/rnaseq/bins.dat')
  subbins<-bins[bins$Start > 54729801 & bins$End < 56176800,]
  subbins<-bins[bins$Start > 55056725 & bins$End < 55254644,]
  
   pdf(paste("/pedigree2/projects/namphuon/tmp/results/turner2017/rnaseq//pdfs/coverage.pdf",sep=""),width=6,height=6)
  p <- ggplot(subbins, aes(x = (Start+End)/2, y = log(Coverage+1), color=Sample))+
    facet_wrap(~Sample,ncol=1)+
    geom_segment(x=55086725,xend=55224644,y=0,yend=0,color='black')+
    annotate('text',x=(55086725+55224644)/2,y=-1,label='EGF')+
    geom_line()    
  print(p);
  dev.off()
  
}

pan_cancer<-function() {
  meta<-read.csv('/pedigree2/projects/namphuon/tmp/data/PANCANCER/metadata/samples.csv')
  counts<-meta[,c('primary_site','rnaseq_size','wgs_size')]
  casted<-cast(melt(counts))  
  casted$Cancer<-paste(casted$primary_site)
  casted$Counts<-casted$rnaseq_size
  casted$Cancer<-factor(casted$Cancer,levels=casted[order(-casted$rnaseq_size),]$Cancer)
  pdf(paste("/pedigree2/projects/namphuon/tmp/data/PANCANCER/pdfs/distribution.pdf",sep=""),width=6,height=6)
  p <- ggplot(casted, aes(x=Cancer, y=Counts))+
    geom_bar(width = 0.5,stat='identity')+    
    xlab('Cancer by tissue type')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))    
  print(p);
  dev.off()
    
  counts<-meta[,c('primary_site','rnaseq_size','wgs_size')]
  casted<-cast(melt(counts),fun=sum)  
  casted$rnaseq_size<-casted$rnaseq_size/1e9
  casted$wgs_size<-casted$wgs_size/1e9

  
}

bret_dove<-function() {
  rates<-read.csv('/pedigree2/projects/namphuon/results/bret_dove/concat/cg.csv')
  subrates<-cast(melt(rates[rates$codon==2,c('species','cg')]),fun=mean)
  rates$species<-factor(rates$species,levels=subrates[order(-subrates$cg),]$species)
  
  subrates<-cast(melt(rates[rates$codon==2,c('gene','cg')]),fun=mean)
  rates$gene<-factor(rates$gene,levels=subrates[order(-subrates$cg),]$gene)
  
  
pdf(paste("/pedigree2/projects/namphuon/results/bret_dove/concat/cg.pdf",sep=""),width=10,height=6)
  p <- ggplot(rates, aes(species, cg))+
    xlab("Species")+
    ylab("% CG")+
    facet_wrap(~codon,ncol=1)+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+    
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=5)    
  print(p);
  dev.off()

av_cg <- mean(subrates$cg)
sd_cg <- sd(subrates$cg)
good_genes = subrates[abs(subrates$cg-av_cg)/sd_cg <= 2,]$gene
write.csv(,file='/pedigree2/projects/namphuon/results/bret_dove/concat/good_genes.csv')

temp_rates<-rates[rates$gene %in% subrates[abs(subrates$cg-av_cg)/sd_cg <= 2,]$gene,]
pdf(paste("/pedigree2/projects/namphuon/results/bret_dove/concat/cg_genes.pdf",sep=""),width=10,height=6)
  p <- ggplot(temp_rates, aes(gene, cg))+
    xlab("Genes")+
    ylab("% CG")+
    facet_wrap(~codon,ncol=1)+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+    
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=5)    
  print(p);
  dev.off()

  gtr<-read.csv("/pedigree2/projects/namphuon/results/bret_dove/concat/gtr.csv",sep=",",header=F);
  names(gtr)<-c('gene','codon','alpha','ac','ag','at','cg','ct','gt')  
  df<-data.frame((gtr[,c('alpha','ac','ag','at','cg','ct')]))
  df$name<-paste(gtr$gene,gtr$codon)
  df$gene<-gtr$gene
  df$codon<-gtr$codon
  
  bad_names<-c()
  for (i in 0:10) {
    for (name in c('alpha','ac','ag','at','cg','ct')) {
      df$outlier <- (df[,name]-mean(df[,name]))/sd(df[,name])
      bad_names<-union(bad_names,df[abs(df$outlier)>10,]$name)
    }
    df<-df[which(!(df$name %in% bad_names)),]
  }
  

  pcaf<-princomp(scale(df[,c('alpha','ac','ag','at','cg','ct')]))
  pca<-data.frame(pcaf$scores[,c(1,2)])    
  pca$codon<-df$codon
  pca$codon<-factor(pca$codon,levels=c(1,2,3))
  
pdf(paste("/pedigree2/projects/namphuon/results/bret_dove/concat/pca.codon.pdf",sep=""),width=6,height=6)
      p<-ggplot(pca,aes(x=Comp.1,y=Comp.2,color=pca$codon))+
        xlab("PCA 1")+
        ylab("PCA 2")+
        scale_color_discrete(name='Codon position')+      
        geom_point(size=1)
      print(p);
      dev.off()
      
      
    
  kres <- c()
  for (i in 2:30) {
    res<-kmeans(pca[,c("Comp.1","Comp.2")],i,nstart=25)
    kres<-rbind(kres,c(i,res$betweenss/res$totss,sum(res$withinss)))
  }
  kres<-data.frame(kres)
  names(kres)<-c('clusters','betweenss_totalss_ratio','within_ss')

  res<-kmeans(pca[,c("Comp.1","Comp.2")],7,nstart=100)
  classes<-data.frame(df$gene,df$codon,res$cluster)
  names(classes)<-c('gene','codon','class')
  write.csv(classes,file='/pedigree2/projects/namphuon/results/bret_dove/concat/class.csv',quote = F,row.names = F)
  
  
}

support<-function() {
  support<-read.csv('/pedigree2/projects/namphuon/results/corona/julie/pasta_3/support.csv',sep='\t')
pdf(paste("/pedigree2/projects/namphuon/results/corona/julie/pasta_3/support.pdf",sep=""),width=6,height=6)
  p <- ggplot(support, aes(Group, Support))+
    xlab("Group")+
    ylab("Support")+
    geom_boxplot()+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=5)    
  print(p);
  dev.off()
  
  
}


asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

intersect_gene<-function(seg_1, seg_2, threshold=0) {
  if (seg_1[1] != seg_2[1]) {
    return(FALSE)
  }
  return(max(as.numeric(seg_1[2]),as.numeric(seg_2[2])) <= min(as.numeric(seg_1[3]),as.numeric(seg_2[3]))+threshold)
}

directionality<-function() {
  all_reads <- read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/statistics/all_samples.csv.tots',sep='\t')
  all_genes <- read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/statistics/all_samples.genes.csv.tots',sep='\t')
  chimeric<-cast(melt(all_reads[,c('Strand','Forward','Backward')]),fun=sum)
  genes<-cast(melt(all_genes[,c('Gene','Gene_direction','Forward','Reverse')]),fun=sum)
  
}

redo_analyses<-function() {
patients<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/patient.activity.5000.50.csv',sep=',')
genes<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/patient.activity.annotation.5000.50.csv',sep=',')

  patients$Name<-paste(patients$Patient,patients$Segment)
  genes$Name<-paste(genes$Patient,genes$Segment)
  
  for (patient in unique(patients$Patient)) {
    vline<-unique(patients[patients$Patient==patient,c('Name','Center')])
    g<-genes[genes$Patient==patient,]
    num=dim(vline)[1]
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/tmp/",patient,".pdf",sep=""),width=6,height=6*num)
  p <- ggplot(patients[patients$Patient == patient,], aes(Position, Coverage))+
            geom_point()+
            geom_line()+
            geom_vline(data=vline, aes(xintercept=Center, linetype="dashed", color = "red"))+
            xlab("Position")+
            ylab("Coverage")+
            geom_segment(data=g,aes(x=g$Start,xend=g$End,y=-g$Position,yend=-g$Position))+
            geom_text(data=g,aes(x = g$Start, y = g$Position, label = g$Type))+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            facet_wrap(~Name,ncol=1,scales='free')          
          print(p);
          dev.off()    
  }
  
    for (patient in unique(patients$Patient)) {
      tmp <-patients[patients$Patient==patient,]
      for (name in unique(tmp$Name)) {
        temp<-tmp[tmp$Name==name,]
        vline<-unique(tmp[tmp$Name==name,c('Name','Center')])
        g<-genes[genes$Name==name,]
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/tmp/",name,".pdf",sep=""),width=6,height=6)
  p1 <- ggplot(temp, aes(Position, Coverage))+
            geom_point()+
            geom_line()+
            ggtitle(name)+
            geom_vline(data=vline, aes(xintercept=Center, linetype="dashed", color = "red"))+
            xlab("Position")+
            ylab("Coverage")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")
  if (dim(g)[1] != 0) {
    g$Position<-1:dim(g)[1]
    p2 <- ggplot(g, aes(Start, Position))+
              geom_segment(data=g,aes(x=g$Start,xend=g$End,y=g$Position,yend=g$Position))+
              geom_text(data=g,aes(x = g$Start, y = g$Position-0.25, label = g$Type))+
              xlim(min(temp$Position),max(temp$Position))+
              theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")
            grid.newpage()
            grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))                         
            print(grid)
  } else {
    print(p1)
  }
          dev.off()            
      }
    }

#all_expression<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_rna_expression_final.csv.bak.nostrand',sep=',')
all_expression<-read.csv('/pedigree2/projects/namphuon/tmp//results/TCGA-CESC/analyses/wgs_rna_expression_final.csv.bak.1',sep=',')

wilcox.test(all_expression$mean_baseline_FPKM,all_expression$Total_FPKM, paired=TRUE)
wilcox.test(all_expression[all_expression$Fusion_mRNA_Present=='True',]$mean_baseline_FPKM,all_expression[all_expression$Fusion_mRNA_Present=='True',]$Total_FPKM, paired=TRUE)
wilcox.test(all_expression[all_expression$Fusion_mRNA_Present=='False',]$mean_baseline_FPKM,all_expression[all_expression$Fusion_mRNA_Present=='False',]$Total_FPKM, paired=TRUE)

wilcox.test(all_expression$LINE_FPKM,all_expression$LINE_mean_baseline_FPKM, paired=TRUE)
wilcox.test(all_expression$LTR_FPKM,all_expression$LTR_mean_baseline_FPKM, paired=TRUE)


melted<-melt(all_expression[,c('Fusion_mRNA_Present', 'fold_change','LTR_fold_change','LINE_fold_change','Gene_fold_change','Oncogene_fold_change')])
all<-melt(all_expression[,c('fold_change','LTR_fold_change','LINE_fold_change','Gene_fold_change','Oncogene_fold_change')])
all_casted<-melt(cast(all[!is.nan(all$value),],fun=mean))
all_casted_sd<-melt(cast(all[!is.nan(all$value),],fun=sd))
all_casted$sd<-all_casted_sd$value.1
all_casted<-all_casted[,c('value','value.1','variable','sd')]
names(all_casted)<-c('Fusion_mRNA_Present','value','variable','sd')
all_casted$Fusion_mRNA_Present<-'All'

casted<-melt(cast(melted[!is.nan(melted$value),],fun=mean))
casted_sd<-melt(cast(melted[!is.nan(melted$value),],fun=sd))
casted$sd<-casted_sd$value
casted$Fusion_mRNA_Present<-factor(casted$Fusion_mRNA_Present,levels=c('Average','True','False'))
casted<-rbind(casted,all_casted)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expression.summary.barplot.fold.pdf",sep=""),width=6,height=6)
  p <- ggplot(casted, aes(x=variable, y=value, fill=Fusion_mRNA_Present))+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), colour="black", width=.1,position=position_dodge(0.9)) +
    geom_bar(position =position_dodge(), stat='identity')+
    xlab("Type")+
    scale_fill_discrete(name='Genomic segment',breaks=c('False','True','All'),labels=c('No fusion mRNA','Fusion mRNA','All'))+    
    scale_x_discrete('Type', breaks=c('fold_change','LINE_fold_change','LTR_fold_change','Gene_fold_change','Oncogene_fold_change'), labels=c('All transcripts','LINE','LTR','Gene','Oncogene'))+
   theme(axis.text.x = element_text(angle = 45, hjust = 1))+   
    ylab("Log2 fold change")
  print(p);
  dev.off()

melted<-melt(all_expression)
melted<-melted[melted$variable %in% c('fold_change','Oncogene_fold_change','Gene_fold_change','LINE_fold_change','LTR_fold_change'),]
melted$Fusion_mRNA_Present<-paste(melted$Fusion_mRNA_Present)
melted$Types<-melted$Fusion_mRNA_Present
all<-melted
all$Types<-'All'
all<-rbind(melted,all)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expression.summary.violin.fold.pdf",sep=""),width=6,height=6)
  p <- ggplot(all, aes(x=variable, y=value, fill=Types))+
    geom_violin()+
    xlab("Type")+
    scale_fill_discrete(name='Genomic segment',breaks=c('False','True','All'),labels=c('No fusion mRNA','Fusion mRNA','All'))+    
    scale_x_discrete('Type', breaks=c('fold_change','LINE_fold_change','LTR_fold_change','Gene_fold_change','Oncogene_fold_change'), labels=c('All transcripts','LINE','LTR','Gene','Oncogene'))+
   theme(axis.text.x = element_text(angle = 45, hjust = 1))+   
    ylab("Log2 fold change")
  print(p);
  dev.off()



#updown<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/up_down_expression.csv.bak',sep='\t')
  updown<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/up_down_expression.10000.csv.bak',sep='\t')
  updown$Direction<-paste(updown$Direction)
  updown[updown$Direction=='Drive upstream',]$Position<- -1*updown[updown$Direction=='Drive upstream',]$Position
  updown[updown$Direction=='Drive upstream' | updown$Direction=='Drive downstream',]$Direction <- 'Simple integration'
  updown[updown$Direction=='Insufficient support',]$Direction <- 'Complex integration'  
  updown[updown$Direction=='Multiple Integrations',]$Direction <- 'Complex integration'
  updown[updown$Direction=='No fusion mRNA',]$Direction <- 'Fusionless mRNA integration'

  updown[updown$Direction=='Complex integration',]$Direction <- 'Complex integration (51)'
  updown[updown$Direction=='Fusionless mRNA integration',]$Direction <- 'Fusionless mRNA integration (107)'
  updown[updown$Direction=='Simple integration',]$Direction <- 'Simple integration (68)'
  
  updown$fold_change = log((updown$FPKM+1)/(updown$AvFPKM+1))/log(2)

  updown$fold_change = (updown$FPKM+1)/(updown$AvFPKM+1)
  
  updown$Position<-factor(updown$Position,levels=unique(updown$Position))
  updown[updown$Fusion>0,]$Fusion<-1
  updown$Fusion<-factor(updown$Fusion,levels=unique(updown$Fusion))
  summary(melt(updown[updown$Position==0,]$Direction))
  #updown$Direction<-factor(updown$Direction,levels=c('t-,V+','t+,V+', 'V+,t-', 'V+,t+','Unknown','Other', 'No fusion mRNA'))
  #updown$Direction<-factor(updown$Direction,levels=c('Unknown', 'Drive upstream', 'Drive downstream', 'Other', 'No fusion mRNA'))
  #updown[updown$Fusion>1,]$Direction<-'Other'
  #updown[updown$Direction %in% c('t-,V+','t+,V+','Unknown'),]$Direction<-'Other'
  #updown[updown$Fusion==0,]$Direction<-'No fusion mRNA'

wilcox.test(updown[updown$Direction !=  'Fusionless mRNA integration (107)' & updown$Position == 0,]$FPKM,updown[updown$Direction !=  'Fusionless mRNA integration (107)' & updown$Position == 0,]$AvFPKM,paired=TRUE)

  casted<-cast(melt(updown[,c('Direction','Position','FPKM','AvFPKM','fold_change')]),fun=mean)
  melted<-melt(casted)
  melted$Position<-as.numeric(paste(melted$Position))

//updown[updown$Position==-99600 & updown$FPKM != 0,]
// melted$Type<-'Virus+ on human positive strand'
// melted[melted$Direction=='Drive upstream',]$Type<-'Virus+ on human negative strand'
// melted[melted$Direction=='Unknown',]$Type<-'Unable to determine'
// melted[melted$Direction=='No fusion mRNA',]$Type<-'No fusion mRNA'
// melted$Type<-'Virus+,Human+ or Human+,Virus+'
// melted[melted$Direction=='negative',]$Type<-'Virus+,Human- or Human-,Virus+'
// melted[melted$Direction=='unknown' & melted$Fusion == 1,]$Type<-'Unable to determine'
// melted[melted$Direction=='unknown' & melted$Fusion == 0,]$Type<-'No fusion mRNA'
#melted$Type<-factor(melted$Type,levels=c('Virus+ on human positive strand','Virus+ on human negative strand', 'Unable to determine', 'No fusion mRNA'))
#melted$Direction<-factor(melted$Direction,levels=c('Simple integration','Complex integration', 'Fusionless mRNA integration'))
melted$Direction<-factor(melted$Direction,levels=c('Simple integration (68)','Complex integration (51)', 'Fusionless mRNA integration (107)'))
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/updown.10000.pdf",sep=""),width=8,height=3)
require(scales)
p <- ggplot(melted[melted$variable=='fold_change',], aes(Position, 2^value))+
          geom_point(size=0.5)+
          xlab("Position (relative to integration point)")+
          ylab("FPKM ratio")+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          facet_wrap(~Direction,ncol=3)
          #scale_x_continuous(labels = comma)
        print(p);
        dev.off()



wgs_expression<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/all_wgs_type_expression.10000.csv',sep='\t')

#Drop zero values
tmp<-cast(wgs_expression[,c('Segment','Type','FPKM')],Segment~Type,fun=sum)
melted<-melt(tmp[(tmp$Gene+tmp$LINE+tmp$LTR+tmp$Oncogene) != 0,])
melted$name<-paste(melted$Segment,melted$Type)
melted<-melted[melted$value != 0,]
wgs_expression$name<-paste(wgs_expression$Segment,wgs_expression$Type)

####
sub<-wgs_expression[wgs_expression$name %in% melted$name,]
wgs_expression<-sub
ltr <-merged[merged$Type=='LTR',c('Segment','Patient','FPKM','fold','mean_FPKM')]
ltr<-ltr[order(-ltr$fold),]
####

real_analysis<-wgs_expression[wgs_expression$Rep=='Real',]
with_integration<-real_analysis[real_analysis$Integration=='True',]
without_integration<-real_analysis[real_analysis$Integration=='False',]
with_integration$fold<-0
casted<-cast(melt(without_integration[,c('name','FPKM')]),fun=mean)
names(casted)<-c('name','mean_FPKM')
merged<-merge(with_integration,casted)
merged<-merged[merged$name %in% melted$name,]

merged$fold<-log((merged$FPKM+1)/(merged$mean_FPKM+1))/log(2)
cast(melt(merged[,c('Type','fold')]),fun=mean)

temp_1<-cast(melt(wgs_expression[!(wgs_expression$Fusion=='True' & wgs_expression$Integration=='False'),c('Segment','Type','Rep','FPKM','Integration','NumGenes')]),fun=mean)
sim_mean<-cast(melt(data.frame(temp_1[temp_1$Rep !='Real',c('Type','FPKM')])),fun=mean)
sim_sd<-cast(melt(data.frame(temp_1[temp_1$Rep !='Real',c('Type','FPKM')])),fun=sd)

real_mean<-cast(melt(data.frame(temp_1[temp_1$Rep=='Real',c('Type','Integration','FPKM','NumGenes')])),fun=mean)
real_sd<-cast(melt(data.frame(temp_1[temp_1$Rep=='Real',c('Type','Integration','FPKM','NumGenes')])),fun=sd)

real_mean$Segment<-'With integration'
real_mean[real_mean$Integration=='False',]$Segment<-'Without integration'
real_mean$sd<-real_sd$FPKM
sim_mean$Segment<-'Random'
sim_mean$sd<-sim_sd$FPKM
combined<-rbind(real_mean[,c('Type','FPKM','Segment','sd')],sim_mean)
combined$Segment<-factor(combined$Segment,levels=c('With integration', 'Without integration', 'Random'))
combined$Type<-factor(combined$Type,levels=c('LINE','LTR','Gene','Oncogene'))

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expression.summary.barplot.pdf",sep=""),width=6,height=6)
  p <- ggplot(combined, aes(x=Type, y=FPKM, fill=Segment))+
    geom_errorbar(aes(ymin=max(FPKM-sd,0), ymax=FPKM+sd), colour="black", width=.1,position=position_dodge(0.9)) +
    geom_bar(position =position_dodge(), stat='identity')+
    xlab("Type")+
    ylab("Mean FPKM")
  print(p);
  dev.off()

casted<-cast(melt(data.frame(merged[,c('Type','fold','Fusion')])),fun=mean)
casted_sd<-cast(melt(data.frame(merged[,c('Type','fold','Fusion')])),fun=sd)
casted$sd<-casted_sd$fold
casted$Fusion<-factor(casted$Fusion,levels=c('True','False'))
casted<-rbind(casted,c('Oncogene','False',NA,NA))
casted$fold<-as.numeric(casted$fold)
casted$sd<-as.numeric(casted$sd)
#casted[is.na(casted$sd),]$sd<-0
casted[is.na(casted$fold),]$fold<-0

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expression.summary.barplot.fold.pdf",sep=""),width=6,height=6)
  p <- ggplot(casted, aes(x=Type, y=fold, fill=Fusion))+
    geom_errorbar(aes(ymin=fold-0.01, ymax=fold+sd), colour="black", width=.1,position=position_dodge(0.9)) +
    scale_fill_discrete(name='Fusion mRNA present',breaks=c('True','False'),labels=c('Yes','No'))+
    geom_bar(position =position_dodge(), stat='identity')+
    xlab("Type")+
    ylab("log2 Fold Change in FPKM")
  print(p);
  dev.off()



real_fold_mean<-cast(melt(data.frame(temp_1[temp_1$Rep=='Real' & temp_1$Integration=='True',c('Type','Segment','Integration','FPKM','NumGenes')])),fun=mean)
tmp<-cast(melt(data.frame(temp_1[temp_1$Rep=='Real' & temp_1$Integration=='False',c('Type','Segment','FPKM')])),fun=mean)
tmp1<-merge(real_fold_mean,tmp,by=c('Type','Segment'))


tmp1$fold<-log((tmp1$FPKM.x+1)/(tmp1$FPKM.y+1))/log(2)
fw <- fitdist(tmp1$fold, "weibull")


tests<-wgs_expression[wgs_expression$Rep=='Real',]
integrations<-tests[tests$Integration=='True',]
casted<-cast(melt(tests[tests$Integration=='False',c('Segment','Type','FPKM')]),fun=mean)
merged<-merge(integrations,casted,by=c('Segment','Type'))
for (type in unique(merged$Type)) {
  print(type)
  print(t.test(merged[merged$Type==type,]$FPKM.x,merged[merged$Type==type,]$FPKM.y,paired=TRUE))
  #print(ks.test(merged[merged$Type==type,]$FPKM.x,merged[merged$Type==type,]$FPKM.y))
  #print(ks.test(log(merged[merged$Type==type,]$FPKM.x+1)/(merged[merged$Type==type,]$FPKM.y+1)/log(2),'pnorm')

#fn <- fitdist(log(merged[merged$Type==type,]$FPKM.x+1)/(merged[merged$Type==type,]$FPKM.y+1)/log(2), "norm")
#pvalues = dnorm(log(merged[merged$Type==type,]$FPKM.x+1)/(merged[merged$Type==type,]$FPKM.y+1)/log(2), mean=fn$estimate['mean'],sd=fn$estimate['sd'])  
  
}



pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expression.summary.barplot.pdf",sep=""),width=6,height=6)
  p <- ggplot(combined, aes(x=Type, y=FPKM, fill=Segment))+
    geom_errorbar(aes(ymin=max(FPKM-sd,0), ymax=FPKM+sd), colour="black", width=.1,position=position_dodge(0.9)) +
    geom_bar(position =position_dodge(), stat='identity')+
    xlab("Type")+
    ylab("Mean FPKM")
  print(p);
  dev.off()

copy_count<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/copy_count.csv',sep=',',header=FALSE)
names(copy_count)<-c('Patient','Copy')
merged$fold<-log((merged$FPKM.x+1)/(merged$FPKM.y+1))/log(2)
merged<-merge(merged,copy_count)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expression.copy.pdf",sep=""),width=6,height=6)
  p <- ggplot(merged[merged$Type=='LINE',], aes(x=Copy, y=fold),shape=Type)+
    geom_point()+
    xlab("Copy Count")+
    ylab("Mean FPKM")
  print(p);
  dev.off()


rna_meta<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/fusion_mrna.csv',sep='\t')
rna_meta[rna_meta$orientation=='H->V',]
#Collapse isoforms
rna_meta$isoforms<-paste(rna_meta$patient,substr(rna_meta$read,1,nchar(paste(rna_meta$read))-3))
rna_meta$read<-NULL

t<-c()
for (isoform in unique(rna_meta$isoforms)) {
  temp<-rna_meta[rna_meta$isoforms == isoform,]
  temp<-temp[order(-temp$bitscore),]
  t<-rbind(t,temp[1,])
}
t$classification<-'No significant hit'
t[which(t$evalue < 1e-2),]$classification<-'Significant protein hit'
t[intersect(which(t$classification == 'Significant protein hit'),grep('predict', t$protein,ignore.case=TRUE)),]$classification<-'Predicted protein hit'
t[intersect(which(t$classification == 'Significant protein hit'),grep('hypoth', t$protein,ignore.case=TRUE)),]$classification<-'Hypothetical protein hit'

rna_direction<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/rna_direction.csv')
rna_direction$total<-rna_direction$vf_hf+rna_direction$vf_hr+rna_direction$hf_vf+rna_direction$hr_vf
rna_direction<-rna_direction[rna_direction$total!=0,]
melted<-melt(rna_direction[,c('patient','rna_integration','vf_hr','vf_hf','hf_vf','hr_vf','gene')])
melted$name<-paste(melted$patient,melted$rna_integration)
rna_direction$name<-paste(rna_direction$patient,rna_direction$rna_integration)
rownames(rna_direction)<-rna_direction$name
melted$value<-melted$value/rna_direction[melted$name,]$total


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/rna_direction.boxplot.pdf",sep=""),width=6,height=6)
  p <- ggplot(melted, aes(variable, value))+
    xlab("Direction")+
    ylab("Proportion of reads supporting")+
    geom_boxplot()+
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=5)    
  print(p);
  dev.off()

best<-c()
for (name in unique(melted$name)) {
  t<-melted[melted$name == name & melted$variable != 'gene',]
  best<-rbind(best,t[order(-t$value),][1,])  
}
rna_direction[rna_direction$gene > 0,]$gene = 1
merged<-merge(best,rna_direction[,c('name','gene')])
merged[merged$variable == 'vf_hr',]$variable<-'vf_hf'
merged[merged$variable == 'hr_vf',]$variable<-'hf_vf'

merged$variable<-factor(merged$variable,levels=c('vf_hf','vf_hr','hr_vf','hf_vf','gene'))
name_list <- list(
  #'vf_hf'="Viral+,Human+",
  'vf_hf'="Human downstream of virus",
  'vf_hr'="Viral+,Human-",
  'hr_vf'="Human-,Viral+",
  #'hf_vf'="Human+,Viral+"
  'hf_vf'="Human upstream of virus"
)

type_labeller <- function(variable,value){
  return(name_list[value])
}


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/rna_direction.summary.histogram.pdf",sep=""),width=5,height=3.5)
  p <- ggplot(merged, aes(x=value,fill=factor(gene)))+
    geom_histogram(bins=10)+
    scale_fill_discrete(name='Gene within\n 10kb',breaks=c(0,1),labels=c('No','Yes'))+
    xlab('Proportion of reads supporting direction')+
    ylab("Number of fusion mRNA")+
    facet_wrap(~variable,labeller=type_labeller)    
  print(p);
  dev.off()


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/rna_direction.summary.barplot.pdf",sep=""),width=6,height=6)
  p <- ggplot(melted, aes(x=variable, y=value,fill=factor(variable)))+
    geom_bar(stat='identity')+
    xlab("Direction")+
    ylab("Number of fusion mRNA")+
    scale_fill_discrete(name='Orientation', breaks=c('vf_hr','vf_hf','hf_vf','hr_vf'), label=c('Virus->-Human', 'Virus->+Human', '+Human->Virus', '-Human->Virus')) 
  print(p);
  dev.off()
  
wgs_fpkm_adv<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_expression.csv',sep='\t')
integrations<-cast(melt(wgs_fpkm_adv[wgs_fpkm_adv$Integration=='True',c('Segment','FPKM','Integration','Fusion')]),fun=mean)
no_integrations<-cast(melt(wgs_fpkm_adv[wgs_fpkm_adv$Integration=='False',c('Segment','FPKM')]),fun=mean)
names(no_integrations)<-c('Segment','mean_FPKM')
merged<-merge(integrations,no_integrations)
merged$fold<-log((merged$FPKM+1)/(merged$mean_FPKM+1))/log(2)
merged<-merged[order(-merged$fold),]

groups<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/groups.csv',sep='\t')
casted<-cast(melt(groups[,c('Patient','HasIntegration','Type','FPKM')]),fun=sum)
cast(melt(data.frame(casted[,c('HasIntegration','Type','FPKM')])),fun=sd)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/groups.pdf",sep=""),width=6,height=6)  
p1 <- ggplot(casted, aes(x=HasIntegration,y=FPKM))+
    geom_violin()+
    ylab("Sum FPKM")+
    facet_wrap(~Type,scales='free')+
    xlab("Integration")
print(p1)
dev.off()
  
wgs_fpkm_adv<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_rna_expression.csv')

tmp<-cast(melt(wgs_fpkm_adv[wgs_fpkm_adv$Integration=='True',c('Segment','FPKM')]))  
wgs_fpkm_adv$Segment<-factor(wgs_fpkm_adv$Segment,levels=tmp[order(tmp$FPKM),]$Segment)

wgs_fpkm_adv$delta<-wgs_fpkm_adv$FPKM
wgs_fpkm_adv$log<-wgs_fpkm_adv$FPKM

integrations<-wgs_fpkm_adv[wgs_fpkm_adv$Integration=='True',]
no_integrations<-wgs_fpkm_adv[wgs_fpkm_adv$Integration!='True',]
no_integrations[no_integrations$Genes > 0,]$Genes<-1
no_integrations[no_integrations$Downstream=='No URR',]$Downstream<-'URR'
casted<-cast(melt(no_integrations[,c('Segment','Integration','Downstream','FPKM','Genes')]),fun=mean)

merged<-merge(casted,integrations[,c('Segment','FPKM')])

for (r in 1:dim(integrations)[1]) {
  row = integrations[r,]
  no_integrations[no_integrations$Segment==row$Segment,]$delta<- row$FPKM-no_integrations[no_integrations$Segment==row$Segment,]$delta  
  no_integrations[no_integrations$Segment==row$Segment,]$Total<-   no_integrations[no_integrations$Segment==row$Segment,]$Total<-row$Total    
  
  no_integrations[no_integrations$Segment==row$Segment,]$log<-log( (row$FPKM+1)/(no_integrations[no_integrations$Segment==row$Segment,]$FPKM+1))
}

casted$log<-0
for (r in 1:dim(integrations)[1]) {
  row = integrations[r,]  
  if (dim(casted[casted$Segment==row$Segment,])[1] != 1) {
    break
  }  
  casted[casted$Segment==row$Segment,]$log<-log2( (row$FPKM+1)/(casted[casted$Segment==row$Segment,]$FPKM+1))
}
casted$Genes<-factor(casted$Genes,levels=c(0,1))

name_list <- list(
  'No Fusion mRNA'="Genomic integration with absence of fusion mRNA",
  'No URR'="Human mRNA to Virus mRNA",
  'URR'="Genomic integration with presence of fusion mRNA"
)

type_labeller <- function(variable,value){
  return(name_list[value])
}



pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.adv.histogram.10000.pdf",sep=""),width=6,height=6)
  p <- ggplot(casted, aes(log,fill=Genes))+
    xlab("log_2(Ratio of mean(FPKM))")+
    ylab("Counts")+
    geom_histogram(bins=10)+
    geom_vline(xintercept=2.321,linetype='dashed',colour='red')+
    scale_fill_discrete(name='Genes present',breaks=c(1,0),labels=c('Gene','No Gene'))+
    #scale_y_continuous(trans=asinh_trans(),breaks=c(0,2,4,8, 16, 32, 64, 128))+        
    facet_wrap(~Downstream,ncol=1,labeller=type_labeller)
    #facet_wrap(Scale~Downstream)+
    #theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.x=element_blank())
    #theme(axis.text.x=element_blank())
  print(p);
  dev.off()  


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.adv.boxplot.10000.pdf",sep=""),width=8,height=4)
  p <- ggplot(no_integrations, aes(Segment, delta))+
    xlab("Segment")+
    ylab("delta FPKM")+
    scale_y_continuous(trans=asinh_trans(),breaks=c(-1000, -100,-10,-1,0,1,10,100, 1000))+    
    geom_boxplot()+
    facet_wrap(~Downstream,ncol=3,scales='free_x',labeller=type_labeller)+
    #facet_wrap(Scale~Downstream,scales='free')+
    #theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.x=element_blank())
    theme(axis.text.x=element_blank())
  print(p);
  dev.off()  

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.all.boxplot.10000.pdf",sep=""),width=8,height=4)
  p <- ggplot(no_integrations, aes(Segment, log))+
    xlab("Segment")+
    ylab("log_2 (FPKM ratio)")+
    #scale_y_continuous(trans=asinh_trans(),breaks=c(-1000, -100,-10,-2,0,2,4,8, 16,32,64,128,256,512,1024,2048))+    
    geom_boxplot()+
    #facet_wrap(~Downstream,ncol=1,scales='free_x',labeller=type_labeller)+
    #facet_wrap(Scale~Downstream,scales='free')+
    #theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.x=element_blank())
    theme(axis.text.x=element_blank())
  print(p);
  dev.off()  

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.fusion.boxplot.10000.pdf",sep=""),width=5,height=4)
  p <- ggplot(no_integrations, aes(Segment, log))+
    xlab("Segment")+
    ylab("log_2 (FPKM ratio)")+
    #scale_y_continuous(trans=asinh_trans(),breaks=c(-1000, -100,-10,-2,0,2,4,8, 16,32,64,128,256,512,1024,2048))+    
    geom_boxplot(outlier.size = 0.1)+
    facet_wrap(~Downstream,ncol=1,scales='free_x',labeller=type_labeller)+
    #facet_wrap(Scale~Downstream,scales='free')+
    #theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.x=element_blank())
    theme(axis.text.x=element_blank())
  print(p);
  dev.off()  


type_expression<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_rna_type_expression.csv')
type_expression$delta<-type_expression$FPKM
type_expression$name<-paste(type_expression$Segment,type_expression$Type)
no_integrations<-type_expression[type_expression$Integration=='False',]
for (name in unique(type_expression$name)) {
  delta <-mean(type_expression[type_expression$name==name & type_expression$Integration=='True',]$FPKM)
  no_integrations[no_integrations$name==name,]$delta<-delta-no_integrations[no_integrations$name==name,]$FPKM
}

tmp<-cast(melt(type_expression[type_expression$Integration=='True',c('Segment','delta')]),fun=mean)
no_integrations$Segment<-factor(no_integrations$Segment,levels=tmp[order(tmp$delta),]$Segment)
no_integrations$Type<-factor(no_integrations$Type,levels=c('LINE','LTR','Gene','Oncogene'))
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.types.pdf",sep=""),width=8,height=4)
  p <- ggplot(no_integrations, aes(Segment, delta))+
    xlab("Segment")+
    ylab("delta FPKM")+
    scale_y_continuous(trans=asinh_trans(),breaks=c(-10000, -1000, -100,-10,-1,0,1,10,100, 1000, 10000))+    
    geom_boxplot()+
    #facet_wrap(~Downstream,ncol=1,scales='free_x',labeller=type_labeller)+
    facet_wrap(~Type,scales='free')+
    #theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.text.x=element_blank())
    theme(axis.text.x=element_blank())
  print(p);
  dev.off()  



all_counts<-read.csv('/pedigree2/projects/namphuon/tmp//results/TCGA-CESC/analyses/counts_type.all1.csv')
#all_counts[all_counts$Type=='LTR',]$Type<-'LINE'
casted<-cast(melt(all_counts[,c('Type','Count','Replicate'),]),fun=sum)
casted$Experiment<-'Expected'
casted[casted$Replicate=='Real',]$Experiment<-'Observed'
sds<-cast(melt(data.frame(casted[,c('Type','Count','Experiment')])),fun=sd)
#sds[is.na(sds$Count),]$Count<-0
means<-cast(melt(data.frame(casted[,c('Type','Count','Experiment')])),fun=mean)
names(sds)<-c('Type','Experiment','SD')
merged<-merge(means,sds)
merged$Experiment<-factor(merged$Experiment,levels=c('Observed', 'Expected'))
merged$Type<-factor(merged$Type, levels=c('LINE','LTR','Gene','Oncogene'))
label.df <- data.frame(Group=c('LINE','LTR','Gene','Oncogene'),Value=c(1850,900,225,6))

# Error bars represent standard error of the mean
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expected.counts.pdf",sep=""),width=8,height=4)
p<-ggplot(merged, aes(x=Type, y=Count, fill=Experiment)) + 
    ylab('Number of Integrations')+
    geom_bar(position=position_dodge(), stat="identity") +    
    geom_errorbar(aes(ymin=Count-SD, ymax=Count+SD),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
print(p)
dev.off()

casted<-cast(melt(all_counts[,c('Type','Count','Replicate'),]),fun=sum)
casted$Experiment<-'Expected'
casted[casted$Replicate=='Real',]$Experiment<-'Observed'
sim<-casted[casted$Experiment=='Expected',]
real<-casted[casted$Experiment=='Observed',]
sim$Type<-factor(sim$Type,levels=c('LINE','LTR','Gene','Oncogene'))
#sim[sim$Count==0,]$Count<-1

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expected.counts.violin.pdf",sep=""),width=8,height=4)  
p1 <- ggplot(sim, aes(x=Type,y=Count))+
    facet_wrap(~Type,scales='free')+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
    geom_point(data=real,aes(x = Type, y = Count,color='red'))+
    ylab("Expected count")+    
    theme(legend.position="none")+        
    #scale_y_log10()+        
    xlab("Type")
print(p1)
dev.off()


for (type in unique(casted$Type)) {
  fn <- fitdist(casted[casted$Experiment=='Expected' & casted$Type == type,]$Count, "norm")
  pvalues = pnorm(casted[casted$Experiment=='Observed' & casted$Type == type,]$Count, mean=fn$estimate['mean'],sd=fn$estimate['sd'])  
  print(paste(type,pvalues))  
}

wgs_rna_distances<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.rna.distances.csv')
wgs_rna_meta<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.rna.meta.csv')

  unassigned<-wgs_rna_distances[is.na(wgs_rna_distances$distance),]
  assigned<-wgs_rna_distances[!is.na(wgs_rna_distances$distance),]
  

b = 10^seq(0,9)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.rna.histogram.distances.pdf",sep=""),width=10,height=4)  
  p <- ggplot(assigned, aes(x=distance+1))+
    geom_histogram(bins=30)+
    xlab('kb from nearest WGS integration')+
    ylab('Counts')+
    #scale_x_log10()
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))
  print(p)
  dev.off()
  
  wgs_fpkm<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.10000.fpkm.csv')
  tmp<-wgs_fpkm[wgs_fpkm$integration=='True',c('segment','fpkm')]
  tmp<-cast(melt(tmp),fun=mean)
  wgs_fpkm$segment<-factor(wgs_fpkm$segment, levels=tmp[order(tmp$fpkm),]$segment)
  #wgs_fpkm[wgs_fpkm$integration=='True',]
   pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.50000.boxplot.pdf",sep=""),width=40,height=6)
  p <- ggplot(wgs_fpkm, aes(segment, fpkm+1, fill=integration))+
    xlab("Integration")+
    ylab("FPKM")+
    scale_y_log10()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p);
  dev.off()
  
  wgs_fpkm_adv<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.directions.advance.csv')
wgs_fpkm_adv<-wgs_fpkm_adv[wgs_fpkm_adv$Distance<9000,]
tmp<-cast(melt(wgs_fpkm_adv[wgs_fpkm_adv$Integration=='True',c('Segment','FPKM')]))  
wgs_fpkm_adv$Segment<-factor(wgs_fpkm_adv$Segment,levels=tmp[order(tmp$FPKM),]$Segment)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.adv.boxplot.pdf",sep=""),width=30,height=6)
  p <- ggplot(wgs_fpkm_adv, aes(Segment, FPKM+1, fill=Integration))+
    xlab("Segment")+
    ylab("FPKM")+
    scale_y_log10()+
    geom_boxplot()+
    facet_wrap(~Downstream,ncol=1,scales='free_x')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p);
  dev.off()

  wgs_fpkm_adv<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.directions.advance.csv')
wgs_fpkm_adv<-wgs_fpkm_adv[wgs_fpkm_adv$Distance<9000,]
tmp<-cast(melt(wgs_fpkm_adv[wgs_fpkm_adv$Integration=='True',c('Segment','FPKM')]))  
wgs_fpkm_adv$Segment<-factor(wgs_fpkm_adv$Segment,levels=tmp[order(tmp$FPKM),]$Segment)


wgs_fpkm_adv$delta<-wgs_fpkm_adv$FPKM
integrations<-wgs_fpkm_adv[wgs_fpkm_adv$Integration=='True',]
no_integrations<-wgs_fpkm_adv[wgs_fpkm_adv$Integration!='True',]

for (r in 1:dim(integrations)[1]) {
  row = integrations[r,]
  no_integrations[no_integrations$Segment==row$Segment,]$delta<- row$FPKM-no_integrations[no_integrations$Segment==row$Segment,]$delta  
  no_integrations[no_integrations$Segment==row$Segment,]$Total<-   no_integrations[no_integrations$Segment==row$Segment,]$Total<-row$Total    
}



pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.adv.boxplot.10000.pdf",sep=""),width=10,height=10)
  p <- ggplot(no_integrations, aes(Segment, delta))+
    xlab("Segment")+
    ylab("FPKM")+
    scale_y_continuous(trans=asinh_trans(),breaks=c(-1000, -100,-10,-1,0,1,10,100, 1000))+    
    geom_boxplot()+
    facet_wrap(~Downstream,ncol=1,scales='free_x')+
    #facet_wrap(Scale~Downstream,scales='free')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p);
  dev.off()  
  
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.adv.boxplot.10000.pdf",sep=""),width=10,height=10)
  p <- ggplot(melted[melted$variable %in% c('delta','Total'),], aes(Segment, value))+
    xlab("Segment")+
    ylab("Value")+
    scale_y_continuous(trans=asinh_trans(),breaks=c(-1000, -100,-10,-1,0,1,10,100, 1000))+    
    geom_boxplot()+
    facet_wrap(variable~Downstream,ncol=2,scales='free_x')+
    #facet_wrap(Scale~Downstream,scales='free')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p);
  dev.off()  
  

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs.fpkm.adv.boxplot.10000.small.pdf",sep=""),width=10,height=6)
  p <- ggplot(no_integrations[no_integrations$delta<100,], aes(Segment, delta))+
    xlab("Segment")+
    ylab("FPKM")+
    #scale_y_continuous(trans=asinh_trans,breaks=c(-100,-50,-10,-1,0,1,10,50,100))+    
    geom_boxplot()+
    facet_wrap(~Downstream,ncol=1,scales='free_x')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p);
  dev.off()  
  
  wgs_fpkm_adv[wgs_fpkm_adv$Segment=='chr8:128737510-128758515' & wgs_fpkm_adv$Integration=='True',]
}

quick_stuff<-function() {
  sim<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/tmp/real.counts.csv', sep='\t',header=FALSE)
  names(sim)<-c('patient','type','Segment','Type','FPKM')
  sim<-sim[sim$FPKM != 0,]
  sim$count<-1
  sim<-unique(sim)
  casted<-cast(melt(sim[,c('Type','count')]))
  casted[casted$Type=='Gene',]$count=6
  casted[casted$Type=='Oncogene',]$count=1
  
  sim<-unique(sim)
  melted<-melt((sim[,c('Type','Segment')]))
  melted$variable=1
  melted<-data.frame(melted[,c('Type','variable')])
  cast<-cast(melt(melted)  )
  cast()
}

viral_other<-function() {
  sim<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/ballgown_noanno/sim.csv', sep='\t',header=FALSE)
  names(sim)<-c('Segment','Patient','Integration','variable','value')
  sim$variable<-paste(sim$variable)
  sim[sim$variable=='ltr',]$variable='LTR'
  sim[sim$variable=='line',]$variable='LINE'
  sim[sim$variable=='gene',]$variable='Gene'
  sim[sim$variable=='onco',]$variable='Oncogene'
  
    
  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/ballgown_noanno/results.csv', sep='\t',header=FALSE)
  names(data)<-c('Segment','Patient','LTR','LINE','Gene','Oncogene','Integration')
  melted<-melt(data)
  
  tmp<-melted[,c('Segment','variable','value')]
  tmp<-melt(cast(tmp,fun=sum))
  tmp<-tmp[tmp$value !=0,]
  tmp$keep<-paste(tmp$Segment,tmp$variable)
  melted$keep<-paste(melted$Segment,melted$variable)
  melted<-melted[melted$keep %in% tmp$keep,]
  
  melted$delta<-melted$value
  integrations<-melted[melted$Integration=='True',]
  for (r in 1:dim(integrations)[1]) {
    row = integrations[r,]
    delta<-mean(integrations[integrations$keep %in% row$keep,]$value)
    melted[melted$keep %in% row$keep,]$delta<-delta-melted[melted$keep %in% row$keep,]$delta
  }
  no_integrations<-melted[melted$Integration=='False',]
  tmp<-cast(melt(no_integrations[,c('Segment','delta')]),fun=mean)
  no_integrations$Segment<-factor(no_integrations$Segment,levels=unique(tmp[order(tmp$delta),]$Segment))
  
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/fpkm.type.pdf",sep=""),width=10,height=4)  
p <- ggplot(no_integrations, aes(x=Segment,y=delta))+
  geom_boxplot()+
  xlab('Segment')+
  ylab('delta FPKM')+
  scale_y_continuous(trans=asinh_trans(),breaks=c(-10000, -1000, -100,-10,-1,0,1,10,100, 1000, 10000))+      
  theme(axis.text.x=element_blank())+  
  facet_wrap(~variable,scales='free')
print(p)
dev.off()
  
}


integration='chr8:129003798-129003989 LTR'
fw <- fitdist(melted[melted$keep==integration & melted$Integration=='False',]$value+0.01, "weibull")
fg <- fitdist(melted[melted$keep==integration & melted$Integration=='False',]$value+0.01, "gamma")
fln <- fitdist(melted[melted$keep==integration & melted$Integration=='False',]$value+0.01, "lnorm")


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/test.pdf",sep=""),width=8,height=8)
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)
dev.off()  

res<-c()
expression<-c()
for (integration in unique(melted$Segment)) {
  for (type in unique(melted$variable)) {
    temp<-melted[melted$Segment==integration & melted$variable==type,]
    if (dim(temp)[1] == 0) {
      next
    }
    if (length(which(temp[temp$Integration=='False',]$value==0))==0) {
      fln <- fitdist(temp[temp$Integration=='False',]$value, "lnorm")
      #fg <- fitdist(temp[temp$Integration=='False',]$value, "gamma")
      #fw <- fitdist(temp[temp$Integration=='False',]$value, "weibull")    
      pvalues = dlnorm(temp[temp$Integration=='True',]$value, meanlog=fln$estimate['meanlog'],sdlog=fln$estimate['sdlog'])
      #pvalues = dweibull(temp[temp$Integration=='True',]$value, shape=fw$estimate['shape'],scale=fln$estimate['scale'])      
    } else {      
      #fn = fitdist(temp[temp$Integration=='False',]$value, "norm")
      #pvalues = dnorm(temp[temp$Integration=='True',]$value, mean=fn$estimate['mean'],sd=fn$estimate['sd'])
      pvalues = dnorm(temp[temp$Integration=='True',]$value, mean=mean(temp[temp$Integration=='False',]$value),sd=sd(temp[temp$Integration=='False',]$value))
    }
    for (value in pvalues) {
      res<-rbind(res, c(integration,type,value))
    }
    expression<-rbind(expression,c(integration,type, mean(temp[temp$Integration=='False',]$value),mean(temp[temp$Integration=='True',]$value)))
  }
}

res<-data.frame(res)
names(res)<-c('Segment','Type','pvalue')
res$pvalue=as.numeric(paste(res$pvalue))  

expression<-data.frame(expression)
names(expression)<-c('Segment','Type','Baseline','Integration')
expression$Baseline=as.numeric(paste(expression$Baseline))  
expression$Integration=as.numeric(paste(expression$Integration))  

casted<-melt(cast(sim[,c('Segment','variable','value')],fun=sum))
casted$keep<-"N"
casted[casted$value != 0,]$keep='Y'
casted$name<-paste(casted$Segment,casted$variable)
sim<-sim[paste(sim$Segment,sim$variable) %in% casted$name,]
sim_res<-c()
sim_expression<-c()
for (integration in unique(sim$Segment)) {
  for (type in unique(sim$variable)) {
    temp<-sim[sim$Segment==integration & sim$variable==type,]
    if (dim(temp)[1] == 0) {
      next
    }
    if (length(which(temp[temp$Integration=='False',]$value==0))==0) {
      fln <- fitdist(temp[temp$Integration=='False',]$value, "lnorm")
      #fg <- fitdist(temp[temp$Integration=='False',]$value, "gamma")
      #fw <- fitdist(temp[temp$Integration=='False',]$value, "weibull")    
      pvalues = dlnorm(temp[temp$Integration=='True',]$value, meanlog=fln$estimate['meanlog'],sdlog=fln$estimate['sdlog'])
      #pvalues = dweibull(temp[temp$Integration=='True',]$value, shape=fw$estimate['shape'],scale=fln$estimate['scale'])      
    } else {      
      #fn = fitdist(temp[temp$Integration=='False',]$value, "norm")
      #pvalues = dnorm(temp[temp$Integration=='True',]$value, mean=fn$estimate['mean'],sd=fn$estimate['sd'])
      pvalues = dnorm(temp[temp$Integration=='True',]$value, mean=mean(temp[temp$Integration=='False',]$value),sd=sd(temp[temp$Integration=='False',]$value))
    }
    for (value in pvalues) {
      sim_res<-rbind(sim_res, c(integration,type,value))
    }
    sim_expression<-rbind(sim_expression,c(integration,type, mean(temp[temp$Integration=='False',]$value),mean(temp[temp$Integration=='True',]$value)))
  }
}

sim_res<-data.frame(sim_res)
names(sim_res)<-c('Segment','Type','pvalue')
sim_res$pvalue=as.numeric(paste(sim_res$pvalue))  
sim_res<-sim_res[sim_res$pvalue != Inf,]

sim_expression<-data.frame(sim_expression)
names(sim_expression)<-c('Segment','Type','Baseline','Integration')
sim_expression$Baseline=as.numeric(paste(sim_expression$Baseline))  
sim_expression$Integration=as.numeric(paste(sim_expression$Integration))  

res$type<-'Real'
sim_res$type<-'Simulated'
joined<-rbind(res[,c('Type','type','pvalue')],sim_res[,c('Type','type','pvalue')])
joined[joined$pvalue>1,]$pvalue=1
for (type in unique(joined$Type)) {
  temp<-joined[joined$Type==type,]
  real<-ecdf(temp[temp$type=='Real',]$pvalue)
  fake<-ecdf(temp[temp$type=='Simulated',]$pvalue)

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pvalues.pdf",sep=""),width=10,height=4)  
p1 <- ggplot(joined[joined$Type!='Oncogene',], aes(x=pvalue+1e-50,color=type))+
  ylim(0,1)+
  xlab('P-value')+
  ylab('CDF')+
  geom_vline(xintercept = 0.05)+
  scale_x_log10()+
  scale_color_discrete(name ='Type')+
  stat_ecdf(geom = "step")+  
  facet_wrap(~Type)
print(p1)
dev.off()
  
}

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pvalues.pdf",sep=""),width=6,height=6)  
p1 <- ggplot(joined, aes(x=type,y=pvalue+1e-300))+
    scale_y_log10()+
    geom_violin()+
    ylab("log(P-value)")+
    facet_wrap(~Type)+
    xlab("Type")
print(p1)
dev.off()

for (type in unique(res$Type)) {  
  temp<-res[res$Type==type,] 
  temp[temp$pvalue <= 1e-100,]$pvalue=1e-100
  temp$Segment<-factor(temp$Segment, levels=temp[order(temp$pvalue),]$Segment)
  t2<-expression[expression$Type==type & expression$Segment %in% temp$Segment,]
  t2$Segment<-factor(t2$Segment, levels=temp[order(temp$pvalue),]$Segment)
#pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/",type,".pvalues.pdf",sep=""),width=10,height=6)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pvalues.pdf",sep=""),width=10,height=6)  
p1 <- ggplot(temp, aes(x=Segment,y=pvalue))+
    scale_y_log10()+
    geom_point()+
    ylab("P-value")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+    
    xlab("Segment")
    melted<-melt(t2)
p2 <- ggplot(melted, aes(x=Segment,y=value,shape=variable,color=variable))+
    geom_point()+
    scale_shape_discrete(name='Type')+
    scale_color_discrete(name='Type')+    
    ylab("mean FPKM")+
    #theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none")+    
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+    
    xlab("Segment")    
  grid.arrange(p1, p2,ncol=1)
  dev.off()
  
}

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pvalues.pdf",sep=""),width=20,height=4)  
p <- ggplot(res, aes(x=Segment,y=pvalue+1e-100))+
    scale_y_log10()+
    facet_wrap(~Type)+        
    ylab("P-value")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+    
    xlab("Segment")

  print(p)
  dev.off()


grid.arrange()

plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)
dev.off()  


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/distribution.pdf",sep=""),width=20,height=4)  
  p <- ggplot(melted, aes(x=Segment,y=value+1))+
    geom_violin()+
    scale_y_log10()+
    facet_grid(.~variable,scales = "free")+        
    ylab("FPKM")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+    
    xlab("Segment")
  print(p)
  dev.off()
}

viral_manuscript<-function() {

  wgs_fpkm_adv<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_expression.csv',sep='\t')
wgs_fpkm_adv$Segment<-paste(wgs_fpkm_adv$Segment)
integrations<-wgs_fpkm_adv[wgs_fpkm_adv$Integration == 'True',]
no_integrations<-wgs_fpkm_adv[wgs_fpkm_adv$Integration == 'False',]
integrations$base_line<-0
for (segment in unique(integrations$Segment)) {
  fpkm = mean(no_integrations[no_integrations$Segment==segment,]$FPKM)
  integrations[integrations$Segment==segment,]$base_line<-fpkm
}


integrations$fold_change<-log((integrations$FPKM+1)/(integrations$base_line+1))/log(2)
mean_int <- (integrations[integrations$Fusion=='True',]$fold_change)
mean_noint <-(integrations[integrations$Fusion=='False',]$fold_change)


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/fpkm_distributions.10000.pdf",sep=""),width=4.5,height=4)
          p <- ggplot(integrations, aes(x=log((FPKM+1)/(base_line+1))/log(2), fill=Fusion)) +
            geom_density(alpha=.3)+
            geom_vline(xintercept = mean(mean_int), size = 1, colour = "skyblue1", linetype = "dashed") +
            geom_vline(xintercept = mean(mean_noint), size = 1, colour = "pink", linetype = "dashed") +               
            xlab(expression(paste("Mean"," ",log[2]-fold," change in FPKM")))+          
            scale_fill_discrete(name='Integration', labels=c('Fusion present', 'Fusionless'), breaks=c('True','False'))+
            theme_bw() +       
            theme(legend.position=c(0.8,0.8), panel.border = element_blank(), panel.grid.major =       
                  element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+                                          
            ylab("Density")
          print(p);
          dev.off()
  
integrations$pvalue<- -1
for (segment in unique(integrations$Segment)) {
  pvalues = 1
  result <- tryCatch(
  {
    fn = fitdist(no_integrations[no_integrations$Segment==segment,]$FPKM, "norm")
    pvalues = min(2*pnorm(integrations[integrations$Segment==segment,]$FPKM, mean=fn$estimate['mean'], sd=fn$estimate['sd'],lower.tail=FALSE),1)
  }, error=function(cond) {
    return(1)
  }
  )
  integrations[integrations$Segment==segment,]$pvalue<-pvalues 
}
integrations$Segment<-factor(integrations$Segment,levels=c(integrations[order(integrations$pvalue),]$Segment))

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/pvalue.10000.pdf",sep=""),width=6,height=4)
          p <- ggplot(integrations, aes(x=Segment, y=log(pvalue), color=Fusion,group=Fusion)) +
            geom_line()+
            geom_point()+
            xlab("Segment")+
            ylab("P-value")+
            scale_color_discrete(name='Integration', labels=c('Fusion present', 'Fusionless'), breaks=c('True','False'))
          print(p);
          dev.off()
  


  times<-read.csv('/pedigree2/projects/namphuon/tmp/results/simulation/times.csv');
  running<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/running.csv', sep=',')  
  
  times<-times[times$Time!=0,]
  times$Integrations<-factor(times$Integrations,levels=c(25,50,100,500,1000))
  times<-rbind(times,c('hg19_SV_viral_int_hpv16_84_95',1000,'ViFi',28820,'Yes'))
  times<-rbind(times,c('hg19_SV_viral_int_hpv16_84_95',1000,'ViFi',28720,'Yes'))
  times<-rbind(times,c('hg19_SV_viral_int_hpv16_84_95',1000,'ViFi',29410,'Yes'))  
  melted<-melt(times[,c('Method','Completed','Time','Integrations')])
  melted$Name<-paste(melted$Method,melted$Completed,sep=',')
  melted[melted$Name=='ViFi,Yes',]$Name = 'ViFi'
  melted[melted$Name=='VERSE,Yes',]$Name = 'VERSE (Completed)'
  melted[melted$Name=='VERSE,No',]$Name = 'VERSE (Failed)'  
  melted$Name<-factor(melted$Name,levels=c('ViFi','VERSE (Completed)','VERSE (Failed)'))
  melted$Time<-as.numeric(melted$Time)
  casted<-cast(melt(melted),fun=mean)
  casted_sd<-cast(melt(melted),fun=sd)
  casted$sd<-casted_sd$Time/60/60
  casted$Integrations<-as.numeric(paste(casted$Integrations))
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/timing.pdf",sep=""),width=6,height=4)
        p <- ggplot(casted, aes(Integrations, Time/60/60,shape=Name,color=Name, fill=Name))+
          geom_errorbar(aes(ymin=Time/60/60-sd, ymax=Time/60/60+sd), width=50) +
          geom_line()+
          geom_point()+
          xlab("Number of Integrations")+
          ylab("Wall clock running time (hr)")
        print(p);
        dev.off()
  


  transcript_count<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/transcript_counts.csv', sep=',')
  
  transcript_count<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_counts.csv', sep=',')
  types<-c('srpRNA','LTR','Satellite','rRNA','DNA','Simple_repeat','scRNA','snRNA','tRNA','RC','RNA','LINE','SINE','Retroposon','Low_complexity','Gene','Oncogene')
  transcript_count<-transcript_count[transcript_count$Type %in% types,]
  transcript_count[grep('RNA',transcript_count$Type),]$Type<-'RNA'
  
  all_counts<-read.csv('/pedigree2/projects/namphuon/tmp//results/TCGA-CESC/analyses/counts_type.all1.csv')
  sub<-all_counts[all_counts$Type=='Oncogene',]
  sub$Replicate<-paste(sub$Replicate)
  sub[sub$Replicate=='Real',]$Replicate<-'Observed'
  merged<-rbind(sub[,c('Type','Count','Replicate')],transcript_count[,c('Type','Count','Replicate')])
  transcript_count<-merged
  casted<-cast(melt(transcript_count[,c('Type','Count','Replicate'),]),fun=sum)
  casted$Experiment<-'Expected'
  casted[casted$Replicate=='Observed',]$Experiment<-'Observed'
  sds<-cast(melt(data.frame(casted[,c('Type','Count','Experiment')])),fun=sd)
  #sds[is.na(sds$Count),]$Count<-0
  means<-cast(melt(data.frame(casted[,c('Type','Count','Experiment')])),fun=mean)
  names(sds)<-c('Type','Experiment','SD')
  merged<-merge(means,sds)

  casted<-cast(melt(transcript_count[,c('Type','Count','Replicate'),]),fun=sum)
  casted$Experiment<-'Expected'
  casted[casted$Replicate=='Observed',]$Experiment<-'Observed'
  sim<-casted[casted$Experiment=='Expected',]
  real<-casted[casted$Experiment=='Observed',]
  types <- c('SINE','LINE','LTR','DNA','Gene','Oncogene')
  sim<-sim[sim$Type %in% types,]  
  sim$Type<-factor(sim$Type,levels=types)

  
  
#pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expected.rna.counts.violin.pdf",sep=""),width=8,height=4)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expected.wgs.counts.violin.pdf",sep=""),width=8,height=4)  
  p1 <- ggplot(sim, aes(x=Type,y=Count))+
      facet_wrap(~Type,scales='free')+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
      geom_point(data=real[real$Type %in% types,],aes(x = Type, y = Count,color='red'))+
      ylab("Number of expressed genomic segments\n of each functional types")+    
      theme(legend.position="none")+        
      #scale_y_log10()+        
      xlab("Type")
  print(p1)
  dev.off()

#pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expected.rna.counts.violin.all.pdf",sep=""),width=8,height=4)  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/expected.wgs.counts.violin.all.pdf",sep=""),width=8,height=4)  
  p1 <- ggplot(sim, aes(x=Type,y=Count))+
      facet_wrap(~Type,scales='free')+
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
      geom_point(data=real,aes(x = Type, y = Count,color='red'))+
      ylab("Number of expressed genomic segments\n of each functional types")+    
      theme(legend.position="none")+        
      #scale_y_log10()+        
      xlab("Type")
  print(p1)
  dev.off()
  
for (type in unique(casted$Type)) {
  fn <- fitdist(casted[casted$Experiment=='Expected' & casted$Type == type,]$Count, "norm")
  pvalues = pnorm(casted[casted$Experiment=='Observed' & casted$Type == type,]$Count, mean=fn$estimate['mean'],sd=fn$estimate['sd'])  
  print(paste(type,pvalues))  
}

  running<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/running.csv', sep=',')
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/running.pdf",sep=""),width=4,height=4)  
  p <- ggplot(running, aes(x=Integrations,y=Time/60/60,group=Method,color=Method,shape=Method))+
    geom_point()+
    scale_color_discrete(name='Methods', breaks=c('ViFi', 'ViFi-NoHMM', 'VERSE'), label=c('ViFi', 'ViFi-NoHMM', 'VERSE'))+    
    scale_shape_discrete(name='Methods', breaks=c('ViFi', 'ViFi-NoHMM', 'VERSE'), label=c('ViFi', 'ViFi-NoHMM', 'VERSE'))+    
    ylab("Running time (hours)")+
    xlab("Number of integrations")
  print(p)
  dev.off()
  
    p <- ggplot(melted, aes(x=Method, y=value,fill=factor(variable)))+
    geom_bar(stat='identity')+
    #scale_fill_discrete(breaks=c('our_pipeline','virus_finder_2', 'vfs'), label=c('Our Method', 'VERSE', 'ViralFusionSeq'))+
    scale_fill_discrete(name='Integration sites', breaks=c('TP','Total'), label=c('True Positive','False Negative'))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Method")+
    scale_x_discrete(labels=c("our_pipeline" = "ViFi", "virus_finder_2" = "VERSE"))+
    facet_wrap(~Model)+    
    ylab("Integration Sites")
  print(p)
  dev.off()

  

  library(ballgown)
  meta_info<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/biological.metadata.info.csv', header=FALSE,sep=',')
  names(meta_info)<-c('patient','integration','genes','exons','oncogenes')
  toremove<-grep('chr2:13303',meta_info$integration)
  meta_info <- meta_info[-toremove,]
  my_list<-data.frame(t(sapply(strsplit(paste(meta_info$integration), ":|-"),c)))
  names(my_list)<-c('chr','start','end')
  my_list$start<-as.numeric(paste(my_list$start))
  my_list$end<-as.numeric(paste(my_list$end))

  meta_info<-cbind(meta_info,my_list)  
  gene_less<-meta_info[meta_info$genes==0,]
  
  dir='/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/ballgown'
  dir='/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/ballgown_noanno_rescored'
  #setwd(dir)
  data_directory = system.file(dir, package='ballgown')
  bg = ballgown(dataDir=dir, samplePattern='TCGA', meas='all')
  save(bg, file='bg_noanno_rescored.rda')
  load('bg.rda')
  
  meta_map <-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/gffcompared.gtf.merged.gtf.tmap',sep='\t')
  class_map <-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/repeats/repeats.types.csv',sep='\t',header=FALSE)
  names(class_map)<-c('ref_gene_id','class')
  merged<-merge(meta_map,class_map)
  
  gene_map <- read.csv('/pedigree2/projects/namphuon/data/references/hg19/annotations/genes.gtf', header=FALSE, sep='\t')
  names(gene_map)<-c('chr','annotation','type','start','end', 'unknown1', 'strand', 'unknown2', 'meta')

  whole_tx_table = texpr(bg, 'all')
  #write.csv(whole_tx_table,file='merged.csv',quote=FALSE ,row.names = FALSE)
  all<-merge(whole_tx_table,merged[,c('qry_gene_id','class')],by.x=c('gene_id'),by.y=c('qry_gene_id'))  
  #all<-whole_tx_table
  #
  
  #Reduce to lines and ltrs
  subset<-all[all$class %in% c('LINE','LTR'),]
  subset$name<-paste(subset$chr,paste(subset$start,subset$end,sep='-'),sep=':')
    
  #Iterate through patients
  df<-c()
  meta_info<-gene_less
  for (patient in unique(meta_info$patient)) {
    print(patient)
    if (paste("FPKM.",patient,sep='') %in% colnames(subset) == FALSE) {
      next
    }
    if (patient %in% unique(df$patient) == TRUE) {
      next
    }
    
     temp<-subset[,c('chr', 'start', 'end', 'class', paste("FPKM.",patient,sep=''))]
     names(temp)<-c('chr', 'start', 'end', 'class', 'FPKM')
#    temp<-subset[,c('chr', 'start', 'end', paste("FPKM.",patient,sep=''))]
#    names(temp)<-c('chr', 'start', 'end', 'FPKM')
    
    temp$patient<-patient
    temp$integration<-'N'
    for (row in which(meta_info$patient==patient)) {
      ints<-which(temp$chr==gene_less[row,]$chr & (temp$start-10000)<=meta_info[row,]$start & (temp$start+10000)>=meta_info[row,]$end)
      if (length(ints) != 0) {
        temp[ints,]$integration<-'Y'      
      }
    }
    df<-rbind(df,temp)
  }
  #df<-df[grep('chr',df$chr),]
  save(df, file='df.10000.rda')
  load('df.10000.rda')

  df$name<-paste(df$chr,paste(df$start,df$end,sep='-'),sep=':')
  temp<-df
  temp$value<-temp$FPKM
  temp$variable<-'FPKM'
  temp$FPKM<-NULL    
  casted<-cast(temp,fun=sum)
  temp<-casted
  
  keeps<-unique(temp[temp$integration=='Y',]$name)
  mini_temp<-temp[temp$name %in% keeps,]
  #mini_temp<-mini_temp[mini_temp$name %in% paste(gene_less$integration),]
  tmp<-cast(melt(as.data.frame(test),id.vars=c('name')),fun=mean)
  mini_temp$name<-factor(mini_temp$name, levels=tmp[order(tmp$FPKM),]$name)
  for (name in unique(mini_temp$name)) {
    if (dim(mini_temp[mini_temp$name == name,])[1] == 0) {
      next
    }
    r=mini_temp[mini_temp$name == name,][1,]
    mini_temp[mini_temp$chr==r$chr & pmax(mini_temp$start,r$start) <= (pmin(mini_temp$end,r$end)+100000),]$name<-r$name    
  }
  
  new_temp<-c()
  for (r in 1:dim(gene_less)[1]) {
    r=gene_less[r,]
    new_temp<-unique(rbind(new_temp,mini_temp[which(mini_temp$chr==r$chr & pmax(mini_temp$start,r$start) <= (pmin(mini_temp$end,r$end)+100000)),]))
  }
  temp<-cast(melt(as.data.frame(new_temp[,c('name','FPKM')])),fun=mean)
  new_temp$name<-factor(new_temp$name, levels=temp[order(temp$FPKM),]$name)
  
    pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeats.fpkm.boxplot.pdf",sep=""),width=12,height=6)
  #pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/noannotation.fpkm.boxplot.pdf",sep=""),width=12,height=6)  
        p <- ggplot(mini_temp, aes(name, FPKM+1,shape=integration,color=integration, fill=integration))+
          #geom_point(position = position_jitter(width = 0.2))+
          xlab("Integration")+
          ylab("FPKM")+
          scale_y_log10()+
          geom_boxplot()+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))+
          facet_wrap(~class,scales = "free")          
        print(p);
        dev.off()
    
  keeper<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/ballgown/test.csv.keep.1',header=FALSE)
  names(keeper)<-c('name')
  mini_df<-df[df$name %in% keeper$name,]
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeats.fpkm.violin.pdf",sep=""),width=16,height=16)  
      p <- ggplot(mini_df, aes(integration, FPKM+1))+
        xlab("Integration")+
        ylab("log(FPKM)")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~class,scales = "free_y")
      print(p);
      dev.off()	  
  
    #Statistical comparison
    for (class in unique(mini_temp$class)) {
      casted<-cast(melt(mini_temp[mini_temp$class==class,c('name','FPKM','integration')]),name~integration,fun=median)    
      print(class)
      print(wilcox.test(x=casted$N,y=casted$Y,paired=TRUE))
    }
    


  mini_df<-df[df$patient=='TCGA-C5-A0TN',]
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeat.fpkm.1.pdf",sep=""),width=16,height=16)  
      p <- ggplot(df, aes(class, FPKM+1, fill=integration,color=integration))+
        #geom_point(position = position_jitter(width = 0.2))+
        xlab("Class")+
        ylab("FPKM")+
        scale_y_log10()+
        geom_boxplot()+
        facet_wrap(~patient,scales = "free_y")
      print(p);
      dev.off()

      
  
  
  #Now identify all lines/ltrs within 200kb of integration, compare expression versus all lines/ltrs within sample
  df<-c()
  for (row in 1:dim(meta_info)[1]) {
    keep<-which(subset$chr==meta_info[row,]$chr & (subset$start-10000)<=meta_info[row,]$start & (subset$start+10000)>=meta_info[row,]$end)
    within<-union(within,keep)
  }
  subset$integration<-'N'
  subset[within,]$integration<-'Y'
  
  temp<-subset[,union(c('class','integration'),colnames(subset)[grep('FPKM',colnames(subset))])]
  

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeat.fpkm.pdf",sep=""),width=16,height=16)  
  
  plotTranscripts('MSTRG.32435', bg, 
    meas='FPKM', colorby='transcript')
  dev.off()  
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeat.fpkm.pdf",sep=""),width=16,height=16)  
      p <- ggplot(temp, aes(class, distance+1))+
        xlab("Type")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~type,scales = "free_y")
      print(p);
      dev.off()	  
    
  
  
repeats<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/repeats/distances.type.csv',header=FALSE,sep=',')
  names(repeats)<-c('type','distance')
  repeats$class<-'type'
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeats.types.violin.pdf",sep=""),width=16,height=16)  
      p <- ggplot(repeats, aes(class, distance+1))+
        xlab("Type")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~type,scales = "free_y")
      print(p);
      dev.off()	  

repeats<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/repeats/distances.class.csv',header=FALSE,sep=',')
  names(repeats)<-c('type','length','distance')
  repeats$class<-'type'
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeats.class.violin.pdf",sep=""),width=16,height=16)  
      p <- ggplot(repeats, aes(type, distance+1))+
        xlab("Type")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()
      print(p);
      dev.off()
      
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/repeats.class.cdf.pdf",sep=""),width=16,height=6)  
      p <- ggplot(repeats, aes(x=distance+1,color=type))+
        xlab("Probability of type with distance")+
        ylab("NT distance")+
        scale_x_log10()+
        stat_ecdf(geom = "step")
      print(p);
      dev.off()	  
      	  
  

  distances<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/biological_repeat.csv.fixed',header=FALSE,sep=',')
  names(distances)<-c('patient','integration','type','distance','experiment','rep')

toremove<-grep('chr2:13303',distances$integration)
distances <- distances[-toremove,]


meta_info<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/biological.metadata.info.csv',header=FALSE,sep=',')
  names(meta_info)<-c('patient','integration','genes','exons','oncogenes')
  toremove<-grep('chr2:13303',meta_info$integration)
  meta_info <- meta_info[-toremove,]
    
  integrations<-unique(meta_info[meta_info$genes==0,]$integration)
  patients<-unique(meta_info[meta_info$oncogenes==0,]$patient)
  temp_meta<-meta_info
  temp_meta$integration<-NULL
  temp_meta<-cast(melt(temp_meta),fun=sum)  
  patients<-unique(temp_meta[temp_meta$oncogenes==0,]$patient)
  distances<-distances[distances$patient %in% patients,]
  
  all_types = unique(distances$type)
  keepers = c('LINE','LTR','ERV')
  distances$rep<-factor(distances$rep, levels=c(unique(distances$rep)))  
  
  res = c()
  for (keep in keepers) {
    res=union(res, grep(keep, all_types))
  }
  subdistances<-distances[!is.na(distances$distance) & distances$type %in% all_types[res],]
  
  #Print all integrations
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.all.types.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(subdistances, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~type,scales = "free_y")
      print(p);
      dev.off()	  

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.all.types.cdf.pdf",sep=""),width=16,height=6)  
      p <- ggplot(subdistances, aes(x=distance+1,color=experiment))+
        xlab("Probability of type with distance")+
        ylab("NT distance")+
        scale_x_log10()+
        stat_ecdf(geom = "step")+
        facet_wrap(~type,scales = "free_x")
      print(p);
      dev.off()	  

  #group by class type
  temp<-subdistances
  temp$class<-""
  keepers = c('LINE','LTR')  
  for (keep in keepers) {
    temp[grep(keep, temp$type),]$class<-keep
  }
  temp$type<-NULL
  casted<-cast(melt(temp),fun=min)
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.all.class.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~class,scales = "free_y")
      print(p);
      dev.off()	  
      
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.all.class.cdf.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(x=distance+1,color=experiment))+
        xlab("Probability of type with distance")+
        ylab("NT distance")+
        scale_x_log10()+
        stat_ecdf(geom = "step")+
        facet_wrap(~class,scales = "free_x")
      print(p);
      dev.off()	  
      
  
  #Now reduce to closest integration to type
  temp<-subdistances
  temp$integration<-NULL
  casted<-cast(melt(temp),fun=min)
  pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.patient.type.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~type,scales = "free_y")
      print(p);
      dev.off()	  
      
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.patient.type.cdf.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(x=distance+1,color=experiment))+
        xlab("Probability of type with distance")+
        ylab("NT distance")+
        scale_x_log10()+
        stat_ecdf(geom = "step")+
        facet_wrap(~type,scales = "free_x")
      print(p);
      dev.off()	  
      
    
  #Now group by class
  temp<-subdistances
  temp$class<-""
  keepers = c('LINE','LTR')  
  for (keep in keepers) {
    temp[grep(keep, temp$type),]$class<-keep
  }
  temp$type<-NULL
  temp$integration<-NULL
  casted<-cast(melt(temp),fun=min)
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.patient.class.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~class,scales = "free_y")
      print(p);
      dev.off()	 
        
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.patient.class.cdf.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(x=distance+1,color=experiment))+
        xlab("Probability of type with distance")+
        ylab("NT distance")+
        scale_x_log10()+
        stat_ecdf(geom = "step")+
        facet_wrap(~class,scales = "free_x")
      print(p);
      dev.off()	  

  temp_meta_info<-meta_info
  temp_meta_info$counter<-1
  integration_count<-cast(melt(temp_meta_info[,c('patient','counter')]))
  merged<-merge(integration_count,casted)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.integrations.patient.class.violin.pdf",sep=""),width=20,height=6)  
      p <- ggplot(merged, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_grid(class~counter,scales = "free_y")
      print(p);
      dev.off()	  
  
  temp$class<-""
  keepers = c('LINE','LTR')  
  for (keep in keepers) {
    temp[grep(keep, temp$type),]$class<-keep
  }
  
  
    
  temp<-subdistances[subdistances$patient %in% patients & subdistances$integration %in% integrations,]
  temp$integration<-NULL
  temp$class<-""
  keepers = c('LINE','ERV')  
  for (keep in keepers) {
    temp[grep(keep, temp$type),]$class<-keep
  }
  temp$type<-NULL
  casted<-cast(melt(temp),fun=min)

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.all.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~type,scales = "free_y")
      print(p);
      dev.off()	  

  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.single.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        scale_y_log10()+
        geom_violin()+
        facet_wrap(~class,scales = "free_y")
      print(p);
      dev.off()	  


  pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.violin.pdf",sep=""),width=16,height=6)  
      p <- ggplot(subdistances, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        geom_violin()+
        scale_y_log10()+
        facet_wrap(~type,scales = "free_y")
      print(p);
      dev.off()	  
      
  subdistances$class<-""
  keepers = c('LINE','LTR')
  
  for (keep in keepers) {
    subdistances[grep(keep, subdistances$type),]$class<-keep
  }
  temp<-subdistances
  temp$type<-NULL
  casted<-cast(melt(temp),fun=min)
   
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.violin.types.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        geom_violin()+
        scale_y_log10()+
        facet_wrap(~class,scales = "free_y")
      print(p);
      dev.off()	  
      
   subdistances$class<-""

  temp$class<-NULL
  casted<-cast(melt(temp, id=c('patient','integration','experiment')),fun=min)

pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/distances.violin.all.pdf",sep=""),width=16,height=6)  
      p <- ggplot(casted, aes(experiment, distance+1))+
        xlab("Integration")+
        ylab("Distance to type")+
        geom_violin()+
        scale_y_log10()
      print(p);
      dev.off()	  
      
   subdistances$class<-""

  keepers = c('LINE','LTR')
  
  for (keep in keepers) {
    subdistances[grep(keep, subdistances$type),]$class<-keep
  }

       
  

  meta_wgs<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_meta.gene.csv',sep=',')


  meta_wgs<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/wgs_meta.gene.csv',sep=',')
  
  #Find number of patients with no integrations
  total_matched_samples = 68
  has_integration = length(unique(meta_wgs$patient))
  has_fusion_mrna_expression = length(unique(meta_wgs[!is.na(meta_wgs$gene),]$patient))
  total_wgs = dim(meta_wgs)[1]
  total_fusion_mrna_expression = dim(meta_wgs[!is.na(meta_wgs$gene),])[1]
  oncohits = meta_wgs[meta_wgs$oncogene=='Y' & !is.na(meta_wgs$oncogene),]
  melt(oncohits[,c('patient', 'gene_id')])
  
  
  mrna_expression<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/mrna_expression.csv',sep=',')
  df<-cast(melt(mrna_expression[,c('patient','coverage')]),fun=sum)  
  rownames(df)<-df$patient  
  mrna_expression$nc<-(mrna_expression$coverage/df[mrna_expression$patient,]$coverage)
  
  df<-cast(melt(mrna_expression[,c('location','nc')]),fun=mean)  
  rownames(df)<-df$location
  temp<-cast(melt(mrna_expression[,c('location','nc')]),fun=sd)
  names(temp)<-c('location','sd')  
  df<-merge(df,temp)  
  rownames(df)<-df$location
  
  mrna_expression$pvalue<-((mrna_expression$nc-df[mrna_expression$location,]$nc)/df[mrna_expression$location,]$sd)
    

  mrna_expression<-mrna_expression[mrna_expression$location %in% rownames(df),]
  mrna_expression$size<-10
  mrna_expression[mrna_expression$is_integration=='Y',]$size<-20
  mrna_expression<-droplevels(mrna_expression)
  
  sorter<-cast(melt(mrna_expression[mrna_expression$is_integration=='Y',c('location','pvalue')]),fun=max)  
  mrna_expression$location<-factor(mrna_expression$location, levels=sorter[order(sorter$pvalue),]$location)
  pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/mrna_expression.pdf",sep=""),width=16,height=6)  
    p <- ggplot(mrna_expression, aes(location, pvalue, shape=is_integration,color=is_integration))+
      xlab("Location")+
      ylab("SD from mean expression")+
      geom_point(position = position_jitter(width = 0.2))+
      scale_shape_discrete(name='Integration Point')+
      scale_colour_discrete(name='Integration Point')+
      #scale_size(breaks=c(10,20), labels=c('N','Y'), name='Integration Point')+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
      #geom_boxplot(colour='black')
    print(p);
    dev.off()	  
    
  mrna_expression$size = 1
  mrna_expression[mrna_expression$is_integration=='Y',]$size = 5
  
  sorter<-cast(melt(mrna_expression[mrna_expression$is_integration=='Y',c('location','nc')]),fun=max)  
  mrna_expression$location<-factor(mrna_expression$location, levels=sorter[order(sorter$nc),]$location)
  pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/mrna_expression.pdf",sep=""),width=16,height=6)  
    p <- ggplot(mrna_expression, aes(location, nc, shape=is_integration,color=is_integration, group=location))+
      xlab("Location")+
      ylab("Normalized coverage")+
      geom_point(position = position_jitter(width = 0.2), aes(size=size))+
      scale_shape_discrete(name='Integration Point')+
      scale_colour_discrete(name='Integration Point')+
      #scale_size(breaks=c(10,20), labels=c('N','Y'), name='Integration Point')+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      geom_boxplot(colour='black')
    print(p);
    dev.off()	  
    

  sorter<-cast(melt(mrna_expression[mrna_expression$is_integration=='Y',c('location','pvalue')]),fun=max)  
  mrna_expression$location<-factor(mrna_expression$location, levels=sorter[order(sorter$pvalue),]$location)
  pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/mrna_expression.pvalue.pdf",sep=""),width=16,height=6)  
    p1 <- ggplot(mrna_expression, aes(location, pvalue, shape=is_integration,color=is_integration, group=location))+
      xlab("Location")+
      ylab("SD from mean expression")+
      geom_point(position = position_jitter(width = 0.2), aes(size=size))+
      scale_shape_discrete(name='Integration Point')+
      scale_colour_discrete(name='Integration Point')+
      #scale_size(breaks=c(10,20), labels=c('N','Y'), name='Integration Point')+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      geom_boxplot(colour='black')
    print(p);
    dev.off()  
    
    p2 <- ggplot(unique(mrna_expression[,c('location','genes_within_100kb')]), aes(location, genes_within_100kb))+
      xlab("Location")+
      ylab("Genes within 100kB")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      geom_point()
    print(p2);
    dev.off() 
    
library(gridExtra)
maxWidth = grid::unit.pmax(p1$widths[2:3], p2$widths[2:3])
grid.newpage()

pdf('/pedigree2/projects/namphuon/tmp/results/simulation/analyses/mrna_expression.test.pdf',
     width=8,
     height=6)
grid.arrange(
    arrangeGrob(p1,p2,nrow=2,heights=c(0,10))
    )
dev.off()

    grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2)))
    
      
  
  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/simulations.csv',sep=",", header=FALSE)
  names(data)<-c('model','coverage','method','threshold','recall','precision')
  
  data$model<-factor(data$model,levels=c("agpv1", "hpv16_123_90", "hpv16_147_95", "hpv16_17_99", "hpv16_29_95", "hpv16_40_99", "hpv16_65_90", "hpv16_84_95", "hpv16_103_95", "hpv16_125_90", "hpv16_15_99", "hpv16_26_90", "hpv16_35_95", "hpv16_60_99", "hpv16_69_99", "hpv16_9_90"))
  
  data$case<-'Brown Howler HPV (65% similarity)'
  data[grep('99',data$model),]$case<-'Easy (99% similarity)'
  data[grep('95',data$model),]$case<-'Medium (95% similarity)'
  data[grep('90',data$model),]$case<-'Hard (90% similarity)'
  data$case<-factor(data$case, levels=c('Easy (99% similarity)','Medium (95% similarity)','Hard (90% similarity)','Brown Howler HPV (65% similarity)'))
  
  data$coverage<-factor(data$coverage,levels=c(25,5))
  melted = melt(data, id=c('model','coverage','method','threshold','case'))
  for (model in unique(data$model)) {
    for (method in unique(data$method)) {
      for (coverage in unique(data$coverage)) {
        if (dim(data[data$model==model & data$method==method,])[1]==0) {
          paste(model,method)
        }    
      }
    }
  }
  casted <- cast(melted[,!(names(melted) %in% c('model'))], fun=mean)
  #matches<-regexpr('case', casted$model)
  #casted$type<-toTitleCase(gsub("_", " ", substr(casted$model,1,matches+3)))
  #casted$type<-factor(casted$type, levels=c('Easy Case', "Med Case", "Hard Case"))
  #casted<-casted[!casted$method=='vfs',]
  casted$recall<-as.numeric(paste(casted$recall))
  casted$precision<-as.numeric(paste(casted$precision))
  
  subcast<-casted[casted$method %in% c('method','method_no_hmm') & casted$threshold == 4,]
  test<-subcast
  test<-rbind(test,casted[casted$method == 'virus_finder_2',])
      pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/simulation.1000.pdf",sep=""),width=6,height=4)  
  p <- ggplot(test[!(test$recall == 0 & test$precision == 0) & test$threshold <=4 & test$coverage==25,], aes(x=precision, y=recall, color=method,shape=method))+
    geom_point()+
    facet_wrap(~case)+
    scale_color_discrete(name='Methods', breaks=c('method','method_no_hmm','virus_finder_2'), label=c('ViFi', 'ViFi-NoHMM', 'VERSE'))+    
    scale_shape_discrete(name='Methods', breaks=c('method','method_no_hmm','virus_finder_2'), label=c('ViFi', 'ViFi-NoHMM', 'VERSE'))+    
    #scale_shape_discrete(name='Coverage', breaks=c(5,25), label=c('5','25'))+        
    xlab("Precision")+
    xlim(0,1)+
    ylim(0,1)+
    ylab("Recall")
  print(p)
  dev.off()


pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/simulation.pdf",sep=""),width=8,height=5)  
  p <- ggplot(casted[!(casted$recall == 0 & casted$precision == 0) & casted$threshold <=4 & casted$case %in% c('Easy (99% similarity)','Brown Howler HPV (65% similarity)'),], aes(x=precision, y=recall,color=method))+
    geom_point()+
    facet_grid(coverage~case)+
    scale_color_discrete(name='Methods', breaks=c('method','method_no_hmm','virus_finder_2'), label=c('Default', 'Default-NoHMM', 'VERSE'))+    
    xlab("Precision")+
    xlim(0,1)+
    ylim(0,1)+
    ylab("Recall")
  print(p)
  dev.off()

  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/HCC/analyses/hcc.csv',sep=",", header=FALSE)
  names(data)<-c('Model','Dataset','Method','TP','Total')
  df <- cast(melt(data[,c('Model','Method','TP','Total')]),fun=sum)
  df$Total<-df$Total-df$TP
  melted<-melt(df)
  melted$variable<-factor(melted$variable,levels=c('Total','TP'))
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/hcc.pdf",sep=""),width=6,height=4)  
  p <- ggplot(melted, aes(x=Method, y=value,fill=factor(variable)))+
    geom_bar(stat='identity')+
    #scale_fill_discrete(breaks=c('our_pipeline','virus_finder_2', 'vfs'), label=c('Our Method', 'VERSE', 'ViralFusionSeq'))+
    scale_fill_discrete(name='Integration sites', breaks=c('TP','Total'), label=c('True Positive','False Negative'))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Method")+
    scale_x_discrete(labels=c("our_pipeline" = "ViFi", "virus_finder_2" = "VERSE"))+
    facet_wrap(~Model)+    
    ylab("Integration Sites")
  print(p)
  dev.off()
  
  venn<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/tcga_results.csv',sep=",", header=TRUE)
  venn$patient<-NULL
  
  df<-cast(melt(venn),fun=sum)
  pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/pdfs/tcga_triple.pdf",sep=""),width=5,height=4)  
  venn.plot <- draw.triple.venn(df$n1, df$n2, df$n3, df$n12, df$n23, df$n13, df$n123, c("ViFi", "VERSE", "Tang et. al Study"),margin=0.05,cat.pos=180);
  dev.off()
  
  overall<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/tcga_results.detailed.csv',sep=",", header=TRUE)
  overall$patient<-NULL  
  abc<-cast(melt(overall),fun=sum)
  melted<-melt(abc)
  melted$method<-factor(melted$method, levels=c('Our Method','VERSE','Larsson'))
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/tcga_overall.pdf",sep=""),width=4,height=4)  
  p <- ggplot(melt(melted), aes(x=method, y=value,fill=variable))+
    geom_bar(stat='identity',position="dodge")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_discrete(name='Within 100kb\n of WGS integration',breaks=c('tp','fp'), labels=c('Yes','No'))+
    xlab("Method")+    
    ylab("RNASeq Inegrations")
  print(p)
  dev.off()  
  
  venn<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/tcga_results.super.csv',sep=",", header=TRUE)
  venn$patient<-NULL
  
  df<-cast(melt(venn),fun=sum)
pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/tcga_quad.pdf",sep=""),width=4,height=4)  
  venn.plot <- draw.quad.venn(df$n1, df$n2, df$n3, df$n4, df$n12, df$n13, df$n14, df$n23, df$n24, df$n34, df$n123, df$n124, df$n134, df$n234, df$n1234, c("Our Method", "VERSE", "Larsson Study", "WGS-supported"));
  dev.off()  
  
  
  time<-data.frame(time=c(41.81,177.23),method=c('Our Method','VERSE'))
pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/tcga_time.pdf",sep=""),width=4,height=4)  
  p <- ggplot(melt(time), aes(x=method, y=value))+
    geom_bar(stat='identity',position="dodge")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Method")+    
    ylab("Wall clock running time (hr)")
  p
  dev.off()  
      
  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/HCC/analyses/tcga.csv.300',sep=" ", header=FALSE)
  names(data)<-c('Cancer','Patient','n1','n2','n3','n12','n23','n13','n123')
  data$Patient<-NULL
  data$Cancer<-NULL
  
  df<-cast(melt(data),fun=sum)
  pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/tcga.300.pdf",sep=""),width=4,height=4)  
  venn.plot <- draw.triple.venn(df$n1, df$n2, df$n3, df$n12, df$n23, df$n13, df$n123, c("WGS Sites", "Larsson Study", "Our Method"));
  dev.off()
  
  data<-read.csv('test.csv',sep=" ", header=FALSE)
  names(data)<-c('read','position','coverage')
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/simulation/analyses/test.pdf",sep=""),width=6,height=4)  
  p <- ggplot(data, aes(x=position, y=coverage))+
    geom_point()+
    xlab("Position")+
    ylab("Coverage")+
    facet_wrap(~read,ncol=1)
  print(p)
  dev.off()
  
  data<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/patients.csv',sep=",", header=FALSE)
  names(data)<-c('patient','integration','distance','true_integration')
  data[is.na(data$distance),]$distance<- 1000000000
  
  
  meta<-read.csv('/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/metainfo.csv',sep=",", header=FALSE)
  names(meta)<-c('patient','integration','distance','type','gene')
  meta$type<-paste(meta$type)
  meta[is.na(meta$distance),]$distance<- 1000000000
  rownames(meta)<-paste(meta$patient,meta$integration)
  
  
  keep = intersect(data$integration,meta$integration)
  subdata<-data[data$integration %in% keep,]
  submeta<-meta[meta$integration %in% keep,]  
  
  merged <- merge(subdata,submeta,by=c('integration','patient'))
  merged[merged$distance.x==1000000000,]$type<-"FP"
  merged$gene_y<-meta[paste(merged$patient,merged$true_integration),]$gene
  merged$same_gene<-merged$gene_y==merged$gene
  
  data[data$distance==1000000000,]$distance<- -500
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/tcga_results.pdf",sep=""),width=6,height=4)  
  p <- ggplot(data, aes(x=distance))+
    geom_histogram(breaks=seq(-500,10000,500))+
    xlab("Distance to nearest WGS integration point")
  print(p)
  dev.off()
  
pdf(paste("/pedigree2/projects/namphuon/tmp/results/TCGA-CESC/analyses/tcga_scatter.pdf",sep=""),width=6,height=4)  
  p <- ggplot(merged, aes(x=distance.x,distance.y,shape=type,color=same_gene))+
    geom_point()+
    scale_shape_discrete(name='Annotation type')+
    scale_color_discrete(name='Gene annoation match')+    
    xlab("Distance to nearest WGS integration point")+
    ylab("Distance to nearest mRNA annotation")+
    scale_x_log10()+
    scale_y_log10()
  print(p)
  dev.off()
  
  
}

venn_diagram<-function() {
  library(VennDiagram)
  library(dict)
  library(grid)
  library(gridBase)
  
  
  data<-read.csv('/pedigree2/projects/namphuon/results/hbv/analyses/tree_analyses/venn.csv',sep=",", header=TRUE)
  for (support in unique(data$support)) {
    df<-data[data$support==support,c('type','n1','n2','n3','n12','n23','n13','n123')]
    test<-list()
    for (t in df$type) {
      casted<-df[df$type==t,]
      #test[[t]] <- draw.triple.venn(casted$n1, casted$n2, casted$n3, casted$n12, casted$n23, casted$n13, casted$n123, c("Genomes_Manual", "Genomes_SRegion", "Genomes_PASTA"),scaled=TRUE,euler.d=TRUE);    
    g<-draw.triple.venn(casted$n1, casted$n2, casted$n3, casted$n12, casted$n23, casted$n13, casted$n123, c("Genomes_Manual", "Genomes_SRegion", "Genomes_PASTA"),scaled=TRUE,euler.d=TRUE,fill = c('red','blue','green'),margin=0.25,main=paste('Subtype:',t," Support:",support,sep=''),main.cex=2);   
      
pdf(paste("/pedigree2/projects/namphuon/results/hbv/analyses/tree_analyses/venn.",t,'.',support,'.pdf',sep=""),width=5,height=5)   
      grid.arrange(gTree(children=g), top=paste('Subtype:',t," Support:",support,sep=''))
      dev.off() 
      
    }
}    
    casted<-cast(melt(df),fun=sum)
      pdf(paste("/pedigree2/projects/namphuon/results/hbv/analyses/tree_analyses/venn.",support,'.pdf',sep=""),width=5,height=5)  
      
    
    
    venn.plot <- draw.triple.venn(casted$n1, casted$n2, casted$n3, casted$n12, casted$n23, casted$n13, casted$n123, c("Genomes_Manual", "Genomes_SRegion", "Genomes_PASTA"),scaled=TRUE,euler.d=TRUE);
     dev.off()
}
    
  layout(matrix(1:8))

pdf(paste("/pedigree2/projects/namphuon/results/hbv/analyses/tree_analyses/venn.",support,'.pdf',sep=""),width=5,height=4)  
    grid.arrange(grobTree(types[['A']]),grobTree(types[['B']]),grobTree(types[['C']]),grobTree(types[['D']]),grobTree(types[['E']]),grobTree(types[['F']]),grobTree(types[['G']]),grobTree(types[['I']]),ncol=4)
    dev.off()
    grid.draw(types[['A']])
    frame()
    grid.draw(types[['B']])
    frame()
    grid.draw(types[['C']])
    frame()
    grid.draw(types[['D']])
    frame()
    grid.draw(types[['E']])
    frame()
    grid.draw(types[['F']])
    frame()
    grid.draw(types[['G']])
    dev.off()

    
  }
  

# Layout of plots - 4 plots of equal size
layout(matrix(1:4, 2, byrow = TRUE))  
  data<-read.csv('genomes.50.csv',sep=",", header=FALSE)
  names(data)<-c('type','support','pasta','manual','combined')
  
  data_90<-read.csv('genomes.90.csv',sep=",", header=FALSE)
  names(data_90)<-c('type','support','pasta','manual','combined')
  
  data<-rbind(data,data_90)  
  data<-data[data$type!='H',]
  types<-c()
  for (support in unique(data$support)) {
    for (type in unique(data$type)) {
      row = data[data$type==type & data$support==support,]
      types[[type]] <- draw.pairwise.venn(row$pasta,row$manual,row$combined,c(paste('Pasta',type),paste('Manual',type)),top=type,main=type)
      
    }
    
    
  layout(matrix(1:8))
    pdf(paste("/pedigree2/projects/namphuon/results/hbv/analyses/pdf/venn.genomes.",support,".pdf",sep=""),width=20,height=10)
    grid.arrange(grobTree(types[['A']]),grobTree(types[['B']]),grobTree(types[['C']]),grobTree(types[['D']]),grobTree(types[['E']]),grobTree(types[['F']]),grobTree(types[['G']]),grobTree(types[['I']]),ncol=4)
    dev.off()
    grid.draw(types[['A']])
    frame()
    grid.draw(types[['B']])
    frame()
    grid.draw(types[['C']])
    frame()
    grid.draw(types[['D']])
    frame()
    grid.draw(types[['E']])
    frame()
    grid.draw(types[['F']])
    frame()
    grid.draw(types[['G']])
    dev.off()

    par(4,2)
    p<-types['A']
    p<-types['B']
    p<-types['C']
    p<-types['D']
    p<-types['E']
    p<-types['F']
    p<-types['G']
    print(p)
    dev.off()	  
  }
  draw.pairwise.venn()
}

kmer<-function() {
  data<-read.csv('matches.csv',sep=",")
  df<-data.frame(threshold=c(numeric(0)),found=c(numeric(0)))
  for (t in 1:24) {
    subdata<-data[data$evalue < 10^(-1*t),]
    founds = dim(subdata[subdata$found == 'Found',])[1]
    total = dim(subdata)[1]
    df<-rbind(df,c(10^(-1*t),founds/total))    
  }
  names(df)<-c('threshold','found')
  pdf(paste("/pedigree2/projects/namphuon/temp/hpv/bitscore.pdf",sep=""),width=6,height=4)  
  p <- ggplot(df, aes(x=threshold, y=found))+
    geom_point()+
    geom_line()+
    scale_x_log10()+
    xlab("Min e-value (log scale)")+
    ylab("TP Rate")
  print(p)
  dev.off()	  
    
  

  data<-data.frame(model=c('easy'),paired_reads=c(462800000))
  data$model<-factor(data$model,levels=c('easy','medium','hard'))
  data<-rbind(data,c('medium',462800000))
  data<-rbind(data,c('hard',462800000))
  data$paired_reads<-as.numeric(data$paired_reads)
  
  hg19_size = 3137161264
  hpv16_size = 7906
  integrations = 1000
  hpv16_total_size = hpv16_size*integrations
  
  data$estimated_hpv_reads = data$paired_reads*(hpv16_total_size/hg19_size)
  data$time<-c(55456/3,58173.64/3,61666./3)    
  kmers<-data.frame(kmer_size=c(21),min_count=c(1),time=c(3556),found=c(2904005),all=c(291193),same=c(57461))
  kmers$model<-'easy'
  kmers<-rbind(kmers,c(21,1,3996,2435117,290624,52186,'medium'))  
  kmers<-rbind(kmers,c(25,1,4670,313512,291193,36616,'easy'))  
  kmers$all<-as.numeric(kmers$all)
  kmers$same<-as.numeric(kmers$same)
  kmers$found<-as.numeric(kmers$found)
  kmers$time<-as.numeric(kmers$time)

  
  
  pdf(paste("/pedigree2/projects/namphuon/temp/hpv/positive.pdf",sep=""),width=6,height=4)  
  p <- ggplot(kmers, aes(x=model, y=same/all, color=kmer_size, shape=kmer_size))+
    geom_point()+
    geom_line()+
    xlab("Model")+
    ylab("TP Rate")
    #scale_shape_discrete(name='System')+  
    #scale_color_discrete(name='System')+      
    #facet_wrap(~method)+
    #theme(panel.margin = unit(2, "lines"))    
  print(p)
  dev.off()	  
    
  pdf(paste("/pedigree2/projects/namphuon/temp/hpv/negative.pdf",sep=""),width=6,height=4)  
  p <- ggplot(kmers, aes(x=model, y=1-(same/found), color=kmer_size, shape=kmer_size))+
    geom_point()+
    geom_line()+
    xlab("Model")+
    ylab("FP Rate")
    #scale_shape_discrete(name='System')+  
    #scale_color_discrete(name='System')+      
    #facet_wrap(~method)+
    #theme(panel.margin = unit(2, "lines"))    
  print(p)
  dev.off()	  
  
}

performance <- function() {


  data<-read.csv('/pedigree2/projects/namphuon/results/performance/stats/time.csv',sep=',',header=T)

  pdf(paste("/pedigree2/projects/namphuon/results/performance/pdf/timing.pdf",sep=""),width=6,height=2)  
  p <- ggplot(data, aes(x=cpu, y=time/60., color=system, shape=system))+
    geom_point()+
    geom_line()+
    xlab("CPUs")+
    ylab("Wall clock time (min)")+  
    scale_shape_discrete(name='System')+  
    scale_color_discrete(name='System')+      
    facet_wrap(~method)+
    theme(panel.margin = unit(2, "lines"))    
  print(p)
  dev.off()	  
  
  res<-c()
  comet<-11/(25+7*(1:10))*3*60
  stampede<-11/(25+7*(1:10))*2*60
  res$comet<-comet
  res$stampede<-stampede
  melted<-melt(res)
  names(res)<-
  
}

score_integration <- function() {
  data<-read.csv('/pedigree2/projects/namphuon/results/cancer_viral_integration/HCC/stats/hcc.integration.csv',sep=',',header=T)
data$integration<-paste(data$integration)

rownames(data)<-data$integration

  integrations<-read.csv('/pedigree2/projects/namphuon/results/cancer_viral_integration/HCC/stats/integrations.csv',sep=',',header=F)
  names(integrations)<-c('sample','read','type','orientation','chr','start','end')  
  splits<-integrations[integrations$type=='split',]
  res<-melt(splits[,c('sample','chr','start','end')])
  res$pos<-paste(res$value)
  res$value<-1
  counts<-cast(melt(res))
  counts$pos<-as.numeric(counts$pos)
  counts$support<-0
  for (i in 1:dim(counts)[1]) {
    row = counts[i,]
    support = integrations[integrations$sample %in% row$sample & integrations$chr %in% row$chr & integrations$type=='chim' & integrations$start <= row$pos & integrations$end >= row$pos ,]
    counts[i,]$support = dim(support)[1]    
  }
  
  subdata<-data[data$sample %in% counts$sample,]
  
  miss <- c()
  viral <- c()
  close <- c()
  for (site in data$integration) {
    if (!data[site,]$sample %in% counts$sample) {
      next;
    }
    row <- data[site,]
    found<-counts[counts$sample %in% row$sample & counts$chr %in% row$chr & counts$pos %in% row$pos,]
    if (dim(found)[1] == 0) {
      closest<-counts[counts$sample %in% row$sample & counts$chr %in% row$chr,]
      if (dim(closest)[1] == 0) {
        miss<-rbind(miss,row)            
        next
      }      
      
      closest$pos<-closest$pos-row$pos
      closest<-closest[order(abs(closest$pos)),][1,]
      if (abs(closest$pos) < 25) {
        row$offset = closest$pos
        close<-rbind(close,row)
      } else {
        miss<-rbind(miss,row)            
      }      
    } else {
      viral<-rbind(viral,row)      
      print(site)
    }
  }
  miss$offset<-99999999999
  viral$offset<-0
  
  validated<-read.csv('/pedigree2/projects/namphuon/data/cancer_viral_integration/data/HCC/metadata/validated.csv',sep='\t',header=T)
  validated<-validated[,c('Library','X5..Chr','Exact.hg19.breakpoint.by.sequencing','Validated.successfully.or.not')]
  names(validated)<-c('sample','chr','pos','validated')  
  validated$chr<-paste('chr',validated$chr,sep='')
  validated<-validated[!is.na(validated$pos),]

  true_counts<-counts
  counts<-true_counts[true_counts$support>=2,]
  
  miss <- c()
  viral <- c()
  close <- c()
  for (i in 1:23) {
    row <- validated[i,]
    found<-counts[counts$sample %in% row$sample & counts$chr %in% row$chr & counts$pos %in% row$pos,]
    if (dim(found)[1] == 0) {
      closest<-counts[counts$sample %in% row$sample & counts$chr %in% row$chr,]
      if (dim(closest)[1] == 0) {
        miss<-rbind(miss,row)            
        next
      }      
      
      closest$pos<-closest$pos-row$pos
      closest<-closest[order(abs(closest$pos)),][1,]
      if (abs(closest$pos) < 200) {
        row$offset = closest$pos
        close<-rbind(close,row)
      } else {
        miss<-rbind(miss,row)            
      }      
    } else {
      viral<-rbind(viral,row)      
      print(site)
    }
  }
  
}

viral_detection <- function() {
  
  data<-read.csv('268T.statistics',sep="\t")
  values<-read.csv('268T.hmmsearch.csv',sep=',',header=F)
  name_map<-read.csv('268T.map',sep='\t',header=F)
  names(name_map)<-c('seq','read_name')
  name_map$read_name=gsub("\\s+.*","",name_map$read_name)  
  name_map$read_name=gsub("@","",name_map$read_name)  
  
  names(values)<-c('drop','seq','evalue','bitscore')
  values$drop<-NULL

  merged<-merge(name_map,values)
  merged<-merge(data,merged)  
  
  merged$human_total<-merged$human_read2_match+merged$human_read1_match
  merged$viral_total<-merged$viral_read2_match+merged$viral_read1_match
  rownames(merged)<-paste(merged$read_name)
  
  submerged<-merged[merged$evalue < 1e-4,]
  
  chimeric<-submerged[submerged$human_read1_reference != 'unmapped' | submerged$human_read2_reference != 'unmapped',]
  chimeric$human<-chimeric$human_read1_reference
  chimeric$start<-chimeric$human_read1_start
  chimeric$end<-chimeric$human_read1_end

  chimeric$vstart<-chimeric$viral_read2_start
  chimeric$vend<-chimeric$viral_read2_end
  
  chimeric[chimeric$human == 'unmapped',]$start<-chimeric[chimeric$human == 'unmapped',]$human_read2_start
  chimeric[chimeric$human == 'unmapped',]$end<-chimeric[chimeric$human == 'unmapped',]$human_read2_end  

  chimeric[chimeric$human == 'unmapped',]$vstart<-chimeric[chimeric$human == 'unmapped',]$viral_read1_start
  chimeric[chimeric$human == 'unmapped',]$vend<-chimeric[chimeric$human == 'unmapped',]$viral_read1_end  
  
  chimeric[chimeric$human == 'unmapped',]$human<-chimeric[chimeric$human == 'unmapped',]$human_read2_reference

  subchimeric <- unique(chimeric[,c('human','start','vstart','end','vend')])
  subchimeric<-subchimeric[order(subchimeric$human,subchimeric$start,subchimeric$vstart),]
  subchimeric[grep('ERR092973',rownames(subchimeric)),]
  subchimeric[grep('ERR093207',rownames(subchimeric)),]
  
  
  
pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/HCC/pdf/test.pdf",sep=""),width=20,height=12)  
    p <- ggplot(subchimeric, aes(x=start, y=end))+
      geom_point()+
      xlab("Start")+
      ylab("End")+
      facet_wrap(~human,scales='free')
    print(p);
    dev.off()	  
  
  
pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/HCC/pdf/test.pdf",sep=""),width=6,height=6)  
    p <- ggplot(submerged, aes(x=log(evalue), y=human_total))+
      geom_point()+
      xlab("E-value")+
      ylab("Human total")
    print(p);
    dev.off()	  
    
  pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/HCC/pdf/test.pdf",sep=""),width=6,height=6)  
    p <- ggplot(merged, aes(x=human_total, y=viral_total,size=bitscore))+
      geom_point()+
      xlab("Match human")+
      ylab("Match virus")
    print(p);
    dev.off()	  

pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/HCC/pdf/test.pdf",sep=""),width=6,height=6)  
    p <- ggplot(submerged, aes(x=log(evalue), y=human_total))+
      geom_point()+
      xlab("E-value")+
      ylab("Human total")
    print(p);
    dev.off()	    
}

coronas <- function() {
  data<-read.csv("/pedigree2/projects/namphuon/data/corona/data/lengths.csv",sep="\t",header=F)
  names(data)<-c('genome','type','length')  
}

plot_distribution <- function() {
  data<-read.csv("/pedigree2/projects/namphuon/results/cancer_viral_integration/integration_sites/analyses/upp/empirical.csv",sep=",",header=F)
  names(data)<-c("read","sample","type","idx","len","char")    
  
  pdf(paste("/pedigree2/projects/namphuon/data/hbv/pdf/distribution_percent.pdf",sep=""),width=6,height=6)  
    p <- ggplot(subdata, aes(type, char/len))+
      facet_wrap(~sample)+
      xlab("Type")+
      ylab("Counts")+
      geom_boxplot()    
    print(p);
    dev.off()	  

  pdf(paste("/pedigree2/projects/namphuon/data/hbv/pdf/distribution_char.pdf",sep=""),width=6,height=6)  
    p <- ggplot(subdata, aes(type, char))+
      facet_wrap(~sample)+
      xlab("Type")+
      ylab("Counts")+
      geom_boxplot()    
    print(p);
    dev.off()	  
    
  pdf(paste("/pedigree2/projects/namphuon/data/hbv/pdf/distribution_length.pdf",sep=""),width=6,height=6)  
    p <- ggplot(subdata, aes(type, len))+
      facet_wrap(~sample)+
      xlab("Type")+
      ylab("Counts")+
      geom_boxplot()    
    print(p);
    dev.off()	  
  
  subdata<-data[data$char/data$len > 0.80 & data$char > 80  ,]  
  write.csv(subdata,file='/pedigree2/projects/namphuon/results/cancer_viral_integration/integration_sites/analyses/upp/reduced.csv',quote = F,row.names = F)
}

plot_hbv_support <- function() {
  data<-read.csv("/pedigree2/projects/namphuon/data/hbv/trees/log.csv",sep=",",header=F)
  names(data)<-c('bipartition','upp_pasta','upp_manual')    
  melted<-melt(data)
pdf(paste("/pedigree2/projects/namphuon/data/hbv/pdf/support.pdf",sep=""),width=6,height=6)  
  p <- ggplot(melted, aes(variable, value))+
    xlab("Alignment")+
    ylab("Support")+
    geom_boxplot()    
  print(p);
  dev.off()	  
  
}

plot_chromo <- function() {
  #repeats<-read.csv("/pedigree2/projects/namphuon/results/cancer_viral_integration/HPV/repeat_masker/repeats.csv",sep=" ",header=T);
  #save(repeats,file='/pedigree2/projects/namphuon/results/cancer_viral_integration/HPV/repeat_masker/repeats.dat')
  #load('/pedigree2/projects/namphuon/results/cancer_viral_integration/HPV/repeat_masker/repeats.dat')
  #dir = '/pedigree2/projects/namphuon/results/cancer_viral_integration/originalTCGA'
  dir = '/pedigree2/projects/namphuon/results/cancer_viral_integration/clinical'
  setwd(dir)
  dirs <- list.dirs()
  #dirs<-gsub("./","",dirs[grep('TCGA',dirs)])
  dirs<-gsub("./","",dirs[dirs != "."])
  scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
  }
  threshold = '10'
  lens = '50000'
  for (d in dirs) {
    virus<-read.csv(paste(dir,"/",d,'/',threshold,'.hpv.',lens,'.viral.csv',sep=''),sep=",",header=F);
    names(virus)<-c('read','chr','region','start','end')
    virus$mid<-(virus$start+virus$end)/2
    counts<-cast(melt(virus[,c('region','start')]))
    names(counts)<-c('region','supporting_reads')
    if (length(unique(virus$region)) == dim(virus)[1]) {
      counts$supporting_reads=1
    }

    #data<-read.csv(paste(dir,"/",d,'/',threshold,'.hpv.',lens,'.csv.fixed',sep=''),sep=",",header=F);
    data<-read.csv(paste(dir,"/",d,'/',threshold,'.hpv.',lens,'.csv',sep=''),sep=",",header=F);
    names(data)<-c('region','pos','type','coverage')  
    data<-merge(data,counts)
    virus<-merge(virus,counts)
    data<-data[grep("GL", data$region,invert=TRUE),]
    virus<-virus[grep("GL", virus$region,invert=TRUE),]
    size = length(unique(virus$region))
    if (size == 0) {
      next
    }
    maxRows <- by(virus, virus$region, function(X) X[which.max(X$mid),])
    res<-do.call("rbind", maxRows)
    
    data$region_name<-paste(data$region,';',data$supporting_reads,' chimeric reads',sep='')
    virus$region_name<-paste(virus$region,';',virus$supporting_reads,' chimeric reads',sep='') 
    data<-data[data$coverage != 0,]
    
    means<-cast(melt(data[,c('type','coverage')]),fun=mean)
#    for (region in unique(means$region)) {
      for (t in unique(means$type)) {
        #data[data$type == t,]$coverage<-data[data$type == t,]$coverage/means[means$type == t,]$coverage      
        #data[data$region == region & data$type == t,]$coverage<-data[data$region == region & data$type == t,]$coverage/means[means$region == region & means$type == t,]$coverage
      }
#    }
             pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/integration_sites/pdf/masked.",d,".",lens,".coverage.pdf",sep=""),width=round(sqrt(size)-1)*2+4,height=round((sqrt(size)-1)*2)+4)
   
    p<-ggplot(data, aes(x = pos,y=coverage,color=type,line=type))+
      #geom_point(size=1)+
      #geom_path(size = 0.5)+      
      geom_smooth()+
      #geom_density(kernel = "rectangular")+
      facet_wrap(~region_name,scales="free_x")+
      xlab('Position')+      
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      #scale_x_continuous(label=scientific_10)+
      expand_limits(y=0)+
      geom_vline(aes(xintercept = mid,color='integration',linetype='longdash'),virus)+
      ylab('Coverage')    
    print(p);
    dev.off()	  
  }
}

plot_virus_sites <-function() {
    data<-read.csv("/pedigree2/projects/namphuon/results/cancer_viral_integration/integration_sites/stats/integration.csv",sep="\t",header=F);
  names(data)<-c('sample','type','read','chromosome','start','end')
  data$uid<-paste(data$sample,data$type)
  
  
  res<-array(0,dim=c(length(unique(data$chromosome)),length(unique(data$uid)),9000))
  dimnames(res)[1]<-list(unique(paste(data$chromosome)))
  dimnames(res)[2]<-list(unique(paste(data$uid)))
  
  for (i in 1:dim(data)[1])  {
#"chr3","C5-A0TN TP"
  #data[data$chromosome=='chr3' & data$uid=='C5-A0TN TP',]
  #i=4    
  #sum(res['chr1',"EK-A2RE TP",])
    res[paste(data[i,]$chromosome),paste(data[i,]$uid),data[i,]$start:data[i,]$end]<-res[paste(data[i,]$chromosome),paste(data[i,]$uid),data[i,]$start:data[i,]$end]+1      
  }
  temp = c()
  for (i in 1:dim(res)[1]) {
    melted<-melt(res[i,,])  
    names(melted)<-c('sample','position',paste(dimnames(res)[[1]][i]))
    if (is.null(dim(temp))) {
      temp<-melted
    } else {
      temp[,paste(dimnames(res)[[1]][i])]<-melted[,paste(dimnames(res)[[1]][i])]
    }
  }
  melted<-melt(temp,id=c('sample','position'))
  #ares<-acast(melted[,c('sample','position','value')],sample~position,fun=sum)  
  #amelted<-melt(ares)
  #names(amelted)<-c('sample','position','value')  
  #amelted<-amelted[amelted$value !=0,]
pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/integration_sites/pdf/density_graph.pdf",sep=""),width=12,height=12)
  p<-ggplot(melted, aes(x = position,weight=value,fill=variable))+
    geom_density(kernel = "rectangular",n=1)+
    facet_wrap(~sample,scales="free_y")+
    xlab('Position')+
    xlim(0,9000)+
    ylab('Coverage')    
  print(p);
  dev.off()	
  
  
pdf(paste("/pedigree2/projects/namphuon/results/cancer_viral_integration/integration_sites/pdf/line_graph.pdf",sep=""),width=12,height=12)
  p<-ggplot(amelted, aes(x = position,y=value))+
    geom_point(size=1)+
    #geom_path(size = 0.5)+
    #geom_density(kernel = "rectangular")+
    facet_wrap(~sample,scales="free_y")+
    xlab('Position')+
    xlim(0,9000)+
    expand_limits(y=0)+
    ylab('Coverage')    
  print(p);
  dev.off()	
  
}

virus_mapping <-function() {
  library(rbamtools)
  bam <- "clean.bam"
  idx <- "clean.bam.bai"
  reader<-bamReader(bam)
  loadIndex(reader,idx)
  indexInitialized(reader)

  align<-getNextAlign(reader)  
  
}

virus_graph <-function() {
  library(igraph)
  data<-read.csv("/pedigree/projects/namphuon/results/cancer_viral_integration/HPV/blast_all/blast.csv",sep="\t",header=F);
  names(data)<-c('query','target','identity','length','mismatch','gap_open','start_query','end_query','start_target','end_target','e_value','bit_score')
  res<-acast(data[,c('query','target','e_value')],query~target,fun=mean)  
  res[is.na(res)] <- 100  
  
  smooth<-res
  smooth[res < 1e-5] = 1
  smooth[res >= 1e-5] = 0
  keep<-colSums(smooth)  
  keep=keep[keep > 10]
    
  clusters=graph.adjacency(res[names(keep),names(keep)],weighted=TRUE)
  
  layout1 <- layout.fruchterman.reingold(clusters)  
  layout1 <- layout.auto(clusters)  
  
    pdf(paste("/pedigree/projects/namphuon/results/cancer_viral_integration/HPV/blast_all/graph.pdf",sep=""),width=12,height=12)
        p<-plot(clusters, layout=layout1,edge.arrow.mode='--')
        print(p);
        dev.off()	
  
}

insects <-function() {
  data<-read.csv("/pedigree2/projects/namphuon/data/julie_insects_1/rates.csv",sep=",",header=F);
  names(data)<-c('gene','codon','alpha','ac','ag','at','cg','ct','gt')
  data<-data[data$codon != 3,]
  data$name<-paste(data$gene,data$codon)
  data$codon<-factor(data$codon,levels=c(1,2,3))
  pcaf<-princomp(scale(subdata[,c('alpha','ac','ag','at','cg','ct')]))
  pca<-data.frame(pcaf$scores[,c(1,2)])    
  pca$codon<-subdata$codon
  
  subdata<-data
  bad_names<-c()
  for (i in 0:10) {
    for (name in c('alpha','ac','ag','at','cg','ct')) {
      subdata$outlier <- (subdata[,name]-mean(subdata[,name]))/sd(subdata[,name])
      bad_names<-union(bad_names,subdata[abs(subdata$outlier)>10,]$name)
    }
    subdata<-subdata[which(!(subdata$name %in% bad_names)),]
  }
  
  subdata<-data[data$name %in% bad_names,]
  subdata<-data
  subdata$gene<-NULL
  subdata$codon<-NULL
  melted<-melt(subdata)  
  for (name in c('alpha','ac','ag','at','cg','ct')) {
    melted[melted$variable==name,]$value <- (melted[melted$variable==name,]$value-mean(melted[melted$variable==name,]$value))/sd(melted[melted$variable==name,]$value)
  }
  melted<-melted[melted$variable != 'gt',]
  
pdf(paste("/pedigree2/projects/namphuon/data/julie_insects_1/rates.sd.pdf",sep=""),width=6,height=6)
      p<-ggplot(melted,aes(variable,value))+
        ylab("Sd from mean")+
        xlab("GTR Parameter")+
        geom_boxplot()+
        scale_fill_discrete(name='GTR Parameter')
      print(p);
      dev.off()
  
  
  pcaf<-princomp(scale(subdata[,c('alpha','ac','ag','at','cg','ct')]))
  pca<-data.frame(pcaf$scores[,c(1,2)])    
  pca$codon<-subdata$codon
  
  
  pdf(paste("/pedigree2/projects/namphuon/data/julie_insects_1/insect_16_taxa.pca.codon.pdf",sep=""),width=6,height=6)
      p<-ggplot(pca,aes(x=Comp.1,y=Comp.2,color=pca$codon))+
        xlab("PCA 1")+
        ylab("PCA 2")+
        scale_color_discrete(name='Codon position')+      
        geom_point(size=1)
      print(p);
      dev.off()
    
  kres <- c()
  for (i in 2:30) {
    res<-kmeans(pca[,c("Comp.1","Comp.2")],i,nstart=25)
    kres<-rbind(kres,c(i,res$betweenss/res$totss,sum(res$withinss)))
  }
  kres<-data.frame(kres)
  names(kres)<-c('clusters','betweenss_totalss_ratio','within_ss')

  res<-kmeans(pca[,c("Comp.1","Comp.2")],10,nstart=25)
  df<-data.frame(subdata$gene,subdata$codon,res$cluster)
  names(df)<-c('gene','codon','class')
  write.csv(df,file='/pedigree2/projects/namphuon/data/julie_insects_1/class.csv',quote = F,row.names = F)

}

se <- function(x) {
  res<-sd(x) /  sqrt(length(x)) 
  return(res)
}

gbm()