# load required packages

library(TCGAbiolinks)
library(tidyr)
library(data.table)
library(SummarizedExperiment)
library(mclust)

# TCGA projects

projects <- data.table(getGDCprojects())
projects[grep("TCGA",id),tumor] -> tumor.types

tumor.types <- "COAD" # for testing purposes

genes.of.interest <- c("TP53","KRAS","APC","BRAF","PIK3CA","SMAD4")


for(tumor in tumor.types) {
  
  load(sprintf("output/%s_data.rda",tumor),verbose=T)
  
  # clinical includes several tables
  
  x <- grep("clinical_patient",names(clinical))
  
  patients <- as.data.table(clinical[[x]])[grep("TCGA",bcr_patient_barcode)]
  mutations <- as.data.table(mutations)
  mutations[,bcr_patient_barcode:=substring(Tumor_Sample_Barcode,1,12)]
  mutations[,HGVSp_Short:=gsub("^p.","",HGVSp_Short)]
  
  patients <- merge(patients,mutations[,list(Nmut=.N),by=bcr_patient_barcode],all.x=T)
  
  # identify if there are clear groups of patiens according to the number of mutations (i.e. a mutator phenotype)
  
  foo <- patients[!is.na(Nmut),list(bcr_patient_barcode,Nmut)]
  foo[,Mut_Cluster := paste0("Cluster_",Mclust(G=1:4,modelNames="E",log2(Nmut+.1))$classification)]
  foo[,N_in_Cluster := .N, by=Mut_Cluster]
  plot1 <- ggplot(foo) + 
    geom_histogram(aes(log10(Nmut+1),fill=Mut_Cluster),bins=round(nrow(foo)/10)) + 
    ggtitle(paste("Mutations in",tumor))
  plot1
  patients <- merge(patients,foo[,list(bcr_patient_barcode,Mut_Cluster,N_in_Cluster)],all.x=T)
  
  # 50 most frequently mutated genes
  
  mutations <- merge(mutations,unique(patients[,list(bcr_patient_barcode,Mut_Cluster,N_in_Cluster)]),all.x=T)
  
  foo <- mutations[Variant_Classification!="Silent",list(N=.N),by=list(Hugo_Symbol,Mut_Cluster,N_in_Cluster)]
  foo[,mut_freq:=N/N_in_Cluster]
  foo[order(mut_freq,decreasing=T)]
  
  x <- dcast(mutations[Hugo_Symbol %in% genes.of.interest],bcr_patient_barcode ~ Hugo_Symbol,value.var = "HGVSp_Short",fun.aggregate = function(x) paste(unique(x) %>% sort,collapse="|"))
  x[x==""] <- "WT"
  patients <- merge(patients,x,all.x=T)
  
  
  # PCA of the methytion
  
  meth <- assays(methylation)[[1]]
  
  
}

