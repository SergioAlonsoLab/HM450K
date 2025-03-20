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
  
  
  
  
  
  
  # 50 most frequently mutated genes
  
  mutations[Variant_Classification!="Silent",list(N=.N),by=Hugo_Symbol][order(N,decreasing=T)][1:50]
  
  
  
  x <- dcast(mutations[Hugo_Symbol %in% genes.of.interest],bcr_patient_barcode ~ Hugo_Symbol,value.var = "HGVSp_Short",fun.aggregate = function(x) paste(unique(x) %>% sort,collapse="|"))
  x[x==""] <- "WT"
  patients <- merge(patients,x,all.x=T)
  
  
  # PCA of the methytion
  
  meth <- assays(methylation)[[1]]
  
  
}

