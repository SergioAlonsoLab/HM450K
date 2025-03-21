# load required packages

#library(TCGAbiolinks)
#library(SummarizedExperiment)

library(tidyr)
library(data.table)
library(mclust)
library(ggplot2)

# Prepared TCGA projects

dir("output/",pattern = "*_data") %>% gsub("_data.rda","",.) -> tumor.types

tumor.types <- "COAD" # for testing purposes

if(!file.exists("misc/hm450k_manifest.csv")) {
  curl::curl_download("https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv",
                      "misc/hm450k_manifest.csv")
  grep("(^cg|^IlmnID)",readLines("misc/hm450k_manifest.csv"),value=T) %>% writeLines("misc/hm450k_manifest.csv")
}

genes.of.interest <- fread("misc/true_sight_genes.txt",header = F)
names(genes.of.interest) <- "hgnc_symbol"

for(tumor in tumor.types) {
  
  load(sprintf("output/%s_data.rda",tumor),verbose=T)
  
  # clinical includes several tables
  
  patients <- as.data.table(clinical[[grep("clinical_patient",names(clinical))]])[grep("TCGA",bcr_patient_barcode)]
  mutations <- as.data.table(mutations)
  mutations[,bcr_patient_barcode:=substring(Tumor_Sample_Barcode,1,12)]
  mutations[,HGVSp_Short:=gsub("^p.","",HGVSp_Short)]
  
  patients <- merge(patients,mutations[,list(Nmut=.N),by=bcr_patient_barcode],all.x=T)
  
  # identify if there are clear groups of patients according to the number of mutations (i.e. a mutator phenotype)
  
  foo <- patients[!is.na(Nmut),list(bcr_patient_barcode,Nmut)]
  foo[,Mut_Cluster := paste0("Cluster_",Mclust(G=1:4,modelNames="E",log2(Nmut+.1))$classification)]
  foo[,N_in_Cluster := .N, by=Mut_Cluster]
  plot.mutations.1 <- ggplot(foo) + 
    geom_histogram(aes(log10(Nmut+1),fill=Mut_Cluster),bins=round(nrow(foo)/10)) + 
    ggtitle(paste("Mutations in",tumor))
  plot.mutations.1
  patients <- merge(patients,foo[,list(bcr_patient_barcode,Mut_Cluster,N_in_Cluster)],all.x=T)
  
  # non silent mutations per gene, considering only the genes in TruSight
  foo <- dcast(mutations[Hugo_Symbol %in% genes.of.interest$hgnc_symbol],bcr_patient_barcode ~ Hugo_Symbol,value.var = "HGVSp_Short",fun.aggregate = function(x) paste(unique(x) %>% sort,collapse="|"))
  foo[foo==""] <- "WT"
  
  # remove genes mutated in less than 5 samples
  foo <- foo[,.SD,.SDcols = which(colSums(foo!="WT") >= 5)]
  names(foo)[-1] <- paste0("Gene:",names(foo)[-1])
  patients <- merge(patients,foo,all.x=T)
  
  # add CNA information 
  
  cna <- as.data.table(cna)
  cna[,bcr_patient_barcode := substr(Sample,1,12)]
  cna[,width := End - Start]
  
  foo <- cna[,list(Ploidy=weighted.mean(Copy_Number,width)),by=bcr_patient_barcode]
  foo[,Ploidy_Cluster := Mclust(Ploidy,G = 4,modelNames = "E")$classification %>% paste0("Cluster_",.)]
  
  patients <- merge(patients,foo,all.x=T)
  

  saveRDS(patients,file=sprintf("output/%s_patient_table.rds",tumor),compress = T)
  
  
}

