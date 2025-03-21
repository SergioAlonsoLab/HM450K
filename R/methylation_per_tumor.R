# load required packages

library(tidyr)
library(data.table)
library(SummarizedExperiment)
library(ggplot2)

# library(eimpute) # to impute missing values in the methylation matrices

# TCGA projects

dir("output/",pattern = "methylation") %>% gsub("_.+","",.) -> tumor.types

# tumor.types <- "COAD" # for testing purposes

# Download HM450K manifest file

if(!file.exists("misc/hm450k_manifest.csv")) {
  curl::curl_download("https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv",
                      "misc/hm450k_manifest.csv")
  grep("(^cg|^IlmnID)",readLines("misc/hm450k_manifest.csv"),value=T) %>% writeLines("misc/hm450k_manifest.csv")
}

hm450.manifest <- fread("misc/hm450k_manifest.csv")
hm450.manifest[Relation_to_UCSC_CpG_Island=="Island",CGI:="Island"]
hm450.manifest[grep("Shore",Relation_to_UCSC_CpG_Island),CGI:="Shore"]
hm450.manifest[grep("Shelf",Relation_to_UCSC_CpG_Island),CGI:="Shelf"]
hm450.manifest[Relation_to_UCSC_CpG_Island=="",CGI:="Open Sea"]

# iterate for all tumor types

for(tumor in tumor.types) {
  
  message("Now analysing ",tumor)
  
  methylation <- readRDS(sprintf("output/%s_methylation_data.rds",tumor))
  
  meth <- assays(methylation)[[1]]

  # select only tumors
  meth <- meth[,substring(colnames(meth),14,15) == "01"] 
  
  # remove sexual chromosomes (they have a BIG effect the PCAs)
  
  probes <- hm450.manifest[!CHR %in% c("X","Y"),IlmnID]
  probes <- intersect(probes,rownames(meth))
  probes <- grep("^cg",probes,value=T)
  
  meth <- meth[probes,]

  # select rows with valid values in at least 90% of the samples (exclude unreliable or poorly evaluated probes)
  
  meth <- meth[(rowSums(!is.na(meth)) / ncol(meth) > 0.9),]
  
  # PCA will fail if NAs are included Substitute NAs by imputing values 
  # I tested eimputed by introducing NAs into a methylation matrix, and it
  # does a extremely accurate imputation. But it takes too much long time.
  # perhaps in the cluster?
  
  nas <- rowSums(is.na(meth))
  
  # substitute by the mean
  
  for(i in which(nas > 0)) {
    
    x <- meth[i,]
    message("row ",i," had ",sum(is.na(x))," NAs. Substituted by ",mean(x,na.rm=T))
    x[is.na(x)] <- mean(x,na.rm=T)

    meth[i,] <- x
    
  }
  
  # average methylation in islands, shores, shelves, and open sea regions
  
  colnames(meth) <- substr(colnames(meth),1,12) # just the patient barcode
  
  M.analsysis <- data.table(bcr_patient_barcode=colnames(meth),
                            M.Mean=colMeans(meth),
                            M.Island=intersect(hm450.manifest[CGI=="Island",Name],rownames(meth)) %>% meth[.,] %>% colMeans,
                            M.Shore=intersect(hm450.manifest[CGI=="Shore",Name],rownames(meth)) %>% meth[.,] %>% colMeans,
                            M.Shelf=intersect(hm450.manifest[CGI=="Shelf",Name],rownames(meth)) %>% meth[.,] %>% colMeans,
                            M.OpenSea=intersect(hm450.manifest[CGI=="Open Sea",Name],rownames(meth)) %>% meth[.,] %>% colMeans)
  

  saveRDS(M.analsysis,file="output/%s_M_Analysis.rds")
  
  # PCA (warning: prcomp is computationally intensive)
  
  # to reduce computational load, we can use irlba
  
  gc()
  try({
    pca1 <- prcomp(t(meth),scale. = F)
    saveRDS(pca1,file=sprintf("output/%s_PCA.rds",tumor))})
  
}

