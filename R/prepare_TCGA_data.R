# load required packages

library(TCGAbiolinks)
library(tidyr)
library(data.table)

# TCGA projects

projects <- data.table(getGDCprojects())
projects[grep("TCGA",id),tumor] -> tumor.types

# Load the queries generated in download_TCGA_data.R

for(tumor in tumor.types) {
  
  load(file = sprintf("output/%s_queries.rda",tumor),verbose=T)
  
  clinical <- GDCprepare(clinical_query)
  #mutations <- GDCprepare(mutations_query)
  #methylation <- GDCprepare(methylation_query)
  #cna <- GDCprepare(cna_query)
  
  save(clinical,mutations,methyalation,cna,file=sprintf("output/%s_data.rda",tumor))
  
}