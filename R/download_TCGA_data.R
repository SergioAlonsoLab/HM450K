# Download the data from TCGA

# install necessary libraries if necessary

ins_packs <- installed.packages()[,"Package"]
if(!"ggplot2" %in% ins_packs) install.packages("ggplot2")
if(!"BiocManager" %in% ins_packs) install.packages("BiocManager")
if(!"tidyr" %in% ins_packs) install.packages("tidyr")
if(!"TCGAbiolinks" %in% ins_packs) BiocManager::install("TCGAbiolinks")

# load required packages

library(ggplot2)
library(TCGAbiolinks)
library(tidyr)
library(data.table)

# Store the data in environments by tumor type

tumor.types <- c("LUAD","")

for(tumor in tumor.types) {
  assign(tumor,new.env()) # create the environments
}

# download methylation data

for(tumor in tumor.types) {
  
  with(get(tumor),{
    
    methylation_query <- GDCquery(project = sprintf("TCGA-%s",tumor), 
                                  data.category = "DNA Methylation",
                                  platform = "Illumina Human Methylation 450",
                                  data.type = "Methylation Beta Value")
    GDCdownload(methylation_query,files.per.chunk = 15,directory = "GDCdata/") 

  })
  
}

# NOTE from Sergio: this strategy will create very large objects eating up the memory 
# perhaps a good option would be to download the data, save it to RDS objects 
# in the large_files_folder, and then empty the environment



