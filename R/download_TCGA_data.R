# Download the data from TCGA

# install necessary libraries if necessary

ins_packs <- installed.packages()[,"Package"]
if(!"BiocManager" %in% ins_packs) install.packages("BiocManager")
if(!"tidyr" %in% ins_packs) install.packages("tidyr")
if(!"TCGAbiolinks" %in% ins_packs) BiocManager::install("TCGAbiolinks")

# load required packages

library(TCGAbiolinks)
library(tidyr)
library(data.table)


# TCGA projects

projects <- data.table(getGDCprojects())
projects[grep("TCGA",id),tumor] -> tumor.types

# download clinical, mutational, CNA, and methylation data

for(tumor in tumor.types) {
  
    project <- paste0("TCGA-",tumor)
    
    methylation_query <- GDCquery(project = sprintf("TCGA-%s",tumor),
                                  data.category = "DNA Methylation",
                                  platform = "Illumina Human Methylation 450",
                                  data.type = "Methylation Beta Value")
    GDCdownload(methylation_query,files.per.chunk = 15,directory = "GDCdata/") 
    
    
    clinical_query <- GDCquery(project = project, 
                               data.category = "Clinical",
                               data.type = "Clinical Supplement",
                               data.format = "BCR Biotab")
    GDCdownload(clinical_query,files.per.chunk = 30,directory = "GDCdata/") 
    
    mutations_query <- GDCquery(project = project,
                                data.category = "Simple Nucleotide Variation",
                                data.type = "Masked Somatic Mutation")
    GDCdownload(mutations_query,files.per.chunk = 15,directory = "GDCdata/") 
    
    cna_query <- GDCquery(project = project,
                          data.category="Copy Number Variation",
                          data.type="Allele-specific Copy Number Segment",
                          workflow.type = "ASCAT3")
    GDCdownload(cna_query,files.per.chunk = 15,directory = "GDCdata/")
    
    
    save(clinical_query,methylation_query,mutations_query,cna_query,file = sprintf("output/%s_queries.rda",tumor))
    
}



length(tumor.types)


