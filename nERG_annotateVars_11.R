### Annotate nERG variants
### prepares VCF file of variants for downstream analysis 
### takes genetic testing results, parses them, updates their position to latest reference transcript, maps to hg38
### Provided as is
### version 1.1
### Kirill Zaslavsky


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")


#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
#BiocManager::install("biomaRt")



#library("myvariant")
##clinVar

library(rentrez)
library(XML)
library(seqinr)
library(Biostrings)
library(biomaRt)
library(GenomicFeatures)
library(ensembldb)
#library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(AnnotationHub)
library(data.table)
library(tidyverse)
source("nERG_function.R")
source("nERG_genomics_function.R")
#source("nERG_sortGeneVars_v0.4.R")

library(httr)
library(jsonlite)
library(xml2)


### Prepare
#####

## server for ensembl REST API
server <- "https://rest.ensembl.org"

## To download latest ens db
## Load the annotation resource. 
ah <- AnnotationHub()

## Query for all available EnsDb databases
#query(ah, "EnsDb")
ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 111))
ahEdb <- ahDb[[1]]

#load bm ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast") #uncomment if loaded


## bruteforce if above doesn't connect 
ensembl <- getMart("genes","hsapiens_gene_ensembl")

# init dump file
mdf2.dump <- data_frame(Variant = character())
 

# split variant into gene, cDNA variant, ptn variant
melted_df[, c("Gene_Name", "cDNA_Variant", "PTN_Variant", "AccNum", "Type") := as.list(splitGeneVariant(Gene_Variant)), by = Gene_Variant] 


# load gene variants to parse
genevars <- melted_df[, Gene_Variant]


######
# create VCF file for all the variants

output_vcf_path <- "vcf_tbl_apr_16_2024.txt" #keep as txt without header _ construct line by line
saved_vcf_file <- fread("vcf_tbl_apr_16_2024.txt", header = TRUE)



result <- CreateVCF(genevars, ensembl, melted_df, mdf2_dump = mdf2.dump, output_file, old_table = saved_vcf_file) ### run if connection to NCBI is interrupted / not starting from scratch
result <- CreateVCF(genevars, ensembl, melted_df, mdf2_dump = mdf2.dump, output_file, old_table = NULL) ### run if starting from scratch

# Use this function to run createVCF even if it crashes due to connection errors
final_result <- run_create_vcf(genevars, ensembl, melted_df, mdf2_dump, output_file)











