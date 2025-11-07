#nERG project - script to categorize variants and interpret cases
#Author Kirill Zaslavsky
#Version 0.1
#provided as is

#merge the initial db.nERG table with the variants, increasing the number of rows as necessary

setDT(db.nERG)[, ID := as.character(ID)]

# Merge db.nERG with vcf_parsed_data_with_ann_filtered_unique
merged_dbNERG_geneVar <- full_join(vcf_parsed_withLargeSVs, db.nERG, 
                                    by = c("INFO_PtID" = "ID"))

merged_dbNERG_geneVar_unique <- merged_dbNERG_geneVar %>%
  distinct(INFO_PtID, CHROM, POS, .keep_all = TRUE)


merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_unique %>%
  dplyr::select(Presumed_Genetic, INFO_PtID, Gene_ID, Feature_ID, HGVS.c, HGVS.p, Annotation, Annotation_Impact, INFO_dbNSFP_clinvar_clnsig, 
         INFO_dbNSFP_clinvar_review, INFO_dbNSFP_REVEL_score, INFO_dbNSFP_AlphaMissense_score, INFO_dbNSFP_MutPred_score, INFO_dbNSFP_VEST4_score,
         SpliceAI_DS_AG, SpliceAI_DS_AL, SpliceAI_DS_DG, SpliceAI_DS_DL,  INFO_dbNSFP_gnomAD_exomes_AF, INFO_dbNSFP_gnomAD_genomes_AF, 
         Pedigree_Inheritance, Disease_cause_AR_AD_XL, `Trans_proven?`, Phenotype, everything())


merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>% add_column(PVS1 = as.character(NA),  #PTC, fs, splice within 1-2 bp
                                               PS1 = as.character(NA),   #same AA change as another path variant
                                               PS2 = as.character(NA),   #de novo (both maternity and paternity confirmed)
                                               PS3 = as.character(NA),   #in vitro / in vivo exp evidence
                                               PS4 = as.character(NA),   #case-control studies OR > 5.0 or no presence in controls if rare
                                               PM1 = as.character(NA),   #mutation hotspot, functional domain without variation
                                               PM2 = as.character(NA),   #absent from controls in pop dbs (i.e. AF < 0.01 AR or 0 if AD) - supporting only?
                                               PM3 = as.character(NA),   #in trans with a pathogenic variant
                                               PM4 = as.character(NA),   #protein length changes due to in-frame deletions
                                               PM5 = as.character(NA),   #novel missense where another change was prev pathogenic
                                               PM6 = as.character(NA),   #assumed de novo
                                               PP1 = as.character(NA),   #co-segregation in pedigree
                                               PP2 = as.character(NA),   #missense variant in gene that has low rate of missense variation and missense vars cause disease
                                               PP3 = as.character(NA),   #in silico (i.e. REVEL thresholds as in Pejaver)
                                               PP4 = as.character(NA),   #phenotype or FHx specific for disease
                                               PP5 = as.character(NA),   #reputable source
                                               .before = "Phenotype"
)

#add final characterization column
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>% add_column(Var_Interpretation = as.character(NA), .before = "PVS1")

# If you need to apply this to multiple columns, you can do it like this:
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(INFO_dbNSFP_gnomAD_exomes_AF = as.numeric(INFO_dbNSFP_gnomAD_exomes_AF),
         INFO_dbNSFP_gnomAD_genomes_AF = as.numeric(INFO_dbNSFP_gnomAD_genomes_AF)) %>%
  mutate(across(c(INFO_dbNSFP_gnomAD_exomes_AF, INFO_dbNSFP_gnomAD_genomes_AF), ~replace_na(.x, 0)))

#check PP3 - REVEL

#process REVEL column to extract number
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(INFO_dbNSFP_REVEL_score = str_extract(INFO_dbNSFP_REVEL_score, "[^,]+") %>% 
           as.numeric())

#process AlphaMissense column to extract number
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(INFO_dbNSFP_AlphaMissense_score = str_extract(INFO_dbNSFP_AlphaMissense_score, "[^,]+") %>% 
           as.numeric())

#process MutPred column to extract number
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(INFO_dbNSFP_MutPred_score = str_extract(INFO_dbNSFP_MutPred_score, "[^,]+") %>% 
           as.numeric())

#process VEST4 column to extract number
merged_dbNERG_geneVar_reordered <- as.data.table(merged_dbNERG_geneVar_reordered %>%
  mutate(INFO_dbNSFP_VEST4_score = str_extract(INFO_dbNSFP_VEST4_score, "[^,]+") %>% 
           as.numeric()))


#####
## optional - merge with old interpretation dataset

## merge old interpretations into the new table
old_interp <- fread("nERG_Apr19_withClinVarPM3PP1_interp.txt", header = TRUE)
slice_old_interp <- old_interp[,c(2,27:48)]
cols_to_update <- names(slice_old_interp)[-c(1,22,23)]
setDT(slice_old_interp)[, INFO_PtID := as.character(INFO_PtID)]
setDT(slice_old_interp)[, POS := as.character(POS)]

#add remaining columns

merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>% add_column(clinvar_pm3_transseg = as.character(NA),  
                                                                                  clinvar_pp1_diseaseseg = as.character(NA),   
                                                                                  clinvar_comments = as.character(NA),   
                                                                                  .before = "Phenotype"
)

# Perform a left join to merge old_data into new_data based on keys
merged_data <- merge(slice_old_interp, merged_dbNERG_geneVar_reordered, by = c("INFO_PtID", "CHROM", "POS"), all.y = TRUE, suffixes = c("", "_old"))

# Update columns in new_data with values from old_data if they are NA
#cols_to_update <- c("PVS1", "PS1", "PS2")

for (col in cols_to_update) {
  merged_data[, (col) := ifelse(is.na(get(col)), get(paste0(col, "_old")), get(col))]
}

# Drop the old data columns
merged_data[, (paste0(cols_to_update, "_old")) := NULL]

#####

output_CVtable_PM3_PP1_path <- "CV_PM3_PM1.txt"
#write.table(merged_dbNERG_clinvar_test, file = output_CVtable_PM3_PP1_path, row.names = FALSE, sep = "\t")

#AnnoCV_PM3_PM1_file <- fread(output_CVtable_PM3_PP1_path, header = TRUE)

#merged_dbNERG_geneVar_reordered <- Annotate_CV_PM3_PP1(merged_dbNERG_geneVar_reordered, existing_tbl = AnnoCV_PM3_PM1_file)
#merged_dbNERG_geneVar_reordered <- as.data.table(merged_dbNERG_geneVar_reordered)
#merged_dbNERG_geneVar_reordered <- cbind(  merged_dbNERG_geneVar_reordered, AnnoCV_PM3_PM1_file[,-(1:3)])


# go to analysis
#####


#first, categorize by REVEL cutoffs as in Pejaver
merged_dbNERG_geneVar_reordered[INFO_dbNSFP_REVEL_score >= 0.932, PP3 := "Strong"]
merged_dbNERG_geneVar_reordered[INFO_dbNSFP_REVEL_score >= 0.773 & INFO_dbNSFP_REVEL_score < 0.932, PP3 := "Moderate"]
merged_dbNERG_geneVar_reordered[INFO_dbNSFP_REVEL_score >= 0.644 & INFO_dbNSFP_REVEL_score < 0.773, PP3 := "Supporting"]

#categorize by VEST cutoffs as in Pejaver for variants without REVEL scores
merged_dbNERG_geneVar_reordered[is.na(INFO_dbNSFP_REVEL_score) & INFO_dbNSFP_VEST4_score >= 0.965, PP3 := "Strong"]
merged_dbNERG_geneVar_reordered[is.na(INFO_dbNSFP_REVEL_score) & INFO_dbNSFP_VEST4_score >= 0.861 & INFO_dbNSFP_VEST4_score < 0.965, PP3 := "Moderate"]
merged_dbNERG_geneVar_reordered[is.na(INFO_dbNSFP_REVEL_score) & INFO_dbNSFP_VEST4_score >= 0.764 & INFO_dbNSFP_VEST4_score < 0.861, PP3 := "Supporting"]


#Use alphaMissense for any remaining variants
merged_dbNERG_geneVar_reordered[is.na(INFO_dbNSFP_REVEL_score) & is.na(INFO_dbNSFP_VEST4_score) & INFO_dbNSFP_AlphaMissense_score > 0.56, PP3 := "Moderate"]

#SpliceAI for Supporting or NA
merged_dbNERG_geneVar_reordered[(PP3 == "Supporting" | PP3 == "") & SpliceAI_DS_AG > 0.2 | SpliceAI_DS_AL > 0.2 | SpliceAI_DS_DG >0.2 | SpliceAI_DS_DL > 0.2, PP3 := "Moderate"]



#check PM2 - gnomad
merged_dbNERG_geneVar_reordered[INFO_dbNSFP_gnomAD_exomes_AF <= 0.01 & INFO_dbNSFP_gnomAD_genomes_AF <= 0.01, PM2 := "Supporting"]


#check PVS1
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(PVS1 = case_when(
    str_detect(Annotation, "stop_gained|frameshift_variant|splice_donor_variant|splice_acceptor_variant|start_lost") ~ "Very Strong"
  ))


#check PS2 - for now the cases need to be checked

#check PM6
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(PM6 = case_when(
    str_detect(Pedigree_Inheritance, "DN") ~ "Moderate"
  ))

#check PP1
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(PP1 = case_when(
    str_detect(Pedigree_Inheritance, "AR|AD|XL") ~ "Supporting"
  ))

#check PM4 - ptn length change
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(PM4 = case_when(
    str_detect(Annotation, "inframe") ~ "Moderate"
  ))


#check PM3 - in trans - requires more thought
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(PM3 = case_when(
    `Trans_proven?` == 1 ~ "Moderate"
  ))

#check PP4 - phenotype - all phenotypes concordant??
merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
  mutate(PP4 = "Supporting")

#check PS1
PTN.changes <- as.data.table(merged_dbNERG_geneVar_unique[,c("Feature_ID", "HGVS.c", "HGVS.p", "INFO_PtID")])

#ensure Gene_ID is properly populated for structural variants
#ensure HGVS.p is properly populated with structural variants 



#assuming all the variants are now categorized as P/LP/VUS
#####



write.table(merged_data, "merged_nERG_annotated_file_May18.txt", sep = "\t", row.names = FALSE)


#solve each presumed_genetic case with P / LP



merged_dbNERG_geneVar_unique <- merged_dbNERG_geneVar %>%
  distinct(INFO_PtID, POS, .keep_all = TRUE)

  write.table(merged_dbNERG_geneVar_reordered, "merged_nERG_annotated_file.txt", sep = "\t", row.names = FALSE)

  
  merged_dbNERG_geneVar_reordered <- fread("merged_nERG_annotated_file.txt", header = TRUE)
  
  
  merged_dbNERG_geneVar_reordered <- merged_dbNERG_geneVar_reordered %>%
    dplyr::select(Presumed_Genetic, INFO_PtID, Gene_ID, Feature_ID, HGVS.c, HGVS.p, Annotation, Annotation_Impact, INFO_dbNSFP_clinvar_clnsig, 
                  INFO_dbNSFP_clinvar_review, INFO_dbNSFP_REVEL_score, INFO_dbNSFP_AlphaMissense_score, INFO_dbNSFP_MutPred_score, INFO_dbNSFP_VEST4_score,
                  SpliceAI_DS_AG, SpliceAI_DS_AL, SpliceAI_DS_DG, SpliceAI_DS_DL,  INFO_dbNSFP_gnomAD_exomes_AF, INFO_dbNSFP_gnomAD_genomes_AF, 
                  Pedigree_Inheritance, Disease_cause_AR_AD_XL, `Trans_proven?`, 
                  Var_Interpretation, PVS1, PS1, PS2, PS3, PS4, PM1, PM2, PM3, PM4, PM5, PM6, PP1, PP2, PP3, PP4, PP5, clinvar_pm3_transseg, clinvar_pp1_diseaseseg, Phenotype, everything())
  
  
  dim(merged_dbNERG_geneVar_reordered[,202:207])
  
  merged_dbNERG_geneVar_reordered[1,-(202:207)]
