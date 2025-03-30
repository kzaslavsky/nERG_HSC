###parse annotated VCF and merge with database
## Orivuded as is
## nERG project
## version 0.1, Mar 9 2024
## Author: Kirill Zaslavsky

#input: vcf file annotated with SnpEff / SnpSift / SpliceAI + large structural variant table

#extract annotations using info in the header for the INFO column
#break down the Annotations in ANN



# script to run
#####
# script to run

vcf_ann_file_path <- "vcf_output_fulltest_sorted_ann_spliceAI_sorted_nsfp.vcf"
#vcf_file_path <- "vcf_output_fulltest_sorted_ann_spliceAI_sorted_nsfp.vcf"

vcf_ann_parsed <- parse_vcf(vcf_ann_file_path)


#vcf_file_path <-vcf_ann_file_path

info_ann <- extract_info_annotations(vcf_ann_file_path)
info_ann_titles <- extract_ann_titles(vcf_ann_file_path)


vcf_parsed_data_interim <- vcf_ann_parsed

# Assuming vcf_parsed_data is the data frame and info_ann_titles contains the column titles
# Apply the function to parse INFO_ANN column
vcf_parsed_data_with_ann <- parse_info_ann(vcf_parsed_data_interim, info_ann_titles)

vcf_parsed_data_with_ann_filtered <- vcf_parsed_data_with_ann %>%
  dplyr::filter(str_starts(Feature_ID, "NM_")) %>%
  dplyr::filter(!str_starts(HGVS.c, pattern = "c\\.-")) %>%
  dplyr::filter(!str_starts(HGVS.c, pattern = "^c\\.\\*")) 

#now just keep unique entry per variant
vcf_parsed_data_with_ann_filtered_unique <- vcf_parsed_data_with_ann_filtered %>%
  distinct(INFO_PtID, POS, .keep_all = TRUE)


vcf_parsed_preSAI <- as.data.table(vcf_parsed_data_with_ann_filtered_unique)
titles <- extract_SpliceAI_titles(vcf_ann_file_path)
# Function to parse INFO_SpliceAI column
vcf_parsed_postSAI <- parse_info_spliceAI(vcf_parsed_preSAI, titles)



# write.table(vcf_parsed_data_with_ann_filtered_unique, file = "vcf_parsed.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#merge parsed SV table with parsed vcf table

vcf_parsed_postSAI <- as.data.table(vcf_parsed_postSAI)
setDT(vcf_parsed_postSAI)[, POS := as.character(POS)]
setDT(largeSVs)[, POS := as.character(POS)]

vcf_parsed_withLargeSVs <- merge(vcf_parsed_postSAI, largeSVs, all = TRUE)

#ensure all PtIDs that are not NA get input into INFO_PtID to prepare full merge
vcf_parsed_withLargeSVs[is.na(INFO_PtID), INFO_PtID := PtID]

#remove unnecessary columns - may need to manually inspect
names(vcf_parsed_withLargeSVs)[57]
vcf_parsed_withLargeSVs <- vcf_parsed_withLargeSVs[,-c(52:63)]



### trouble shooting
# Merge db.nERG with vcf_parsed_data_with_ann_filtered_unique
setDT(melted_df)[, PtID := as.character(PtID)]

merged_mdf_geneVar <- right_join(vcf_parsed_postSAI, melted_df, 
                                 by = c("INFO_PtID" = "PtID"))

merged_mdf_geneVar.full <- full_join(vcf_parsed_withLargeSVs, melted_df, 
                                 by = c("INFO_PtID" = "PtID"), 
                                 relationship = "many-to-many")

merged_mdf_geneVar.full <- full_join( melted_df,vcf_parsed_postSAI, 
                                     by = c("PtID" = "INFO_PtID"))


merged_mdf_geneVar_unique.full <- merged_mdf_geneVar.full %>%
  distinct(INFO_PtID, Gene_Variant, .keep_all = TRUE)


merged_mdf_geneVar_unique <- merged_mdf_geneVar.full %>%
  distinct(INFO_PtID, POS, .keep_all = TRUE)

#sanity checks
setdiff(merged_mdf_geneVar_unique[,INFO_PtID], melted_df[,PtID])
setdiff(merged_mdf_geneVar_unique[,Gene_Variant], melted_df[,Gene_Variant])



write.table(merged_mdf_geneVar_unique, "merged_troubleshoot_geneset_wSpliceAI.txt", sep = "\t", row.names = FALSE)






