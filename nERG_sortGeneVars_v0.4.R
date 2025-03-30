## script to sort genetic variants
## part of nERG project
## provided as is
## Version 0.1 
## Author: Kirill Zaslavsky


## algorithm

#expand table per ID for every genetic variant + keeps columns necessary for variant classification alter


# Example data frame
df <- (db.nERG[,.SD, .SDcols = c("ID", "GENE", "Mutation_cDNAseq", "Mutation_PTNseq", "RefSeq_accnum", "Type",
                                 "Pedigree_Inheritance", "Trans_proven?", "Presumed_Genetic", "Phenotype")])

#split the gene and variants into unique columns
n_vars <- dim(df[, t(tstrsplit(GENE, ",", names = FALSE))])[2]
gene_var_names <- paste0("Gene", 1:n_vars)
variant_var_names <- paste0("Variant", 1:n_vars)
PTNvariant_var_names <- paste0("PTNVariant", 1:n_vars)
ACCNUMvariant_var_names <- paste0("ACCNUMVariant", 1:n_vars)
Zygosity_var_names <- paste0("Type", 1:n_vars)



genes <- df[, t(tstrsplit(GENE, ",", names = FALSE))]
vars <- df[, t(tstrsplit(Mutation_cDNAseq, ",", names = FALSE))]
PTNvars <- df[, t(tstrsplit(Mutation_PTNseq, ",", names = FALSE))]
ACCNUMvars <- df[, t(tstrsplit(RefSeq_accnum, ",", names = FALSE))]
Zygosity <- df[, t(tstrsplit(Type, ",", names = FALSE))]


df[, (gene_var_names) := (tstrsplit(GENE, ",", names = FALSE))]
df[, (variant_var_names) := (tstrsplit(Mutation_cDNAseq, ",", names = FALSE))]

df[, (PTNvariant_var_names[1:dim(PTNvars)[2]]) := (tstrsplit(Mutation_PTNseq, ",", names = FALSE))]
df[, (PTNvariant_var_names[dim(PTNvars)[2]:length(PTNvariant_var_names)]) := NA]

df[, (ACCNUMvariant_var_names[1:dim(ACCNUMvars)[2]]) := (tstrsplit(RefSeq_accnum, ",", names = FALSE))]
df[, (ACCNUMvariant_var_names[dim(ACCNUMvars)[2]:length(ACCNUMvariant_var_names)]) := NA]

df[, (Zygosity_var_names[1:dim(Zygosity)[2]]) := (tstrsplit(Type, ",", names = FALSE))]
df[, (Zygosity_var_names[dim(Zygosity)[2]:length(Zygosity_var_names)]) := NA]


# Find the unique suffixes (e.g., "1", "2", "3")
suffixes <- unique(gsub("^Gene|Variant", "", names(df)[grep("^Gene[0-9]+$|^Variant[0-9]+$", names(df))]))

# Loop through the unique suffixes
for (suffix in suffixes) {
  # Create the new column name
  new_col_name <- paste("Gene_Variant", suffix, sep = "")
  
  # Merge the corresponding columns into the new column
  df[, (new_col_name) := paste(get(paste0("Gene", suffix)), 
                               get(paste0("Variant", suffix)), 
                               get(paste0("PTNVariant", suffix)), 
                               get(paste0("ACCNUMVariant", suffix)),
                               get(paste0("Type", suffix)))]
  
  # Remove the original columns
  df[, c(paste0("Gene", suffix), paste0("Variant", suffix), paste0("PTNVariant", suffix), paste0("ACCNUMVariant", suffix), paste0("Type", suffix)) := NULL]
}


df[,c("GENE","Mutation_cDNAseq", "Mutation_PTNseq", "RefSeq_accnum", "Type") := NULL]

#melt so that each variant is on one row
melted_df <- melt(df, id.vars = grep("^Gene_Variant[0-9]+$", names(df), invert = TRUE), 
                  measure.vars = grep("^Gene_Variant[0-9]+$", names(df)), 
                  variable.name = "Gene_Variant_to_rm", value.name = "Gene_Variant")
setorder(melted_df, ID)

melted_df[,Gene_Variant_to_rm := NULL]
melted_df <- melted_df[Gene_Variant != "NA NA NA NA NA"]

names(melted_df)[1] <- "PtID"

write.table(melted_df, "melted_df.txt", row.names = FALSE, sep = "\t")


