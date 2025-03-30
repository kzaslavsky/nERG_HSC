# prepare for SpliceAI
# part of nERG project
# provided as is
# Version 0.1
# Author: Kirill Zaslavsky

# Read all lines from the VCF file
all_lines <- readLines("/Users/kirillzaslavsky12/snpEff/vcf_output_fulltest_sorted_ann.vcf")

# Identify header lines (starting with ##) and data lines
header_lines <- all_lines[grepl("^##", all_lines)]
column_titles <- all_lines[grepl("^#CHROM", all_lines)]  # Column headers
data_lines <- all_lines[!grepl("^##", all_lines) & !grepl("^#CHROM", all_lines)]

# Convert data lines to a data frame
vcf_data <- fread(text = c(column_titles, data_lines), data.table = TRUE)

if ("#CHROM" %in% names(vcf_data)) {
  names(vcf_data)[which(names(vcf_data) == "#CHROM")] <- "CHROM"
}


# Filter out structural variants
VCF_structuralVariants <- vcf_data[grep("<DEL>|<DUP>", vcf_data$ALT), ]
# Filter out entries containing <DEL> or <DUP>
filtered_data <- vcf_data[-grep("<DEL>|<DUP>", ALT)]

# Save structural variants for later merging
fwrite(VCF_structuralVariants, "VCF_structuralVariants.vcf", col.names = FALSE, sep = "\t", quote = FALSE)

# Save filtered data to a new VCF file
output_file <- "VCF_for_spliceAI.vcf"
writeLines(c(header_lines, column_titles), output_file)
fwrite(filtered_data, file = output_file, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)



# Read the processed file back
processed_data <- fread("VCF_for_spliceAI_annotated.vcf", data.table = TRUE)

# get updated header
all_lines_newheader <- readLines("VCF_for_spliceAI_annotated.vcf")
header_lines <- all_lines_newheader[grepl("^##", all_lines_newheader)]


# Merge processed data with structural variants (previously filtered)
full_vcf <- rbindlist(list(processed_data, VCF_structuralVariants), use.names = TRUE, fill = TRUE)

# Identify header lines (starting with ##) and data lines

# Save the full reconstituted VCF
final_output_file <- "vcf_output_fulltest_sorted_ann_spliceAI.vcf"
writeLines(c(header_lines, column_titles), final_output_file)
fwrite(full_vcf, file = final_output_file, sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE)