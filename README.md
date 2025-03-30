# nERG_HSC
Scripts for the nERG analysis, provided as is
Source data may be made available only with a signed DTA; contact corresponding author of manuscript

nERG_Analysis_01.3.R will reproduce the analyses from the paper
-- requires nERG_function.R to run

nERG_annotateVars_11.R will parse variants in many formats and map them to hg38 to facilitate annotation
-- require nERG_function.R and nERG_genomics_function.R to run

nERG_prepForSpliceAI.R - ensures VCF file passed on SpliceAI can be handled by it
nERG_sortGeneVars_v0.4.R - sorts variants after annotation with SnpEff/SnpSift and SpliceAI

nERG_sortVCF.R - merges annotated variants with main database
nERG_CatVariants.R - attempts basic categorization for some ACMG criteria (e.g., PP3 based on REVEL cutoffs in Pejaver et al.), but does require manual expert review
