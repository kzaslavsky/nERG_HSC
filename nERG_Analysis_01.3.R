# Title: Data Processing and Analysis for nERG Data
# Author: Kirill Zaslavsky
# Version: 0.13
# Date: March 29, 2025

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory to directory of file (R-studio only)
#if not using R-studio, set working directory manually
#using setwd()

#####
# Load Necessary Packages - 
###
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(sjPlot)
library(dplyr)
library(scales)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(tidyr)
library(lubridate)
library(qwraps2)
library(arsenal)
source("nERG_function.R")

dir.create("Output")
dir.create("Output/Plots")
dir.create("Output/Stats")


##### 

## Prep

#load data for remapping gene variants
db.nERG.toremap <- fread("Data/CONTRASTnERG_master_Apr29_2024_cleaned.txt", header=TRUE)
#go through pipeline to remap and recategorize
#after this, the new db


#Load full database
db.nERG <- fread("merged_data_May18_analysis.txt", header = TRUE)

#Load variant-level database
db.nERG.variantlevel <- fread("db_nERG_variantlevel.txt", header = TRUE)
db.nERG.variantlevel.toexp <- db.phenMerge[Final_Case_GeneDx != "unsolved", .SD, .SDcols = c("INFO_PtID", "Age", "SEX", "Final_Case_GeneDx","Gene_ID", "Feature_ID", "HGVS.c", "HGVS.p","Var_Interpretation", "Type","Mutation_PTNseq",
                                                                                             "Phenotype_unified3","Systemic_Issues_unified3")]
### Supp Table 5 - Genetic variants in solved cases
write.table(db.nERG.caselevel.toexp, "db_variants_foreachcase_withsolution.txt", sep = "\t", row.names = FALSE)

#load case-level database with phenotypic summary
db.nERG.caselevel.phenSum <- fread("phenlevelSumm_nERG.txt", header = TRUE) 



##### Data for Supp Tables 1-4
##extract known gene variants

#Variant-level - solved cases only
# Make a full HGVS string for each variants
db.nERG[,GeneVar := paste0(Feature_ID, ":", HGVS.c, " (", HGVS.p, ")")]

# Make an annotations data.table

dt.anno <- db.nERG[,.(Gene_ID, GeneVar, Var_Interpretation, INFO_dbNSFP_clinvar_id, clinvar_ids, INFO_dbNSFP_clinvar_clnsig, INFO_dbNSFP_clinvar_review)]


# Filter variants for solved cases only
db.nERG.varlevel.solved <- db.nERG.variantlevel[Final_Case_GeneDx == Gene_ID,] %>%
  filter(Final_Case_GeneDx == Gene_ID) %>%
  filter(Presumed_Genetic == "yes") %>% setDT %>%
  filter(Final_Case_GeneDx != "unsolved") %>% setDT

# Make a full HGVS string for each variants
db.nERG.varlevel.solved[,GeneVar := paste0(Feature_ID, ":", HGVS.c, " (", HGVS.p, ")")]


# Variants in ClinVar in solved cases
df.known <- as.data.frame(sort(table(db.nERG.varlevel.solved[!is.na(INFO_dbNSFP_clinvar_id),GeneVar])))
colnames(df.known) <- c("GeneVar", "Count")
df.known <- df.known[order(-df.known$Count), ]

df <- merge(df.known, dt.anno, by = "GeneVar", all.x = TRUE)
df <- df %>% distinct(GeneVar, .keep_all = TRUE) %>% setDT

df <- df[,.(Gene_ID, GeneVar, Count, Var_Interpretation, INFO_dbNSFP_clinvar_id, clinvar_ids, INFO_dbNSFP_clinvar_clnsig, INFO_dbNSFP_clinvar_review)]

df <- df[order(-df$Count), ]

write.table(df, file = "vars_INclinVar.txt", sep = "\t", row.names = FALSE)


# Variants not in ClinVar in solved cases
df2.novel <- data.frame(sort(table(db.nERG.varlevel.solved[is.na(INFO_dbNSFP_clinvar_id),GeneVar])))
colnames(df2.novel) <- c("GeneVar", "Count")
df2.novel <- df2.novel[order(-df2.novel$Count), ]

df2 <- merge(df2.novel, dt.anno, by = "GeneVar", all.x = TRUE)
df2 <- df2 %>% distinct(GeneVar, .keep_all = TRUE) %>% setDT

df2 <- df2[order(-df2$Count), ]


df2 <- df2[,.(Gene_ID, GeneVar, Count, Var_Interpretation, INFO_dbNSFP_clinvar_id, clinvar_ids, INFO_dbNSFP_clinvar_clnsig, INFO_dbNSFP_clinvar_review)]

write.table(df2, file = "vars_NOTINclinVar.txt", sep = "\t", row.names = FALSE)

# Variants in unsolved
df.unsolved <- as.data.frame(sort(table(db.nERG[Final_Case_GeneDx == "unsolved" & Gene_ID != "",GeneVar])))

colnames(df.unsolved) <- c("GeneVar", "Count")
df.unsolved <- df.unsolved[order(-df.known$Count), ]

df3 <- merge(df.unsolved, dt.anno, by = "GeneVar", all.x = TRUE)
df3 <- df3 %>% distinct(GeneVar, .keep_all = TRUE) %>% setDT

df3 <- df3[,.(Gene_ID, GeneVar, Count, Var_Interpretation, INFO_dbNSFP_clinvar_id, clinvar_ids, INFO_dbNSFP_clinvar_clnsig, INFO_dbNSFP_clinvar_review)]

df3 <- df3[order(-df3$Count), ]

write.table(df3, file = "vars_unsolved.txt", sep = "\t", row.names = FALSE)


##### Summary statistics on cases and variants - for Supplementary Figure 1
# Genetic pts
num_pts_genetic <- length(unique(db.nERG.caselevel[Presumed_Genetic == "yes", INFO_PtID])) #number of presumed genetic cases
num_pts <- length(unique(db.nERG.caselevel[, INFO_PtID])) #number of pts with genetic testing
num_pts_gentest.avail <- length(unique(db.nERG.caselevel[Gene_ID != "NA", INFO_PtID])) #number of pts with genetic testing

# Proportion with available tests
prop_gentest.avail <- num_pts_gentest.avail / num_pts_genetic

# Number of genes with variants reported
num_rep_genes <- length(unique(db.nERG.caselevel[,Gene_ID])) #67 if not counting NA

# Find how many potients had multiple genes in report
tableGene = (table(db.nERG.caselevel[,Gene_ID, by = INFO_PtID]))
indices = integer()
for (i in 1:dim(tableGene)[1])
{
  if (sum(tableGene[i,] > 2))
  {
    indices = c(indices, i)
  }
}

tableGene2[indices[1],] %>% filter_all(any_vars(. >0))
dim(tableGene2[1,] %>% select_if(. > 0))

tableGene2 <- as.data.frame.matrix(tableGene)
indices = integer()
for (i in 1:dim(tableGene)[1])
{
  if (dim(tableGene2[i,] %>% select_if(. > 0))[2] > 1)
  {
    indices = c(indices, i)
  }
}

for (i in 1:length(indices))
{
  print(tableGene2[indices[i],][which(tableGene2[indices[i],]>0)]  )
 # print(dim(indices[i]))
}
moreThanOneGene <- length(indices)




#####

## Figure 1

## Genetic spectrum - pie chart

#summary stats - genetic spectrum
gene.count.table <- table(db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", Final_Case_GeneDx])
gene.prop.table <- prop.table(gene.count.table)
rare_cutoff <- 0.025
rare.genes <- names(gene.prop.table[gene.prop.table < rare_cutoff])
db.nERG.caselevel.phenSum[,gene_pie := Final_Case_GeneDx]
db.nERG.caselevel.phenSum[gene_pie %in% rare.genes, gene_pie := "Rare"]



#plot for ALL
db.nERG.caselevel.phenSum[,gene_pie := Final_Case_GeneDx]
db.nERG.caselevel.phenSum[gene_pie %in% rare.genes, gene_pie := "Rare"]
gene.count.table <- table(db.nERG.caselevel.phenSum[Final_Case_GeneDx !="unsolved", gene_pie])
gene.prop.table <- as.data.frame(prop.table(gene.count.table))
gene.prop.table$Freq <- round(gene.prop.table$Freq*100, digits = 1)
setorder(gene.prop.table, -Freq)
setDT(gene.prop.table)
names(gene.prop.table)[1] <- "Gene"
gene.prop.table[,Gene := factor(Gene, levels = c("Rare", "CACNA1F", "RS1", "TRPM1", "NYX", "IDUA", "CABP4"), ordered = TRUE)]
setorder(gene.prop.table, Gene)

#set the color mapping

#Define  factor levels
factor_levels <- c("Rare", "CACNA1F", "RS1", "TRPM1", "NYX", "IDUA", "CABP4")
palette <- brewer.pal(length(factor_levels), "Set3")
color_mapping <- setNames(palette, factor_levels)

#positions
df2 <- gene.prop.table %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

gene_pie.plot.all <- ggplot(gene.prop.table, aes(x = "", y = Freq, fill = Gene)) %>% exp_theme() +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y", direction = -1) + 
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(Gene, ", ", Freq, "%")),
                   size = 4.5, nudge_x = 0.85, nudge_y=1, show.legend = FALSE) +
  labs(fill = "Gene") + scale_fill_manual(values = color_mapping) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(size = 30),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) 

ggsave("gene_pie_all.pdf", plot = gene_pie.plot.all, device = "pdf", width = 6, height = 6)


#plot for CSNB
gene.count.table <- table(db.nERG.caselevel.phenSum[Final_Case_GeneDx !="unsolved" & Phenotype_unified3 %in% c("CSNB", "CSNB with fundus changes"), gene_pie])
gene.prop.table <- as.data.frame(prop.table(gene.count.table))
gene.prop.table$Freq <- round(gene.prop.table$Freq*100, digits = 1)
#setorder(gene.prop.table, -Freq)
setDT(gene.prop.table)
names(gene.prop.table)[1] <- "Gene"
gene.prop.table[,Gene := factor(Gene, levels = c("Rare", "CACNA1F", "RS1", "TRPM1", "NYX", "IDUA", "CABP4"), ordered = TRUE)]
setorder(gene.prop.table, Gene)

# Get the positions
df2 <- gene.prop.table %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

gene_pie.plot.CSNB <- ggplot(gene.prop.table, aes(x = "", y = Freq, fill = Gene)) %>% exp_theme() +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y", direction = -1) + 
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(Gene, ", ", Freq, "%")),
                   size = 4.5, nudge_x = 0.85, nudge_y=1, show.legend = FALSE) +
  labs(fill = "Gene") + scale_fill_manual(values = color_mapping) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(size = 30),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) 

ggsave("gene_pie_CSNB.pdf", plot = gene_pie.plot.CSNB, device = "pdf", width = 6, height = 6)


#####

## CSNB Refractive Error - Fig 3


## Set up Variables:

## Group by gene & ensure only those with Age > 5 are left
groups <- db.nERG.caselevel.phenSum[!is.na(SphOU) & Age > 5, groups.gene.RefErr := .N,by = gene_pie]

# Distribution of refractive error in CSNB
dens_SphOU <- ggplot(data = db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified2 == "CSNB" & Age > 5],
                     aes(x=SphOU)) %>% exp_theme() + geom_density(alpha = 0.5, fill = "firebrick3", color = "firebrick3") +
  labs(x = "Refractive Error") +
  scale_x_continuous(breaks = seq(-20,10,5))+
  theme(legend.position = "none")

ggsave("dens_SphOU.pdf", dens_SphOU, width = 5, height= 6)

d <- density( db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified2 == "CSNB" & Age > 5, SphOU])
dens_df <- data.frame(x = d$x, y = d$y)

#fine peaks of the distribution
find_peaks <- function(dens) {
  dens %>%
    mutate(
      # Look ahead and behind within y
      y_prev = lag(y),
      y_next = lead(y)
    ) %>%
    filter(
      # Local maxima condition
      y_prev < y & y_next < y
    )
}

peaks_df <- find_peaks(dens_df)
peaks_df


# Set up variables - synaptic vs post-synaptic CSNB genes
db.nERG.caselevel.phenSum[Final_Case_GeneDx %in% c("NYX", "TRPM1", "GPR179"), CSNB.gene.syn := "CSNB_PostSynaptic_Gene"]
db.nERG.caselevel.phenSum[Final_Case_GeneDx %in% c("GNB3", "CACNA1F", "CABP4", "RDH5", "GRK1"), CSNB.gene.syn := "CSNB_NonPostSynaptic_Gene"]

# Set up groups by individual CSNB gene
db.nERG.caselevel.phenSum[Final_Case_GeneDx %in% c("CACNA1F", "NYX", "TRPM1", "GPR179", "GNB3", "CABP4", "RDH5", "GRK1"), CSNB.gene := "CSNB_Gene"]
db.nERG.caselevel.phenSum[Final_Case_GeneDx %in% c("CACNA1F", "NYX", "TRPM1"), CSNB.gene := Final_Case_GeneDx]
db.nERG.caselevel.phenSum[Final_Case_GeneDx %in% c("GPR179", "GNB3", "CABP4", "RDH5", "GRK1"), CSNB.gene := "CSNB_Others"]

db.nERG.caselevel.phenSum[, CSNB.gene := as.factor(CSNB.gene) ]
db.nERG.caselevel.phenSum[, CSNB.gene.syn := factor(CSNB.gene.syn, levels = c("CSNB_PostSynaptic_Gene", "CSNB_NonPostSynaptic_Gene"), ordered = TRUE) ]
db.nERG.caselevel.phenSum[, CSNB.gene := factor(CSNB.gene, levels = c("NYX", "TRPM1", "CACNA1F", "CSNB_Others"), ordered = TRUE) ]


# Calculate means and standard errors
summary_dt.syn <- db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified3 == "CSNB", 
                                        .(mean = mean(SphOU), se = sd(SphOU) / sqrt(.N)), by = CSNB.gene.syn]
summary_dt.csnb <- db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified3 == "CSNB", 
                                            .(mean = mean(SphOU), se = sd(SphOU) / sqrt(.N)), by = CSNB.gene]


# plot by localization
SphOU.gene.syn <- ggplot(data = db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified2 == "CSNB" & Age > 5],
       aes(x = CSNB.gene.syn, y = SphOU, fill = CSNB.gene.syn)) %>% exp_theme + geom_point(color = "black", shape = 21, size = 5, 
                                                                 position = position_jitter(width=0.15), alpha = 0.6) + 
   geom_errorbar(data = summary_dt.syn, aes(y = mean, ymin = mean - se, ymax = mean + se), 
                 width = 0.1, size = 1, color = "black") + theme(legend.position="none") +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", color = "black",width = 0.3)+
  scale_fill_manual(values = c("firebrick3", "lightskyblue4"))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(size = 30),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  labs(x = "Gene Function", y = "Refractive Error") +
  scale_y_continuous(limits = c(-20,10), breaks = seq(-20,10,5)) + 
  scale_x_discrete(labels = c("CSNB_PostSynaptic_Gene" = "Postsynaptic \n ON-Bipolar Cells", "CSNB_NonPostSynaptic_Gene" = "Other")) 

ggsave("SphOU_synaptic_gene_relabeled.pdf", plot = SphOU.gene.syn, device = "pdf", width = 5, height = 6)
  
# pt plot by gene
SphOU.gene.cnsb <- ggplot(data = db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified2 == "CSNB" & Age > 5],
       aes(x = CSNB.gene, y = SphOU, fill = CSNB.gene)) %>% exp_theme + geom_point(color = "black", shape = 21, size = 5, 
                                                                                           position = position_jitter(width=0.15), alpha = 0.4) + 
  geom_errorbar(data = summary_dt.csnb, aes(y = mean, ymin = mean - se, ymax = mean + se), 
                width = 0.2, size = 1, color = "black") + theme(legend.position="none") +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", color = "black",width = 0.3)+
  scale_fill_manual(values = c("firebrick1","firebrick4","lightskyblue2", "lightskyblue4"))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(size = 30),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  labs(x = "Gene", y = "Refractive Error") +
  scale_y_continuous(limits = c(-20,10), breaks = seq(-20,10,5)) + 
  scale_x_discrete(labels = c("NYX" = "NYX", 
                              "TRPM1" = "TRPM1",
                              "CACNA1F" = "CACNA1F",
                              "CSNB_Others" = "Others")) #+
  #coord_flip()

ggsave("SphOU_csnb_gene.pdf", plot = SphOU.gene.cnsb, device = "pdf", width = 7, height = 6)


## Linear Regression to predict impact on SpH based on gene function in CSNB - Fig 3 and Table 4
db.nERG.caselevel.phenSum[, CSNB.gene.syn := factor(CSNB.gene.syn, levels = c("CSNB_NonPostSynaptic_Gene", "CSNB_PostSynaptic_Gene"), ordered = TRUE) ]

# prediction of SpH in CSNB cases only
model.sphOU <- lm(SphOU ~ CSNB.gene.syn + Age + SEX + OU_logMAR + DA_3.0_a_uV + b_a_ratio_DA_3 +  Flicker_Peak_uV_fixed, 
                  data = db.nERG.caselevel.phenSum[!is.na(CSNB.gene) & !is.na(groups.gene.RefErr) & Phenotype_unified2 == "CSNB"& Age >=5 ,])
tab_model(model.sphOU)

#Therefore, rod ON bipolar dysfunction > myopia


#summary table 1


db.nERG.caselevel.phenSum[,Phenotype_unified3:= factor(Phenotype_unified3, levels = c("CSNB", "Photoreceptor/RPE dystrophy", "Retinoschisis", "Others"), ordered = TRUE)]
db.nERG.caselevel.phenSum[,gene_pie := factor(gene_pie, c("CACNA1F", "TRPM1", "NYX", "CABP4", "RS1", "IDUA", "Rare", "unsolved"), ordered = TRUE)]


#factor phenotypes
db.nERG.caselevel.phenSum[,Phenotype_unified3 := factor(Phenotype_unified3, 
                                                       levels = c("CSNB",
                                                                  "Retinoschisis",
                                                                  "Photoreceptor/RPE dystrophy",
                                                                  "Others"
                                                                  ), 
                                                       ordered = TRUE)]

#stratify based on presence of systemic disease
db.nERG.caselevel.phenSum[,Phenotype_unified4 := factor(Phenotype_unified4, 
                                                        levels = c("CSNB-nonMS",
                                                                   "CSNB-MS",
                                                                   "Retinoschisis","PRD-nonMS",
                                                                   "PRD-MS",
                                                                   
                                                                   "Others"
                                                        ), 
                                                        ordered = TRUE)]


tableby_obj <- tableby(Phenotype_unified3 ~ Age + SEX + OU_logMAR + SphOU + a_TYPE_unified + DA_0.01_b_uV_fixed + DA_3.0_a_uV +
                         DA_3.0_b_uV + b_a_ratio_combined + LA_3.0_a_uV_fixed + LA_3.0_a_ms_fixed +
                         LA_3.0_b_uV_fixed + Flicker_Peak_uV_fixed + Flicker_Peak_ms_fixed, data = db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved"],
                       digits = 1,
                       digits.count=0,
                       digits.pct = 1,
                       digits.n = 0,
                       digits.p = 1)

tableby_obj %>% set_labels(list(Age = "Age, yrs", OU_logMAR = "BCVA (logMAR)", 
                                SphOU = "Refractive Error (Diopters)", a_TYPE_unified = "A-wave amplitude", 
                                DA_0.01_b_uV_fixed = "DA 0.01 b-wave amplitude",
                                DA_3.0_a_uV = "DA 3.0 a-wave amplitude",
                                b_a_ratio_combined = "B:A ratio",
                                LA_3.0_a_uV_fixed = "LA 3.0 a-wave amplitude",
                                LA_3.0_a_ms_fixed = "LA 3.0 a-wave implicit time",
                                LA_3.0_b_uV_fixed = "LA 3.0 b-wave amplitude",
                                Flicker_Peak_uV_fixed = "30 Hz Flicker Peak Amplitude",
                                Flicker_Peak_ms_fixed = "30 Hz Flicker Implicit Time")) %>% write2word(file = "Table1_nERG_177.docx", digits = 1, test = FALSE)


# Table 2 - genes associations - Part 1

db.nERG.caselevel.phenSum[, gene_pie := factor(gene_pie, levels =c("CACNA1F", "TRPM1", "NYX", "CABP4", "RS1", "IDUA", "Rare"), ordered = TRUE)]

table2_genes_part1 <- tableby(Phenotype_unified4 ~ gene_pie, data = db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved"], test = FALSE,
                              digits = 1,
                              digits.count=0,
                              digits.pct = 1,
                              digits.n = 0,
                              digits.p = 1)
summary(table2_genes_part1)
write2word(table2_genes_part1, file = "Table2_part1.docx")

ctrl.nopct <- tableby.control(cat.stats=c("countN"))

# Geneorder
geneNames <- unique(db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved",Final_Case_GeneDx])
geneNames.rest <- geneNames[!(geneNames %in% c("CACNA1F", "TRPM1", "NYX", "CABP4", "RS1", "IDUA"))]
geneNames.rest <- sort(geneNames.rest)
geneNames.factor <- factor(geneNames, levels =c("CACNA1F", "TRPM1", "NYX", "CABP4", "RS1", "IDUA", geneNames.rest), ordered = TRUE)

db.nERG.caselevel.phenSum[, Final_Case_GeneDx.factor := Final_Case_GeneDx]
db.nERG.caselevel.phenSum[, Final_Case_GeneDx.factor := factor(Final_Case_GeneDx.factor, levels =c("CACNA1F", "TRPM1", "NYX", "CABP4", "RS1", "IDUA", geneNames.rest), ordered = TRUE)]


# Table 2 - genes assoc - part 2 - rare - Supp Table 6
table2_genes_part2 <- tableby(Phenotype_unified4 ~ Final_Case_GeneDx.factor, data = db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved"], test = FALSE, N = TRUE)
summary(table2_genes_part2)
write2word(table2_genes_part2, file = "Table2_part2.docx")



# Fisher exact tests for enrichment of of systemic disease

# Enrichment of systemic disease in PR dystrophies
table(db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", Sys_bin], db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", PR_bin])

tab <- table(db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", Sys_bin], db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", PR_bin])
tab

fisher.test(tab)


# Enrichment in CSNB
tab <- table(db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", Sys_bin], db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", CSNB_bin])
tab
fisher.test(tab)

# Enrichment in RS
tab <- table(db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", RS_bin], db.nERG.caselevel.phenSum[Final_Case_GeneDx != "unsolved", RS_bin])
tab
fisher.test(tab)



# Analysis for RS1 - Table 4 and Supp Table 8


#logMAR
library(lme4)
m <- glm(OU_logMAR ~ Age + SphOU + DA_3.0_a_uV + b_a_ratio_DA_3_fixed + Flicker_Peak_uV_fixed, data = db.nERG.caselevel.phenSum[Final_Case_GeneDx == "RS1" & Age > 5,])
summary(m)
tab_model(m)

m <- glm(b_a_ratio_DA_3_fixed ~ Age + SphOU + DA_3.0_a_uV + DA_3.0_b_uV + Flicker_Peak_uV_fixed, data = db.nERG.caselevel.phenSum[Final_Case_GeneDx == "RS1" &  Final_Case_GeneDx != "unsolved" & Age > 5,])
summary(m)
tab_model(m)

