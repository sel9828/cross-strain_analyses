# Project: Microalgae Cross-strain Growth Medium Reuse Experiment
# Author: Sarah Loftus, sarah.e.loftus@gmail.com

# This code carries out differential abundance analyses on bacteria community MetaCyc pathways data generated from picrust2.
# Notes on Deseq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Input data and metadata are available on Figshare: https://doi.org/10.6084/m9.figshare.7831913.v2

# Load packages
library(DESeq2) # version 1.14.1
library(tidyverse)
library(gplots)

# Requirements for DESeq dataset:
# 1. cts = Counts matrix in absolute numbers (integers, not normalized) (function pathways as rows, samples as columns)
# 2. coldata = Table of column data; first column row names must match columns names of cts; additional columns describe the sample (e.g., treatment, replicate)
# 3. design = a column of coldata that is being used to compare the samples (e.g., treatment)

# Counts matrix (including Spirochaete OTUs)

cts_table_pi2 <- read.delim("input_data/Data6_Pathways.tsv", 
             header = T, row.names = 1)  # Pathways table

cts_pi2 <- round(cts_table_pi2[ , -c(1, 8, 15)]) # remove inoculum samples, since only comparing treatment and control on final day. Round to nearest integer

# Split cts matrix by algae 
cts_C323_pi2 <- cts_pi2[,1:6]   # Staurosira sp. C323
cts_D046_pi2 <- cts_pi2[,7:12]  # Chlorella sp. D046
cts_Navi_pi2 <- cts_pi2[,13:18] # Navicula sp. SFP

# Remove pathways that have 0 counts in all samples
cts_C323_pi2 <- cts_C323_pi2[rowSums(cts_C323_pi2)!=0, ]
cts_D046_pi2 <- cts_D046_pi2[rowSums(cts_D046_pi2)!=0, ]
cts_Navi_pi2 <- cts_Navi_pi2[rowSums(cts_Navi_pi2)!=0, ]

# Column data 
reps <- c("A", "B", "C", "D", "E", "F") # biological replicates
treatment <- c(rep("Fresh", 3), rep("Reused", 3))

coldata_C323 <- data.frame(reps, treatment)
coldata_D046 <- data.frame(reps, treatment)
coldata_Navi <- data.frame(reps, treatment)

rownames(coldata_C323) <- colnames(cts_C323_pi2)
rownames(coldata_D046) <- colnames(cts_D046_pi2)
rownames(coldata_Navi) <- colnames(cts_Navi_pi2)

# DESeq dataset
dds_C323_pi2 <- DESeqDataSetFromMatrix(countData = cts_C323_pi2, 
                                   colData = coldata_C323,
                                   design = ~ treatment)

dds_D046_pi2 <- DESeqDataSetFromMatrix(countData = cts_D046_pi2, 
                                   colData = coldata_D046,
                                   design = ~ treatment)

dds_Navi_pi2 <- DESeqDataSetFromMatrix(countData = cts_Navi_pi2, 
                                   colData = coldata_Navi,
                                   design = ~ treatment)

# Specify the control treatment
dds_C323_pi2$treatment <- relevel(dds_C323_pi2$treatment, ref = "Fresh")
dds_D046_pi2$treatment <- relevel(dds_D046_pi2$treatment, ref = "Fresh")
dds_Navi_pi2$treatment <- relevel(dds_Navi_pi2$treatment, ref = "Fresh")

# MLE
DESeq_C323_MLE_pi2 <- DESeq(dds_C323_pi2, modelMatrixType="standard", betaPrior=FALSE)
DESeq_D046_MLE_pi2 <- DESeq(dds_D046_pi2, modelMatrixType="standard", betaPrior=FALSE)
DESeq_Navi_MLE_pi2 <- DESeq(dds_Navi_pi2, modelMatrixType="standard", betaPrior=FALSE)

a <- 0.05

# MLE results 
res_C323_MLE_pi2 <- results(DESeq_C323_MLE_pi2, name="treatment_Reused_vs_Fresh", alpha = a)
res_D046_MLE_pi2 <- results(DESeq_D046_MLE_pi2, name="treatment_Reused_vs_Fresh", alpha = a)
res_Navi_MLE_pi2 <- results(DESeq_Navi_MLE_pi2, name="treatment_Reused_vs_Fresh", alpha = a)

# Summaries of adjp < a
summary(res_C323_MLE_pi2)
summary(res_D046_MLE_pi2)
summary(res_Navi_MLE_pi2)

# Tables of sig results 
# Add pathway descriptions to significant results table & just keep log2foldChange and p-adj values
pathway_descrip <- read.csv("input_data/pathway_descriptions.csv")

as.data.frame(res_C323_MLE_pi2) %>% 
  rownames_to_column("pathway") %>% 
  filter(padj < a) %>% 
  select(pathway, log2FoldChange, padj) %>% 
  inner_join(pathway_descrip, by = "pathway") -> res_C323_sig_pi2
  
as.data.frame(res_D046_MLE_pi2) %>% 
  rownames_to_column("pathway") %>% 
  filter(padj < a) %>% 
  select(pathway, log2FoldChange, padj) %>% 
  inner_join(pathway_descrip, by = "pathway") -> res_D046_sig_pi2

as.data.frame(res_Navi_MLE_pi2) %>% 
  rownames_to_column("pathway") %>% 
  filter(padj < a) %>% 
  select(pathway, log2FoldChange, padj) %>% 
  inner_join(pathway_descrip, by = "pathway") -> res_Navi_sig_pi2

# Export significant results 
write.csv(res_C323_sig_pi2, file="DESeq2_results/C323_DESeq2_picrust2_results.csv")
write.csv(res_D046_sig_pi2, file="DESeq2_results/D046_DESeq2_picrust2_results.csv")
write.csv(res_Navi_sig_pi2, file="DESeq2_results/Navi_DESeq2_picrust2_results.csv")

# Combine significant results among algae strains into one table

res_C323_sig_pi2 %>% 
  full_join(res_D046_sig_pi2, by = c("pathway", "description"), suffix = c("_C323", "_D046")) %>% 
  full_join(res_Navi_sig_pi2, by = c("pathway", "description")) %>% 
  rename(padj_Navi = padj, log2FoldChange_Navi = log2FoldChange) -> all_res_sig_pi2

# Significant pathways unique to each alga

all_res_sig_pi2 %>% 
  filter(is.na(padj_D046)) %>% 
  filter(is.na(padj_Navi)) -> C323_unique_sig_pi2

all_res_sig_pi2 %>% 
  filter(is.na(padj_C323)) %>% 
  filter(is.na(padj_Navi)) -> D046_unique_sig_pi2

all_res_sig_pi2 %>% 
  filter(is.na(padj_C323)) %>% 
  filter(is.na(padj_D046)) -> Navi_unique_sig_pi2

#### Perform same analyses on predicted metabolic pathways from OTU table that excluded the Spirochaete OTUS (#16 and 27) ####

# Counts matrix (excluding Spirochaete OTUs)

cts_table_pi2_noSpiro <- read.delim("input_data/Data7_Pathways_noSpiro.tsv", 
                            header = T, row.names = 1)  # Pathways table

cts_pi2_noSpiro <- round(cts_table_pi2_noSpiro[ , -c(1, 8, 15)]) # remove inoculum samples, since only comparing treatment and control on final day. Round to nearest integer

# Split cts matrix by algae 
cts_C323_pi2_noSpiro <- cts_pi2_noSpiro[,1:6]   # Staurosira sp. C323
cts_D046_pi2_noSpiro <- cts_pi2_noSpiro[,7:12]  # Chlorella sp. D046
cts_Navi_pi2_noSpiro <- cts_pi2_noSpiro[,13:18] # Navicula sp. SFP

# Remove pathways that have 0 counts in all samples
cts_C323_pi2_noSpiro <- cts_C323_pi2_noSpiro[rowSums(cts_C323_pi2_noSpiro)!=0, ]
cts_D046_pi2_noSpiro <- cts_D046_pi2_noSpiro[rowSums(cts_D046_pi2_noSpiro)!=0, ]
cts_Navi_pi2_noSpiro <- cts_Navi_pi2_noSpiro[rowSums(cts_Navi_pi2_noSpiro)!=0, ]

# Column data 
# reps and treatment variables are same as above.

coldata_C323_noSpiro <- data.frame(reps, treatment)
coldata_D046_noSpiro <- data.frame(reps, treatment)
coldata_Navi_noSpiro <- data.frame(reps, treatment)

rownames(coldata_C323_noSpiro) <- colnames(cts_C323_pi2_noSpiro)
rownames(coldata_D046_noSpiro) <- colnames(cts_D046_pi2_noSpiro)
rownames(coldata_Navi_noSpiro) <- colnames(cts_Navi_pi2_noSpiro)

# DESeq dataset
dds_C323_pi2_noSpiro <- DESeqDataSetFromMatrix(countData = cts_C323_pi2_noSpiro, 
                                       colData = coldata_C323_noSpiro,
                                       design = ~ treatment)

dds_D046_pi2_noSpiro <- DESeqDataSetFromMatrix(countData = cts_D046_pi2_noSpiro, 
                                       colData = coldata_D046_noSpiro,
                                       design = ~ treatment)

dds_Navi_pi2_noSpiro <- DESeqDataSetFromMatrix(countData = cts_Navi_pi2_noSpiro, 
                                       colData = coldata_Navi_noSpiro,
                                       design = ~ treatment)

# Specify the control treatment
dds_C323_pi2_noSpiro$treatment <- relevel(dds_C323_pi2_noSpiro$treatment, ref = "Fresh")
dds_D046_pi2_noSpiro$treatment <- relevel(dds_D046_pi2_noSpiro$treatment, ref = "Fresh")
dds_Navi_pi2_noSpiro$treatment <- relevel(dds_Navi_pi2_noSpiro$treatment, ref = "Fresh")

# MLE
DESeq_C323_MLE_pi2_noSpiro <- DESeq(dds_C323_pi2_noSpiro, modelMatrixType="standard", betaPrior=FALSE)
DESeq_D046_MLE_pi2_noSpiro <- DESeq(dds_D046_pi2_noSpiro, modelMatrixType="standard", betaPrior=FALSE)
DESeq_Navi_MLE_pi2_noSpiro <- DESeq(dds_Navi_pi2_noSpiro, modelMatrixType="standard", betaPrior=FALSE)

# a is same as above

# MLE results 
res_C323_MLE_pi2_noSpiro <- results(DESeq_C323_MLE_pi2_noSpiro, name="treatment_Reused_vs_Fresh", alpha = a)
res_D046_MLE_pi2_noSpiro <- results(DESeq_D046_MLE_pi2_noSpiro, name="treatment_Reused_vs_Fresh", alpha = a)
res_Navi_MLE_pi2_noSpiro <- results(DESeq_Navi_MLE_pi2_noSpiro, name="treatment_Reused_vs_Fresh", alpha = a)

# Summaries of adjp < a
summary(res_C323_MLE_pi2_noSpiro)
summary(res_D046_MLE_pi2_noSpiro)
summary(res_Navi_MLE_pi2_noSpiro)

# Tables of sig results 
# Add pathway descriptions to significant results table & just keep log2foldChange and p-adj values
# pathway_descrip is same as above

as.data.frame(res_C323_MLE_pi2_noSpiro) %>% 
  rownames_to_column("pathway") %>% 
  filter(padj < a) %>% 
  select(pathway, log2FoldChange, padj) %>% 
  inner_join(pathway_descrip, by = "pathway") -> res_C323_sig_pi2_noSpiro

as.data.frame(res_D046_MLE_pi2_noSpiro) %>% 
  rownames_to_column("pathway") %>% 
  filter(padj < a) %>% 
  select(pathway, log2FoldChange, padj) %>% 
  inner_join(pathway_descrip, by = "pathway") -> res_D046_sig_pi2_noSpiro

as.data.frame(res_Navi_MLE_pi2_noSpiro) %>% 
  rownames_to_column("pathway") %>% 
  filter(padj < a) %>% 
  select(pathway, log2FoldChange, padj) %>% 
  inner_join(pathway_descrip, by = "pathway") -> res_Navi_sig_pi2_noSpiro

# Export significant results 
write.csv(res_C323_sig_pi2_noSpiro, file="DESeq2_results/C323_DESeq2_picrust2_noSpiro_results.csv")
write.csv(res_D046_sig_pi2_noSpiro, file="DESeq2_results/D046_DESeq2_picrust2_noSpiro_results.csv")
write.csv(res_Navi_sig_pi2_noSpiro, file="DESeq2_results/Navi_DESeq2_picrust2_noSpiro_results.csv")

# Combine significant results among algae strains into one table

res_C323_sig_pi2_noSpiro %>% 
  full_join(res_D046_sig_pi2_noSpiro, by = c("pathway", "description"), suffix = c("_C323", "_D046")) %>% 
  full_join(res_Navi_sig_pi2_noSpiro, by = c("pathway", "description")) %>% 
  rename(padj_Navi = padj, log2FoldChange_Navi = log2FoldChange) -> all_res_sig_pi2_noSpiro

write.csv(all_res_sig_pi2_noSpiro, file="DESeq2_results/all_DESeq2_picrust2_noSpiro_results.csv")

# Significant pathways unique to each alga

all_res_sig_pi2_noSpiro %>% 
  filter(is.na(padj_D046)) %>% 
  filter(is.na(padj_Navi)) -> C323_unique_sig_pi2_noSpiro

all_res_sig_pi2_noSpiro %>% 
  filter(is.na(padj_C323)) %>% 
  filter(is.na(padj_Navi)) -> D046_unique_sig_pi2_noSpiro

all_res_sig_pi2_noSpiro %>% 
  filter(is.na(padj_C323)) %>% 
  filter(is.na(padj_D046)) -> Navi_unique_sig_pi2_noSpiro

