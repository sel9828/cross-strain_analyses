# Project: Microalgae Cross-strain Growth Medium Reuse Experiment
# Author: Sarah Loftus, sarah.e.loftus@gmail.com

# This code carries out differential abundance analyses on bacteria community function pathways dataset generated from picrust2.
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
#"input_data/
cts_table_pi2 <- read.delim("/Users/sarahloftus/Documents/PICRUSt/2019_cross-strain_experiment/cross_strain_all/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv", 
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

coldata <- data.frame(reps, treatment)

rownames(coldata) <- colnames(cts_C323_pi2)
rownames(coldata) <- colnames(cts_D046_pi2)
rownames(coldata) <- colnames(cts_Navi_pi2)

# DESeq dataset
dds_C323_pi2 <- DESeqDataSetFromMatrix(countData = cts_C323_pi2, 
                                   colData = coldata,
                                   design = ~ treatment)

dds_D046_pi2 <- DESeqDataSetFromMatrix(countData = cts_D046_pi2, 
                                   colData = coldata,
                                   design = ~ treatment)

dds_Navi_pi2 <- DESeqDataSetFromMatrix(countData = cts_Navi_pi2, 
                                   colData = coldata,
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
pathway_descrip <- read.csv("/Users/sarahloftus/Box Sync/Home Folder sel28/Sync/Johnson_Lab/Experiments_and_Data/Medium_Recycling/2018 Cross-strain Medium Reuse Experiment/Data_Analyses/picrust2_output_analyses/pathway_descriptions.csv")

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
write.csv(res_C323_sig_pi2, file="DESeq2_results/C323_DESeq2_results.csv")
write.csv(res_D046_sig_pi2, file="DESeq2_results/D046_DESeq2_results.csv")
write.csv(res_Navi_sig_pi2, file="DESeq2_results/Navi_DESeq2_results.csv")

